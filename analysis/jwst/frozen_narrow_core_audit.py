"""Uniform fixed-narrow-FWHM broad Halpha refits across all 44 primary-sample sources
(29 LRDs and 15 comparison broad-line AGN) plus A2744-QSO1 under three narrow-core rules:
[O III] lambda5007-anchored, universal FWHM_narrow = 150 km s^-1, and universal 30 km s^-1.

Requires external DJA spectra fetched via `python scripts/fetch_jwst.py --download`.
Writes data/results/frozen_fixed_narrow_audit.json.
"""
import collections
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18

import os
sys.path.insert(0, str(Path(__file__).resolve().parent))
import narrow_core_refit as t2
from lrdwidth import fit_batch as tfit
from lrdwidth import linefit as core
from lrdwidth import stats as S
from lrdwidth.paths import DATA, RESULTS, ROOT, WIDTHS

SPECTRA_DIR = Path(os.environ.get("LRDWIDTH_SPECTRA", str(ROOT / "work" / "jwst_spectra")))


def extract_fit_cov_row(fit_obj, z_fit: float, cls_label: str, source_id: str, mu: float = 1.0, mu_err: float = 0.0):
    """Build regression row (x, y, 2x2 cov in log10 L_Ha and log10 FWHM) from a refitted single-Gaussian."""
    names = fit_obj["names"]
    p = np.array(fit_obj["params"], dtype=float)
    cov_full = np.array(fit_obj["covariance_conditional"], dtype=float)
    ii = [names.index("F_broad1"), names.index("W_broad1")]
    sub = cov_full[np.ix_(ii, ii)]
    jac = np.diag(1.0 / (math.log(10.0) * p[ii]))
    log_cov = jac @ sub @ jac.T
    if mu_err > 0.0 and mu > 0.0:
        # L_delens = L_app / mu => Var(log10 L_delens) += (mu_err / (mu * ln 10))^2
        log_cov[0, 0] += (mu_err / (mu * math.log(10.0))) ** 2
    eig_min = float(np.linalg.eigvalsh(log_cov).min())
    assert np.isfinite(log_cov).all() and eig_min > 0.0, f"Non-SPD cov for {source_id}: {log_cov}"
    dl_cm = Planck18.luminosity_distance(z_fit).to_value("cm")
    L_erg_s = float(fit_obj["broad_flux_erg_s_cm2"]) * 4.0 * math.pi * dl_cm ** 2 / mu
    return {
        "source_id": source_id,
        "cls": cls_label,
        "x": math.log10(L_erg_s) - 43.0,
        "y": math.log10(float(fit_obj["FWHM_intrinsic_kms"])),
        "FWHM_kms": float(fit_obj["FWHM_intrinsic_kms"]),
        "SNR_broad": float(fit_obj["broad_SNR_conditional"]),
        "delta_BIC_null": float(fit_obj["delta_BIC_null"]),
        "reduced_chi2": float(fit_obj["reduced_chi2"]),
        "absorption_adopted": bool(fit_obj["absorption"]),
        "bound_parameters": list(fit_obj["bound_parameters"]),
        "cov": log_cov.tolist(),
    }


def run_pair_regression(lrd_rows, oth_rows=None):
    def fit_one(rows):
        x = np.array([r["x"] for r in rows], dtype=float)
        y = np.array([r["y"] for r in rows], dtype=float)
        C = np.array([r["cov"] for r in rows], dtype=float)
        b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
        ps = p.sum(0)
        sq = [float(v) for v in S.quant(s, ps)]
        bq = [float(v) for v in S.quant(b, p.sum(1))]
        p95 = float(S.quant(s, ps, qs=(0.95,))[0])
        co = np.polyfit(x, y, 1)
        ols_sd = float(np.std(y - np.polyval(co, x), ddof=1))
        return {
            "N": len(rows),
            "a_mean": float((p * ah).sum()),
            "slope_q": bq,
            "sigma_q": sq,
            "sigma_p95_onesided": p95,
            "raw_ols_sd": ols_sd,
            "_s": s,
            "_ps": ps,
        }

    rL = fit_one(lrd_rows)
    out = {k: v for k, v in rL.items() if not k.startswith("_")}
    if oth_rows is not None:
        rO = fit_one(oth_rows)
        dq = [float(v) for v in S.diffq(rL["_s"], rL["_ps"], rO["_ps"])]
        out["other"] = {k: v for k, v in rO.items() if not k.startswith("_")}
        out["diff_LRD_minus_other_q"] = dq
    return out


def main():
    targets = {t["file_id"]: t for t in json.loads((DATA / "jwst" / "targets.json").read_text())}
    targets_by_grp = collections.defaultdict(list)
    for fid, t in targets.items():
        targets_by_grp[t["physical_source_group"]].append(t)

    frozen_rows = list(csv.DictReader((DATA / "sample" / "frozen_sources.tsv").open(), delimiter="\t"))
    assert len(frozen_rows) == 44

    audit_44 = []
    rows_oiii_rule = []
    rows_w150_rule = []
    rows_w30_rule = []

    for r in frozen_rows:
        sid = r["source_id"]
        fid = r["file_id"]
        grp = r["physical_source_group"]
        cls_label = r["lrd_class_stage2"]
        z_fit = float(r["z_fit"])
        fits_path = SPECTRA_DIR / targets[fid]["path"]

        work = json.loads((T2_DIR / "working" / f"{fid}.json").read_text())
        ha = next(v for v in work["lines"] if v["line"] == "Halpha")
        case1 = next(c for c in ha["cases"] if abs(float(c["lsf_FWHM_scale"]) - 1.0) < 1e-6)
        orig_main = case1["main"]
        orig_null = case1["null"]
        orig_p = dict(zip(orig_main["names"], orig_main["params"]))

        # Measure [O III] 5007 on same spectrum first, then sibling grating spectra
        oiii_same = t2.fit_oiii_doublet(fits_path, z_fit)
        oiii_used = oiii_same
        oiii_origin = "same_spectrum"
        if not (oiii_same.get("in_passband") and oiii_same.get("SNR_5007", 0.0) >= 3.0 and oiii_same.get("delta_BIC_5007", -99.0) >= 6.0):
            for t_alt in targets_by_grp[grp]:
                if t_alt["file_id"] == fid:
                    continue
                cand = t2.fit_oiii_doublet(T2_DIR / t_alt["input_path"], z_fit)
                if cand.get("in_passband") and cand.get("SNR_5007", 0.0) >= 3.0 and cand.get("delta_BIC_5007", -99.0) >= 6.0:
                    oiii_used = cand
                    oiii_origin = f"sibling:{t_alt['file_id']}"
                    break

        if oiii_used.get("in_passband") and oiii_used.get("SNR_5007", 0.0) >= 3.0 and oiii_used.get("delta_BIC_5007", -99.0) >= 6.0:
            w_adopted = float(oiii_used["W_OIII_kms"])
            w_src = f"OIII_{oiii_origin}(W={w_adopted:.1f},SNR={oiii_used['SNR_5007']:.1f})"
        else:
            w_adopted = 150.0
            w_src = "fixed_150"

        d_ha = core.load_line(fits_path, z_fit, "Halpha", 1.0)
        fit_ad, _ = t2.fit_halpha_fixed_wnarrow(d_ha, w_adopted, orig_null["BIC"])
        fit_150, _ = t2.fit_halpha_fixed_wnarrow(d_ha, 150.0, orig_null["BIC"]) if abs(w_adopted - 150.0) > 1e-4 else (fit_ad, None)
        fit_30, _ = t2.fit_halpha_fixed_wnarrow(d_ha, 30.0, orig_null["BIC"]) if abs(w_adopted - 30.0) > 1e-4 else (fit_ad, None)

        row_ad = extract_fit_cov_row(fit_ad, z_fit, cls_label, sid)
        row_150 = extract_fit_cov_row(fit_150, z_fit, cls_label, sid)
        row_30 = extract_fit_cov_row(fit_30, z_fit, cls_label, sid)
        rows_oiii_rule.append(row_ad)
        rows_w150_rule.append(row_150)
        rows_w30_rule.append(row_30)

        orig_L = float(r["main_Halpha_broad_apparent_model_luminosity_Planck18_erg_s"])
        orig_w = float(r["main_FWHM_intrinsic_kms"])
        audit_44.append({
            "source_id": sid,
            "cls": cls_label,
            "grating": r["grating"],
            "z_fit": z_fit,
            "orig_W_narrow_kms": float(orig_p["W_narrow"]),
            "orig_F_narrow_over_F_broad": float(orig_p["F_narrow"] / max(orig_p["F_broad1"], 1e-12)),
            "orig_FWHM_kms": orig_w,
            "orig_log10_LHa": math.log10(orig_L),
            "OIII_W_kms": oiii_used.get("W_OIII_kms"),
            "OIII_SNR": oiii_used.get("SNR_5007"),
            "adopted_W_narrow_kms": w_adopted,
            "adopted_W_narrow_source": w_src,
            "new_FWHM_kms": row_ad["FWHM_kms"],
            "new_log10_LHa": row_ad["x"] + 43.0,
            "dlogFWHM_new_minus_orig": row_ad["y"] - math.log10(orig_w),
            "dlogLHa_new_minus_orig": (row_ad["x"] + 43.0) - math.log10(orig_L),
            "new_SNR_broad": row_ad["SNR_broad"],
            "new_bound_parameters": row_ad["bound_parameters"],
            "W150_FWHM_kms": row_150["FWHM_kms"],
            "W150_log10_LHa": row_150["x"] + 43.0,
        })

    # Also build QSO1 regression row under both w_adopted=150 and w=30
    qso_json = json.loads((DATA_DIR / "task8_20261112_v1" / "S4_QSO1_fit.json").read_text())
    d_qso = core.load_line(CODE_DIR / qso_json["spectrum"], float(qso_json["z"]), "Halpha", 1.0)
    qso_fit_150, _ = t2.fit_halpha_fixed_wnarrow(d_qso, 150.0)
    qso_fit_30, _ = t2.fit_halpha_fixed_wnarrow(d_qso, 30.0)
    qso_row_150 = extract_fit_cov_row(qso_fit_150, float(qso_json["z"]), "LRD", "A2744-QSO1", mu=7.5, mu_err=1.5)
    qso_row_30 = extract_fit_cov_row(qso_fit_30, float(qso_json["z"]), "LRD", "A2744-QSO1", mu=7.5, mu_err=1.5)

    lrd_ad = [r for r in rows_oiii_rule if r["cls"] == "LRD"]
    oth_ad = [r for r in rows_oiii_rule if r["cls"] != "LRD"]
    lrd_150 = [r for r in rows_w150_rule if r["cls"] == "LRD"]
    oth_150 = [r for r in rows_w150_rule if r["cls"] != "LRD"]
    lrd_30 = [r for r in rows_w30_rule if r["cls"] == "LRD"]
    oth_30 = [r for r in rows_w30_rule if r["cls"] != "LRD"]

    assert len(lrd_ad) == 29 and len(oth_ad) == 15
    top3_ids = {r["source_id"] for r in sorted(lrd_ad, key=lambda x: x["x"], reverse=True)[:3]}
    lrd_ad_notop3 = [r for r in lrd_ad if r["source_id"] not in top3_ids]
    assert len(lrd_ad_notop3) == 26

    reg_oiii_29 = run_pair_regression(lrd_ad, oth_ad)
    reg_oiii_26_notop3 = run_pair_regression(lrd_ad_notop3, oth_ad)
    reg_oiii_30_with_qso = run_pair_regression(lrd_ad + [qso_row_150], oth_ad)
    reg_oiii_30_with_qso30 = run_pair_regression(lrd_ad + [qso_row_30], oth_ad)
    reg_w150_29 = run_pair_regression(lrd_150, oth_150)
    reg_w30_29 = run_pair_regression(lrd_30, oth_30)

    # Hard self-checks: all regression posteriors finite and normalized
    for nm, rr in [
        ("oiii_29", reg_oiii_29),
        ("oiii_26_notop3", reg_oiii_26_notop3),
        ("oiii_30_qso150", reg_oiii_30_with_qso),
        ("w150_29", reg_w150_29),
        ("w30_29", reg_w30_29),
    ]:
        assert 0.0 < rr["sigma_q"][2] < rr["sigma_p95_onesided"] < 0.20, f"Unexpected sigma in {nm}: {rr}"

    # Compute QSO1 and 10 critical failed LRD residuals relative to the consistent fixed-W_narrow 29-LRD relation!
    a_fix = reg_oiii_29["a_mean"]
    b_fix = reg_oiii_29["slope_q"][2]
    qso_res_vs_fixed_rel_150 = qso_row_150["y"] - (a_fix + b_fix * qso_row_150["x"])
    qso_res_vs_fixed_rel_30 = qso_row_30["y"] - (a_fix + b_fix * qso_row_30["x"])

    t2_json = json.loads((RESULTS / "narrow_degeneracy_refit.json").read_text())
    crit11_vs_fixed_rel = []
    for row11 in t2_json["critical_10_plus_QSO1"]:
        rad = row11["refit_adopted"]
        logL_val = rad.get("log10_LHa_erg_s", rad.get("log10_LHa_delensed_erg_s"))
        res_eq3_val = rad.get("residual_Eq3_dex", rad.get("residual_Eq3_delensed_dex"))
        x_new = float(logL_val) - 43.0
        y_new = math.log10(float(rad["FWHM_kms"]))
        crit11_vs_fixed_rel.append({
            "source_id": row11["source_id"],
            "orig_W_narrow_kms": row11["orig_W_narrow_kms"],
            "adopted_W_narrow_kms": row11["adopted_W_narrow_kms"],
            "adopted_FWHM_kms": rad["FWHM_kms"],
            "adopted_SNR_broad": rad["SNR_broad"],
            "resid_vs_free_Eq3_dex": res_eq3_val,
            "resid_vs_fixed29_rel_dex": y_new - (a_fix + b_fix * x_new),
        })

    out = {
        "reg_frozen_OIII_anchored_29LRD_15oth": reg_oiii_29,
        "reg_frozen_OIII_anchored_dropTop3_26LRD_15oth": reg_oiii_26_notop3,
        "reg_frozen_plus_QSO1_w150_30LRD": reg_oiii_30_with_qso,
        "reg_frozen_plus_QSO1_w30_30LRD": reg_oiii_30_with_qso30,
        "reg_frozen_universal_W150_29LRD_15oth": reg_w150_29,
        "reg_frozen_universal_W30_29LRD_15oth": reg_w30_29,
        "QSO1_residual_vs_fixed_narrow_29LRD_relation": {
            "fixed_relation_a": a_fix,
            "fixed_relation_b": b_fix,
            "QSO1_W150_resid_dex": qso_res_vs_fixed_rel_150,
            "QSO1_W30_resid_dex": qso_res_vs_fixed_rel_30,
        },
        "crit11_residuals_vs_fixed_narrow_29LRD_relation": crit11_vs_fixed_rel,
        "frozen_44_audit": audit_44,
    }
    (RESULTS / "frozen_fixed_narrow_audit.json").write_text(json.dumps(out, indent=2))

    print("=== frozen fixed-narrow audit Frozen-44 Fixed-W_narrow Audit ===")
    for a in sorted([x for x in audit_44 if x["cls"] == "LRD"], key=lambda x: x["orig_log10_LHa"], reverse=True)[:8]:
        print(
            f"  {a['source_id']:18s} origLogL={a['orig_log10_LHa']:.2f}->{a['new_log10_LHa']:.2f} "
            f"| Wn: {a['orig_W_narrow_kms']:5.1f} -> {a['adopted_W_narrow_kms']:5.1f} ({a['adopted_W_narrow_source']}) "
            f"| FWHM: {a['orig_FWHM_kms']:6.1f} -> {a['new_FWHM_kms']:6.1f} (dlogW={a['dlogFWHM_new_minus_orig']:+.3f})"
        )
    print("-- Regressions --")
    for name, rr in [
        ("OIII-anchored 29 LRD", reg_oiii_29),
        ("OIII-anchored 26 LRD (no top3)", reg_oiii_26_notop3),
        ("OIII-anchored 29 + QSO1(W150)", reg_oiii_30_with_qso),
        ("OIII-anchored 29 + QSO1(W30)", reg_oiii_30_with_qso30),
        ("Universal W150 29 LRD", reg_w150_29),
        ("Universal W30 29 LRD", reg_w30_29),
    ]:
        sq = rr["sigma_q"]
        bq = rr["slope_q"]
        dq = rr["diff_LRD_minus_other_q"]
        print(
            f"  {name:32s}: a={rr['a_mean']:.4f}, b={bq[2]:.4f} [{bq[0]:.4f},{bq[4]:.4f}], "
            f"sigma_int={sq[2]:.4f} [{sq[0]:.4f},{sq[4]:.4f}] (95%<{rr['sigma_p95_onesided']:.4f}), "
            f"dSigma={dq[2]:+.4f} [{dq[0]:+.4f},{dq[4]:+.4f}]"
        )
    print(
        f"  QSO1 residual vs fixed-narrow 29-LRD relation: W150 -> {qso_res_vs_fixed_rel_150:+.4f} dex, "
        f"W30 -> {qso_res_vs_fixed_rel_30:+.4f} dex"
    )
    print("-- 11 Critical Sources vs Consistent Fixed-Narrow 29-LRD Relation --")
    for c in crit11_vs_fixed_rel:
        print(
            f"  {c['source_id']:18s} origWn={c['orig_W_narrow_kms']:5.1f} -> fixWn={c['adopted_W_narrow_kms']:5.1f} "
            f"| FWHM={c['adopted_FWHM_kms']:6.1f} | vs_Eq3={c['resid_vs_free_Eq3_dex']:+.3f} "
            f"| vs_fix29={c['resid_vs_fixed29_rel_dex']:+.3f}"
        )


if __name__ == "__main__":
    main()
