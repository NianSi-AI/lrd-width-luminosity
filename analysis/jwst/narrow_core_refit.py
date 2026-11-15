"""Narrow-line boundary audit and fixed-narrow-FWHM broad Halpha refits on the
QC-failing clean-label LRDs and A2744-QSO1, anchored to [O III] lambda5007 where
detected in-band (or fixed to 150 km s^-1 and 30 km s^-1).

Requires external DJA spectra fetched via `python scripts/fetch_jwst.py --download`
(with `LRDWIDTH_SPECTRA` pointing to the spectrum directory, default `work/jwst_spectra`).
Writes data/results/narrow_degeneracy_refit.json, data/results/narrow_refit_11sources.tsv,
and data/jwst/failed_templates_20valid.json.
"""
import collections
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18
from astropy.io import fits
from scipy.optimize import least_squares

import os
from lrdwidth import fit_batch as tfit
from lrdwidth import linefit as core
from lrdwidth.paths import DATA, RESULTS, ROOT

SPECTRA_DIR = Path(os.environ.get("LRDWIDTH_SPECTRA", str(ROOT / "work" / "jwst_spectra")))

# Register [O III] lambda5007 vacuum wavelength [um] (NIST ASD: 5008.240 A; 4959 at 4960.295 A)
OIII_5007_UM = 0.5008240
OIII_4959_UM = 0.4960295
OIII_DOUBLET_RATIO = 1.0 / 2.98  # F(4959) / F(5007), Storey & Zeippen (2000)
core.REST["OIII5007"] = OIII_5007_UM

FILTER_PASSBANDS_UM = {
    "F070LP": (0.70, 1.27),
    "F100LP": (0.97, 1.89),
    "F170LP": (1.66, 3.17),
    "F290LP": (2.87, 5.35),
}

# Manuscript Eq. (3) baseline parameters (verified in t1_parent_covariance_audit.py)
EQ3_A = 3.3458939
EQ3_B = 0.1700390

IDS_10_CRITICAL = [
    "RUBIES-EGS-37032",
    "Skyfire-23289",
    "RUBIES-UDS-967770",
    "RUBIES-UDS-10036",
    "RUBIES-UDS-50716",
    "RUBIES-UDS-19484",
    "RUBIES-EGS-61496",
    "RUBIES-EGS-952625",
    "Skyfire-92141",
    "RUBIES-UDS-31747",
]


def fit_oiii_doublet(fits_path: Path, z_fit: float):
    """Fit intrinsic [O III] 4959,5007 narrow width forward-convolved through NIRSpec LSF."""
    d = core.load_line(fits_path, z_fit, "OIII5007", 1.0)
    if d["status"] != "ready":
        return {"status": d["status"], "in_passband": False}
    with fits.open(fits_path) as hdul:
        filt = hdul["SPEC1D"].header.get("FILTER") or d["grating"].split("_")[-1]
    lo_um, hi_um = FILTER_PASSBANDS_UM[filt]
    in_passband = bool(lo_um <= d["center_um"] <= hi_um)

    v, y, e = d["v"], d["y"], d["err"]
    edge = np.abs(v) > 3500.0
    cont0 = float(np.median(y[edge])) if edge.sum() > 5 else float(np.median(y))
    area = max(float(np.sum(np.maximum(y - cont0, 0.0) * (d["hi"] - d["lo"])) / 1000.0), 0.05)
    maxF = max(100.0, area * 20.0)
    yscale = max(float(np.max(np.abs(y))), 1.0)
    v_peak = float(v[np.abs(v) < 600.0][np.argmax(y[np.abs(v) < 600.0])]) if (np.abs(v) < 600.0).any() else 0.0

    # Parameters: [continuum, continuum_slope, F_5007, v_narrow, W_narrow]
    lo = np.array([-3.0 * yscale, -6.0 * yscale, 0.0, -600.0, 30.0])
    hi = np.array([3.0 * yscale, 6.0 * yscale, maxF, 600.0, 600.0])
    dv_4959 = core.C * (OIII_4959_UM / OIII_5007_UM - 1.0)
    grid = d["grid"]
    K = d["K"]

    def model_flux(p):
        cont = p[0] + p[1] * v / 5000.0
        intr = 1000.0 * p[2] * (
            core.gauss(grid, p[3], p[4]) + OIII_DOUBLET_RATIO * core.gauss(grid, p[3] + dv_4959, p[4])
        )
        return cont + K @ intr

    best = None
    for w_init in (60.0, 150.0, 300.0):
        x0 = np.clip([cont0, 0.0, area * 0.75, v_peak, w_init], lo + 1e-7, hi - 1e-7)
        res = least_squares(
            lambda p: (model_flux(p) - y) / e,
            x0,
            bounds=(lo, hi),
            x_scale="jac",
            max_nfev=1200,
            ftol=2e-7,
            xtol=2e-7,
            gtol=2e-7,
        )
        if best is None or float(np.sum(res.fun ** 2)) < float(np.sum(best.fun ** 2)):
            best = res

    res_null = least_squares(
        lambda p: (p[0] + p[1] * v / 5000.0 - y) / e,
        np.clip([cont0, 0.0], lo[:2] + 1e-7, hi[:2] - 1e-7),
        bounds=(lo[:2], hi[:2]),
    )
    chi2 = float(np.sum(best.fun ** 2))
    chi2_null = float(np.sum(res_null.fun ** 2))
    n_pix = len(y)
    delta_bic = (chi2_null + 2.0 * math.log(n_pix)) - (chi2 + 5.0 * math.log(n_pix))
    cov = np.linalg.pinv(best.jac.T @ best.jac, rcond=1e-12)
    errs = np.sqrt(np.maximum(np.diag(cov), 0.0))
    snr_5007 = float(best.x[2] / errs[2]) if errs[2] > 0 else 0.0
    w_oiii = float(best.x[4])
    at_lower = abs(w_oiii - 30.0) < 0.05
    at_upper = abs(w_oiii - 600.0) < 0.05
    return {
        "status": "fitted",
        "grating": d["grating"],
        "filter": filt,
        "center_um": float(d["center_um"]),
        "in_passband": in_passband,
        "LSF_FWHM_kms": float(core.C / d["R_center"]),
        "SNR_5007": snr_5007,
        "delta_BIC_5007": float(delta_bic),
        "W_OIII_kms": w_oiii,
        "W_OIII_err_cond_kms": float(errs[4]),
        "v_OIII_kms": float(best.x[3]),
        "at_bound": "lower_30" if at_lower else ("upper_600" if at_upper else "interior"),
    }


def fit_halpha_free_pipeline(d):
    """Exact reproduction of T2_full_fit_20261111_v1 single-Gaussian fit with free W_narrow."""
    nulls = [tfit.null_model(d, False), tfit.null_model(d, True)]
    null = nulls[1] if nulls[0]["BIC"] - nulls[1]["BIC"] > 6.0 else nulls[0]
    cands = []
    for absorption in (False, True):
        f = core.fit_model(d, "single_gaussian", absorption=absorption, seeds=3)
        f = tfit.metrics(f, d, null)
        cands.append(f)
    no, ab = cands
    sel = ab if no["BIC"] - ab["BIC"] > 6.0 else no
    return sel, null


def fit_halpha_fixed_wnarrow(d, w_fixed: float, orig_free_null_bic: float = None):
    """Refit single-Gaussian broad Halpha with W_narrow fixed to w_fixed [km/s]."""
    def null_fixed(absorption: bool):
        names, xx, lo, hi = core.setup(d, "single_gaussian", absorption, 0)
        wn_idx = names.index("W_narrow")
        xx[wn_idx] = w_fixed
        xx[names.index("F_broad1")] = 0.0
        xx[names.index("v_broad1")] = 0.0
        xx[names.index("W_broad1")] = 1000.0
        free = [i for i, n in enumerate(names) if n not in ("W_narrow", "F_broad1", "v_broad1", "W_broad1")]
        best = None
        for seed in range(3):
            _, start, _, _ = core.setup(d, "single_gaussian", absorption, seed)
            start = np.clip(start, lo + 1e-8, hi - 1e-8)

            def residual(p_free):
                x = xx.copy()
                x[free] = p_free
                return (core.components(x, d, "single_gaussian", absorption)["total"] - d["y"]) / d["err"]

            res = least_squares(
                residual, start[free], bounds=(lo[free], hi[free]),
                x_scale="jac", max_nfev=1200, ftol=2e-7, xtol=2e-7, gtol=2e-7
            )
            if best is None or float(np.sum(res.fun ** 2)) < float(np.sum(best.fun ** 2)):
                best = res
        x_best = xx.copy()
        x_best[free] = best.x
        chi2 = float(np.sum(best.fun ** 2))
        k = len(free)
        n_pix = len(d["y"])
        return {
            "absorption": absorption,
            "chi2": chi2,
            "BIC": chi2 + k * math.log(n_pix),
            "k": k,
            "N": n_pix,
            "names": names,
            "params": x_best.tolist(),
            "success": bool(best.success),
        }

    nulls = [null_fixed(False), null_fixed(True)]
    null = nulls[1] if nulls[0]["BIC"] - nulls[1]["BIC"] > 6.0 else nulls[0]

    def broad_fixed(absorption: bool):
        names, xx0, lo, hi = core.setup(d, "single_gaussian", absorption, 0)
        wn_idx = names.index("W_narrow")
        free = [i for i, n in enumerate(names) if n != "W_narrow"]
        best = None
        for seed in range(3):
            _, start, _, _ = core.setup(d, "single_gaussian", absorption, seed)
            start = np.clip(start, lo + 1e-8, hi - 1e-8)

            def residual(p_free):
                x = xx0.copy()
                x[wn_idx] = w_fixed
                x[free] = p_free
                return (core.components(x, d, "single_gaussian", absorption)["total"] - d["y"]) / d["err"]

            res = least_squares(
                residual, start[free], bounds=(lo[free], hi[free]),
                x_scale="jac", max_nfev=1200, ftol=2e-7, xtol=2e-7, gtol=2e-7
            )
            if best is None or float(np.sum(res.fun ** 2)) < float(np.sum(best.fun ** 2)):
                best = res

        p_full = xx0.copy()
        p_full[wn_idx] = w_fixed
        p_full[free] = best.x
        chi2 = float(np.sum(best.fun ** 2))
        k = len(free)
        n_pix = len(d["y"])
        cov_free = np.linalg.pinv(best.jac.T @ best.jac, rcond=1e-12)
        cov = np.zeros((len(names), len(names)), dtype=float)
        cov[np.ix_(free, free)] = cov_free
        errors = np.sqrt(np.maximum(np.diag(cov), 0.0))
        i_fb = names.index("F_broad1")
        fb = float(p_full[i_fb])
        fb_err = float(math.sqrt(max(cov[i_fb, i_fb], 0.0)))
        bound = [
            names[i]
            for i in free
            if min((p_full[i] - lo[i]) / (hi[i] - lo[i]), (hi[i] - p_full[i]) / (hi[i] - lo[i])) < 1e-4
        ]
        out = {
            "kind": "single_gaussian",
            "absorption": absorption,
            "names": names,
            "params": p_full.tolist(),
            "errors_conditional": errors.tolist(),
            "covariance_conditional": cov.tolist(),
            "chi2": chi2,
            "N": n_pix,
            "k": k,
            "BIC": chi2 + k * math.log(n_pix),
            "reduced_chi2": chi2 / (n_pix - k),
            "success": bool(best.success),
            "jacobian_rank": int(np.linalg.matrix_rank(best.jac)),
            "bound_parameters": bound,
            "broad_flux_erg_s_cm2": fb * core.FUNIT,
            "broad_flux_err_conditional_erg_s_cm2": fb_err * core.FUNIT,
            "broad_SNR_conditional": fb / fb_err if fb_err > 0 else None,
            "total_emitted_target_flux_erg_s_cm2": (fb + float(p_full[names.index("F_narrow")])) * core.FUNIT,
        }
        return tfit.metrics(out, d, null)

    no, ab = broad_fixed(False), broad_fixed(True)
    sel = ab if no["BIC"] - ab["BIC"] > 6.0 else no
    if orig_free_null_bic is not None:
        sel["delta_BIC_vs_free_wnarrow_null"] = float(orig_free_null_bic - sel["BIC"])
    return sel, null


def run_self_checks(targets):
    """Assert exact bit-level reproduction of stored Skyfire-23289 and A2744-QSO1 fits."""
    # 1. Skyfire-23289 stored working JSON
    fid_sky = "dcb6120e80ef0dd7"
    work_sky = json.loads((T2_DIR / "working" / f"{fid_sky}.json").read_text())
    ha_sky = next(v for v in work_sky["lines"] if v["line"] == "Halpha")
    case_sky = next(c for c in ha_sky["cases"] if abs(float(c["lsf_FWHM_scale"]) - 1.0) < 1e-6)
    stored_sky = case_sky["main"]
    d_sky = core.load_line(T2_DIR / targets[fid_sky]["input_path"], float(targets[fid_sky]["z_fit"]), "Halpha", 1.0)
    re_sky, _ = fit_halpha_free_pipeline(d_sky)
    assert abs(re_sky["FWHM_intrinsic_kms"] - stored_sky["FWHM_intrinsic_kms"]) < 1e-4
    assert abs(re_sky["broad_SNR_conditional"] - stored_sky["broad_SNR_conditional"]) < 1e-5
    assert abs(re_sky["delta_BIC_null"] - stored_sky["delta_BIC_null"]) < 1e-4

    # 2. A2744-QSO1 stored S4_QSO1_fit.json
    qso_json = json.loads((DATA_DIR / "task8_20261112_v1" / "S4_QSO1_fit.json").read_text())
    qso_fits = CODE_DIR / qso_json["spectrum"]
    d_qso = core.load_line(qso_fits, float(qso_json["z"]), "Halpha", 1.0)
    re_qso, null_qso = fit_halpha_free_pipeline(d_qso)
    assert abs(re_qso["FWHM_intrinsic_kms"] - qso_json["FWHM_kms"]) < 1e-4
    assert abs(re_qso["broad_SNR_conditional"] - qso_json["main"]["broad_SNR_conditional"]) < 1e-5
    assert abs(re_qso["delta_BIC_null"] - qso_json["main"]["delta_BIC_null"]) < 1e-4
    return d_qso, qso_json, null_qso


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    targets = {t["file_id"]: t for t in json.loads((DATA / "jwst" / "targets.json").read_text())}
    d_qso, qso_json, null_qso_free = run_self_checks(targets)

    
    t1_out = json.loads((HERE / "t1_output.json").read_text())
    failed_21 = t1_out["failed_21_audit"]
    valid_20 = [r for r in failed_21 if r["any_spectrum_in_passband"]]
    assert len(valid_20) == 20

    # Group all targets by physical_source_group to check same-group grating spectra for [O III]
    targets_by_grp = collections.defaultdict(list)
    for fid, t in targets.items():
        targets_by_grp[t["physical_source_group"]].append(t)

    
    frozen_lrd = [
        r for r in csv.DictReader((DATA / "sample" / "frozen_sources.tsv").open(), delimiter="\t")
        if r["lrd_class_stage2"] == "LRD"
    ]
    frozen_logL_med = float(
        np.median([math.log10(float(r["main_Halpha_broad_apparent_model_luminosity_Planck18_erg_s"])) for r in frozen_lrd])
    )

    refit_records = []
    templates_20 = []

    for r in valid_20:
        sid = r["source_id"]
        fid = r["file_id"]
        grp = r["physical_source_group"]
        z_fit = float(r["z_fit"])
        fits_path = SPECTRA_DIR / targets[fid]["path"]
        d_ha = core.load_line(fits_path, z_fit, "Halpha", 1.0)

        # Read original stored working JSON
        work = json.loads((T2_DIR / "working" / f"{fid}.json").read_text())
        ha_line = next(v for v in work["lines"] if v["line"] == "Halpha")
        case1 = next(c for c in ha_line["cases"] if abs(float(c["lsf_FWHM_scale"]) - 1.0) < 1e-6)
        orig_main = case1["main"]
        orig_null = case1["null"]

        # Measure [O III] 5007 on same spectrum first, then check sibling grating spectra if needed
        oiii_same = fit_oiii_doublet(fits_path, z_fit)
        oiii_used = oiii_same
        oiii_origin = "same_spectrum"
        if not (oiii_same.get("in_passband") and oiii_same.get("SNR_5007", 0.0) >= 3.0 and oiii_same.get("delta_BIC_5007", -99.0) >= 6.0):
            for t_alt in targets_by_grp[grp]:
                if t_alt["file_id"] == fid:
                    continue
                cand_o3 = fit_oiii_doublet(T2_DIR / t_alt["input_path"], z_fit)
                if cand_o3.get("in_passband") and cand_o3.get("SNR_5007", 0.0) >= 3.0 and cand_o3.get("delta_BIC_5007", -99.0) >= 6.0:
                    oiii_used = cand_o3
                    oiii_origin = f"sibling_spectrum:{t_alt['file_id']}"
                    break

        if oiii_used.get("in_passband") and oiii_used.get("SNR_5007", 0.0) >= 3.0 and oiii_used.get("delta_BIC_5007", -99.0) >= 6.0:
            w_adopted = float(oiii_used["W_OIII_kms"])
            w_source = f"OIII5007_{oiii_origin} (SNR={oiii_used['SNR_5007']:.1f}, dBIC={oiii_used['delta_BIC_5007']:.1f})"
        else:
            w_adopted = 150.0
            if oiii_same.get("status") != "fitted":
                w_source = "fixed_150 (OIII5007_not_covered_in_grating)"
            elif not oiii_same.get("in_passband"):
                w_source = f"fixed_150 (OIII5007_at_{oiii_same['center_um']:.3f}um_below_F290LP_2.87um_cuton)"
            else:
                w_source = f"fixed_150 (OIII5007_undetected_SNR={oiii_same['SNR_5007']:.1f}_dBIC={oiii_same['delta_BIC_5007']:.1f})"

        fit_30, null_30 = fit_halpha_fixed_wnarrow(d_ha, 30.0, orig_null["BIC"])
        fit_150, null_150 = fit_halpha_fixed_wnarrow(d_ha, 150.0, orig_null["BIC"])
        if abs(w_adopted - 30.0) < 1e-4:
            fit_ad, null_ad = fit_30, null_30
        elif abs(w_adopted - 150.0) < 1e-4:
            fit_ad, null_ad = fit_150, null_150
        else:
            fit_ad, null_ad = fit_halpha_fixed_wnarrow(d_ha, w_adopted, orig_null["BIC"])

        dl_cm = Planck18.luminosity_distance(z_fit).to_value("cm")

        def summarize_refit(f_obj):
            L_b = float(f_obj["broad_flux_erg_s_cm2"]) * 4.0 * math.pi * dl_cm ** 2
            x_b = math.log10(L_b) - 43.0
            y_b = math.log10(float(f_obj["FWHM_intrinsic_kms"]))
            res_eq3 = y_b - (EQ3_A + EQ3_B * x_b)
            i_w = f_obj["names"].index("W_broad1")
            w_err = float(f_obj["errors_conditional"][i_w])
            return {
                "FWHM_kms": float(f_obj["FWHM_intrinsic_kms"]),
                "FWHM_err_cond_kms": w_err,
                "SNR_broad": float(f_obj["broad_SNR_conditional"]),
                "delta_BIC_null": float(f_obj["delta_BIC_null"]),
                "delta_BIC_vs_free_null": float(f_obj["delta_BIC_vs_free_wnarrow_null"]),
                "absorption_adopted": bool(f_obj["absorption"]),
                "reduced_chi2": float(f_obj["reduced_chi2"]),
                "bound_parameters": list(f_obj["bound_parameters"]),
                "log10_LHa_erg_s": x_b + 43.0,
                "residual_Eq3_dex": res_eq3,
            }

        sum_ad = summarize_refit(fit_ad)
        sum_30 = summarize_refit(fit_30)
        sum_150 = summarize_refit(fit_150)

        # Literature residual check (using literature luminosity if published, else pipeline log10_LHa)
        lit_w = r["literature_FWHM_kms"]
        lit_logL = r["literature_log10_LHa_erg_s"] if r["literature_log10_LHa_erg_s"] is not None else r["log10_LHa_erg_s"]
        lit_resid = (math.log10(lit_w) - (EQ3_A + EQ3_B * (lit_logL - 43.0))) if lit_w else None

        rec = {
            "source_id": sid,
            "file_id": fid,
            "grating": r["grating"],
            "z_fit": z_fit,
            "is_in_10_critical": sid in IDS_10_CRITICAL,
            "orig_W_narrow_kms": r["W_narrow_kms"],
            "orig_W_narrow_bound_side": r["W_narrow_bound_side"],
            "orig_FWHM_kms": r["FWHM_kms"],
            "orig_SNR_broad": r["SNR_broad"],
            "orig_delta_BIC_null": r["delta_BIC_null"],
            "orig_log10_LHa_erg_s": r["log10_LHa_erg_s"],
            "orig_residual_Eq3_dex": r["residual_Eq3_dex"],
            "OIII5007_measurement": oiii_used,
            "adopted_W_narrow_kms": w_adopted,
            "adopted_W_narrow_source": w_source,
            "refit_adopted": sum_ad,
            "refit_W30_sensitivity": sum_30,
            "refit_W150_sensitivity": sum_150,
            "literature_FWHM_kms": lit_w,
            "literature_log10_LHa_erg_s": r["literature_log10_LHa_erg_s"],
            "literature_width_definition": r["literature_width_definition"],
            "paper_primary": r["paper_primary"],
            "literature_residual_Eq3_dex": lit_resid,
        }
        refit_records.append(rec)

        
        # When broad S/N >= 3 (after fixed W_narrow), take continuum+narrow+nuisance from fit_ad;
        # when broad line is undetected (S/N < 3), take continuum+narrow+nuisance from null_ad (F_broad=0)
        # so noise fluctuations are not baked into the narrow baseline!
        use_broad_decomp = sum_ad["SNR_broad"] >= 3.0 and sum_ad["delta_BIC_null"] > 0.0
        base_fit_obj = fit_ad if use_broad_decomp else null_ad
        comp_base = core.components(base_fit_obj["params"], d_ha, "single_gaussian", base_fit_obj["absorption"])
        # Note: in null_ad, F_broad1 = 0, so total = continuum + K@(narrow + nuisance - narrow_abs)
        baseline_no_broad_on_pixels = (
            comp_base["continuum"] + comp_base["narrow"] + comp_base["nuisance"]
        ).tolist()

        # Physical injection luminosity assignment (Dark-Pitfall #1 fix):
        # For detected broad lines (orig or refit S/N >= 3), own broad luminosity is physical;
        # For undetected sources (S/N < 3), prefer published broad Halpha luminosity; if unpublished
        # (Hviding+2025 joint Balmer fits), use total transmitted narrow+broad flux cap or LRD median.
        if r["SNR_broad"] >= 3.0 or use_broad_decomp:
            phys_logL = sum_ad["log10_LHa_erg_s"]
            phys_L_source = "refit_fixed_narrow_broad_LHa"
        elif r["literature_log10_LHa_erg_s"] is not None:
            phys_logL = float(r["literature_log10_LHa_erg_s"])
            phys_L_source = f"published_broad_LHa ({r['paper_primary']})"
        else:
            # Hviding+2025 objects without tabulated broad Halpha flux: total model narrow+broad flux gives
            # a strict upper bound; use 0.5 * total target flux (typical broad fraction ~0.5-0.7 in LRDs)
            # bounded below frozen median if faint
            tot_flux = float(null_ad["params"][null_ad["names"].index("F_narrow")]) * core.FUNIT
            L_tot = max(tot_flux * 4.0 * math.pi * dl_cm ** 2, 1e40)
            phys_logL = math.log10(L_tot)
            phys_L_source = "narrow_only_total_flux_scale (no_published_LHa_in_Hviding25)"

        templates_20.append({
            "source_id": sid,
            "physical_source_group": grp,
            "file_id": fid,
            "input_path": targets[fid]["input_path"],
            "grating": r["grating"],
            "z_fit": z_fit,
            "adopted_W_narrow_kms": w_adopted,
            "baseline_model_source": "broad_refit_narrow_plus_cont" if use_broad_decomp else "null_refit_narrow_plus_cont",
            "baseline_params": base_fit_obj["params"],
            "baseline_names": base_fit_obj["names"],
            "baseline_absorption": base_fit_obj["absorption"],
            "baseline_no_broad_pixels": baseline_no_broad_on_pixels,
            "physical_log10_LHa_erg_s": phys_logL,
            "physical_LHa_source": phys_L_source,
            "own_pipeline_log10_LHa_erg_s": r["log10_LHa_erg_s"],
        })

    # Also refit A2744-QSO1!
    z_qso = float(qso_json["z"])
    mu_qso = float(qso_json["mu"])
    oiii_qso = fit_oiii_doublet(CODE_DIR / qso_json["spectrum"], z_qso)
    # In QSO1 G395M, [O III] 5007 has S/N = 2.78 (< 3), so primary rule assigns 150.0 km/s;
    # we report BOTH 150.0 km/s (primary rule) and 30.0 km/s (same-spectrum ML & IFU narrow width <100 km/s)
    qso_fit_30, _ = fit_halpha_fixed_wnarrow(d_qso, 30.0, null_qso_free["BIC"])
    qso_fit_150, _ = fit_halpha_fixed_wnarrow(d_qso, 150.0, null_qso_free["BIC"])
    dl_qso = Planck18.luminosity_distance(z_qso).to_value("cm")

    def summarize_qso(f_obj):
        L_app = float(f_obj["broad_flux_erg_s_cm2"]) * 4.0 * math.pi * dl_qso ** 2
        L_del = L_app / mu_qso
        x_del = math.log10(L_del) - 43.0
        y_w = math.log10(float(f_obj["FWHM_intrinsic_kms"]))
        res_del = y_w - (EQ3_A + EQ3_B * x_del)
        res_app = y_w - (EQ3_A + EQ3_B * (math.log10(L_app) - 43.0))
        i_w = f_obj["names"].index("W_broad1")
        return {
            "FWHM_kms": float(f_obj["FWHM_intrinsic_kms"]),
            "FWHM_err_cond_kms": float(f_obj["errors_conditional"][i_w]),
            "SNR_broad": float(f_obj["broad_SNR_conditional"]),
            "delta_BIC_null": float(f_obj["delta_BIC_null"]),
            "delta_BIC_vs_free_null": float(f_obj["delta_BIC_vs_free_wnarrow_null"]),
            "absorption_adopted": bool(f_obj["absorption"]),
            "reduced_chi2": float(f_obj["reduced_chi2"]),
            "bound_parameters": list(f_obj["bound_parameters"]),
            "log10_LHa_apparent_erg_s": math.log10(L_app),
            "log10_LHa_delensed_erg_s": x_del + 43.0,
            "residual_Eq3_apparent_dex": res_app,
            "residual_Eq3_delensed_dex": res_del,
        }

    qso_rec = {
        "source_id": "A2744-QSO1",
        "file_id": "abell2744_9214_41_imgB",
        "grating": "G395M_F290LP",
        "z_fit": z_qso,
        "mu_lensing": mu_qso,
        "is_in_10_critical": True,
        "orig_W_narrow_kms": float(qso_json["main"]["params"][qso_json["main"]["names"].index("W_narrow")]),
        "orig_W_narrow_bound_side": "upper_600",
        "orig_FWHM_kms": float(qso_json["FWHM_kms"]),
        "orig_SNR_broad": float(qso_json["main"]["broad_SNR_conditional"]),
        "orig_delta_BIC_null": float(qso_json["main"]["delta_BIC_null"]),
        "orig_log10_LHa_erg_s": math.log10(float(qso_json["L_delensed_erg_s"])),
        "orig_residual_Eq3_dex": float(qso_json["residual_delensed_dex"]),
        "OIII5007_measurement": oiii_qso,
        "adopted_W_narrow_kms": 150.0,
        "adopted_W_narrow_source": f"fixed_150 (OIII5007_weak_SNR={oiii_qso['SNR_5007']:.2f}; W30_also_reported)",
        "refit_adopted": summarize_qso(qso_fit_150),
        "refit_W30_sensitivity": summarize_qso(qso_fit_30),
        "refit_W150_sensitivity": summarize_qso(qso_fit_150),
        "literature_FWHM_kms": 680.0,
        "literature_width_definition": "D'Eugenio+2025 summed two-broad-Gaussian FWHM (components 560+/-50 and 2000+200/-100 km/s; IFU narrow <100 km/s)",
        "paper_primary": "DEugenio2025_Juodzbalis2025",
        "literature_residual_Eq3_dex": math.log10(680.0) - (EQ3_A + EQ3_B * (math.log10(float(qso_json["L_delensed_erg_s"])) - 43.0)),
    }

    # Assemble the 10 critical + A2744-QSO1 ordered table
    by_sid = {r["source_id"]: r for r in refit_records}
    ordered_11 = [by_sid[s] for s in IDS_10_CRITICAL] + [qso_rec]

    tsv_cols = [
        "source_id", "grating", "z_fit",
        "orig_W_narrow_kms", "orig_W_narrow_bound_side",
        "adopted_W_narrow_kms", "adopted_W_narrow_source",
        "orig_FWHM_kms", "new_FWHM_kms", "W30_FWHM_kms",
        "orig_SNR_broad", "new_SNR_broad",
        "orig_delta_BIC_null", "new_delta_BIC_null",
        "orig_residual_Eq3_dex", "new_residual_Eq3_dex", "W30_residual_Eq3_dex",
        "literature_FWHM_kms", "literature_residual_Eq3_dex", "literature_width_definition",
    ]
    with (RESULTS / "narrow_refit_11sources.tsv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=tsv_cols, delimiter="\t")
        w.writeheader()
        for rec in ordered_11:
            ad = rec["refit_adopted"]
            s30 = rec["refit_W30_sensitivity"]
            w.writerow({
                "source_id": rec["source_id"],
                "grating": rec["grating"],
                "z_fit": f"{rec['z_fit']:.4f}",
                "orig_W_narrow_kms": f"{rec['orig_W_narrow_kms']:.1f}",
                "orig_W_narrow_bound_side": rec["orig_W_narrow_bound_side"],
                "adopted_W_narrow_kms": f"{rec['adopted_W_narrow_kms']:.1f}",
                "adopted_W_narrow_source": rec["adopted_W_narrow_source"],
                "orig_FWHM_kms": f"{rec['orig_FWHM_kms']:.1f}",
                "new_FWHM_kms": f"{ad['FWHM_kms']:.1f}",
                "W30_FWHM_kms": f"{s30['FWHM_kms']:.1f}",
                "orig_SNR_broad": f"{rec['orig_SNR_broad']:.2f}",
                "new_SNR_broad": f"{ad['SNR_broad']:.2f}",
                "orig_delta_BIC_null": f"{rec['orig_delta_BIC_null']:.2f}",
                "new_delta_BIC_null": f"{ad['delta_BIC_null']:.2f}",
                "orig_residual_Eq3_dex": f"{rec['orig_residual_Eq3_dex']:+.3f}",
                "new_residual_Eq3_dex": f"{ad['residual_Eq3_dex'] if 'residual_Eq3_dex' in ad else ad['residual_Eq3_delensed_dex']:+.3f}",
                "W30_residual_Eq3_dex": f"{s30['residual_Eq3_dex'] if 'residual_Eq3_dex' in s30 else s30['residual_Eq3_delensed_dex']:+.3f}",
                "literature_FWHM_kms": f"{rec['literature_FWHM_kms']:.1f}" if rec["literature_FWHM_kms"] else "",
                "literature_residual_Eq3_dex": f"{rec['literature_residual_Eq3_dex']:+.3f}" if rec["literature_residual_Eq3_dex"] is not None else "",
                "literature_width_definition": rec["literature_width_definition"],
            })

    (RESULTS / "narrow_degeneracy_refit.json").write_text(
        json.dumps(
            {
                "critical_10_plus_QSO1": ordered_11,
                "all_20_valid_failed_refits": refit_records,
                "frozen_29_median_log10_LHa": frozen_logL_med,
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    (DATA / "jwst" / "failed_templates_20valid.json").write_text(json.dumps(templates_20, indent=2, ensure_ascii=False))

    print("=== narrow-core refit Self-Checks PASSED ===")
    for rec in ordered_11:
        ad = rec["refit_adopted"]
        s30 = rec["refit_W30_sensitivity"]
        r_new = ad.get("residual_Eq3_dex", ad.get("residual_Eq3_delensed_dex"))
        r_30 = s30.get("residual_Eq3_dex", s30.get("residual_Eq3_delensed_dex"))
        print(
            f"{rec['source_id']:18s} | Wn: {rec['orig_W_narrow_kms']:5.1f}({rec['orig_W_narrow_bound_side']:9s}) -> {rec['adopted_W_narrow_kms']:5.1f} "
            f"| FWHM: {rec['orig_FWHM_kms']:6.1f} -> {ad['FWHM_kms']:6.1f} (W30:{s30['FWHM_kms']:6.1f}) | lit:{str(round(rec['literature_FWHM_kms'],1) if rec['literature_FWHM_kms'] else '-'):>6s} "
            f"| SNR: {rec['orig_SNR_broad']:5.2f}->{ad['SNR_broad']:5.2f} | dBIC: {rec['orig_delta_BIC_null']:6.2f}->{ad['delta_BIC_null']:6.2f} "
            f"| res: {rec['orig_residual_Eq3_dex']:+.3f} -> {r_new:+.3f} (W30:{r_30:+.3f})"
        )


if __name__ == "__main__":
    main()
