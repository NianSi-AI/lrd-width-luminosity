"""Full 49-template parent-sample injection-recovery suite across all 49 in-band
clean-label LRD spectral noise templates (29 primary templates plus 20 QC-failing
in-band parent templates), plus paired common-random-number width-offset pass-rate grids.

Requires external DJA spectra fetched via `python scripts/fetch_jwst.py --download` (directory
`LRDWIDTH_SPECTRA`, default `work/jwst_spectra`) and the dispersion tables of `scripts/fetch_calibration.sh`.
Repository inputs: data/jwst/injection_parameters.json (released main-fit parameters and spectrum/error-array
hashes of the 29 primary templates, read through templates.py), data/sample/frozen_sources.tsv (gratings),
data/jwst/failed_templates_20valid.json (written by analysis/jwst/narrow_core_refit.py), data/jwst/targets.json
(spectrum hashes), data/jwst/parent49_criteria.json, and the stored reference results named below.
Checkpoints go to --work (default work/parent49_injection). Writes data/results/parent49_injection.json.
"""
import os
for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_k] = "1"

import argparse
import collections
import concurrent.futures
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18

from lrdwidth import fit_batch as fit
from lrdwidth import linefit as c
from lrdwidth import stats as S
from lrdwidth.paths import DATA, RESULTS, ROOT, WIDTHS

SPECTRA_DIR = Path(os.environ.setdefault("LRDWIDTH_SPECTRA", str(ROOT / "work" / "jwst_spectra")))
from templates import TEMPLATES, template  # noqa: E402  (reads LRDWIDTH_SPECTRA, set above)
r1_ref = json.loads((RESULTS / "width_luminosity.json").read_text())
r2_ref = json.loads((RESULTS / "robustness.json").read_text())
# Stored 29-template sloped injection suites; the generating script is not part of this repository.
stage2_ref = json.loads((WIDTHS / "sloped_injection_recovery.json").read_text())
prereg = json.loads((DATA / "jwst" / "parent49_criteria.json").read_text())
assert prereg["frozen_before_simulation"] is True
assert prereg["parent_sample_definition"]["N_parent_valid"] == 49

OBS_SIG_MED = float(r1_ref["width_single_gaussian"]["LRD"]["sigma"][2])  # 0.0540268026 dex
OBS_OLS_SD = float(r1_ref["width_single_gaussian"]["LRD"]["ols_sd"])     # 0.0667985521 dex


def _canonical_ab(defn_key: str):
    w1_in = json.loads((WIDTHS / "regression_inputs.json").read_text())
    rows = [r for r in w1_in[defn_key] if r["cls"] == "LRD"]
    assert len(rows) == 29
    x = np.array([r["x"] for r in rows], dtype=float)
    y = np.array([r["y"] for r in rows], dtype=float)
    C = np.array([r["cov"] for r in rows], dtype=float)
    b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
    return float((p * ah).sum()), float(S.quant(b, p.sum(1))[2]), float(S.quant(s, p.sum(0))[2])


A_GAUSS, B_GAUSS, S_GAUSS = _canonical_ab("single_gaussian")
A_EXP, B_EXP, S_EXP = _canonical_ab("exponential")
assert abs(A_GAUSS - stage2_ref["canonical_relations_29LRD"]["single_gaussian"]["a_mean"]) < 1e-12
assert abs(B_GAUSS - stage2_ref["canonical_relations_29LRD"]["single_gaussian"]["b_median"]) < 1e-12
assert abs(A_EXP - stage2_ref["canonical_relations_29LRD"]["exponential"]["a_mean"]) < 1e-12
assert abs(B_EXP - stage2_ref["canonical_relations_29LRD"]["exponential"]["b_median"]) < 1e-12

# The 44 frozen sources in their original (seed-defining) order, as in injection_parameters.json, with the
# grating from frozen_sources.tsv (same order).
_FROZEN_TSV = list(csv.DictReader((DATA / "sample" / "frozen_sources.tsv").open(), delimiter="\t"))
assert [t["source"]["file_id"] for t in TEMPLATES] == [r["file_id"] for r in _FROZEN_TSV]
FROZEN_ALL = [dict(t["source"], grating=r["grating"]) for t, r in zip(TEMPLATES, _FROZEN_TSV)]
LRD_INDICES = [i for i, r in enumerate(FROZEN_ALL) if r["lrd_class_stage2"] == "LRD"]
assert len(LRD_INDICES) == 29

FAILED20_RAW = json.loads((DATA / "jwst" / "failed_templates_20valid.json").read_text())
assert len(FAILED20_RAW) == 20
TARGET_SHA256 = {t["file_id"]: t["sha256"] for t in json.loads((DATA / "jwst" / "targets.json").read_text())}


def require_spectra():
    """Stop before any fitting if one of the 49 template spectra is missing or differs from the frozen version."""
    expected = {SPECTRA_DIR / TEMPLATES[i]["input_path"]: TEMPLATES[i]["input_sha256"] for i in LRD_INDICES}
    expected.update({SPECTRA_DIR / r["input_path"]: TARGET_SHA256[r["file_id"]] for r in FAILED20_RAW})
    missing = [p for p in expected if not p.is_file()]
    hint = ("download them with `python scripts/fetch_jwst.py --download --out work/jwst_spectra` "
            "or point LRDWIDTH_SPECTRA to an existing copy")
    if missing:
        raise SystemExit(f"{len(missing)} of {len(expected)} template spectra missing (first: {missing[0]}); {hint}")
    changed = [p for p, sha in expected.items() if hashlib.sha256(p.read_bytes()).hexdigest() != sha]
    if changed:
        raise SystemExit(f"{len(changed)} template spectra differ from the frozen versions (first: {changed[0]}); {hint}")

# Worker-local cache of loaded 49 spectral templates + noiseless exponential->Gaussian calibration
_WORKER_TEMPLATES = None


def build_all_49_templates():
    """Load and cache all 49 valid parent LRD templates (0..28 frozen29, 29..48 failed20)."""
    global _WORKER_TEMPLATES
    if _WORKER_TEMPLATES is not None:
        return _WORKER_TEMPLATES

    templates = []
    # 1. 29 frozen passing LRD templates (as in the stored 29-template suites of sloped_injection_recovery.json):
    #    released LSF x1 main-fit parameters; templates.template() checks the spectrum and error-array hashes.
    for t_id, orig_idx in enumerate(LRD_INDICES):
        r = FROZEN_ALL[orig_idx]
        d, p = template(orig_idx)
        x_phys = float(r["x_logL43"])
        templates.append({
            "template_id": t_id,
            "template_group": "frozen29",
            "seed_idx": orig_idx,  # Preserves SeedSequence([20261113, sc_id, rep, orig_idx])!
            "source_id": r["source_id"],
            "grating": r["grating"],
            "z_fit": float(r["z_fit"]),
            "d": d,
            "base_pars": p,
            "x_phys": x_phys,
            "F_broad_phys": float(p["F_broad1"]),
            "x_own": x_phys,
            "F_broad_own": float(p["F_broad1"]),
        })

    # 2. 20 valid in-passband QC-failing LRD templates (Dark-Pitfall #1 & #2 clean baselines)
    for j, r20 in enumerate(FAILED20_RAW):
        t_id = 29 + j
        seed_idx = 100 + j
        z_fit = float(r20["z_fit"])
        d = c.load_line(SPECTRA_DIR / r20["input_path"], z_fit, "Halpha", 1.0)
        assert d["status"] == "ready", f"Template not ready: {r20['source_id']}"
        p_raw = dict(zip(r20["baseline_names"], r20["baseline_params"]))
        # For null_refit_narrow_plus_cont (undetected broad line), center broad line on systemic narrow Halpha
        v_b = float(p_raw["v_narrow"]) if r20["baseline_model_source"] == "null_refit_narrow_plus_cont" else float(p_raw["v_broad1"])
        dl_cm = Planck18.luminosity_distance(z_fit).to_value("cm")
        area_4pi_dl2_funit = 4.0 * math.pi * (dl_cm ** 2) * c.FUNIT

        logL_phys = float(r20["physical_log10_LHa_erg_s"])
        logL_own = float(r20["own_pipeline_log10_LHa_erg_s"])
        Fb_phys = (10.0 ** logL_phys) / area_4pi_dl2_funit
        Fb_own = max((10.0 ** logL_own) / area_4pi_dl2_funit, 1e-6)

        p_clean = {
            "continuum": float(p_raw["continuum"]),
            "continuum_slope": float(p_raw["continuum_slope"]),
            "F_narrow": float(p_raw["F_narrow"]),
            "v_narrow": float(p_raw["v_narrow"]),
            "W_narrow": float(p_raw["W_narrow"]),
            "F_broad1": Fb_phys,
            "v_broad1": v_b,
            "W_broad1": 1500.0,
            "F_NII6585": float(p_raw.get("F_NII6585", 0.0)),
        }
        templates.append({
            "template_id": t_id,
            "template_group": "failed20",
            "seed_idx": seed_idx,
            "source_id": r20["source_id"],
            "grating": r20["grating"],
            "z_fit": z_fit,
            "d": d,
            "base_pars": p_clean,
            "x_phys": logL_phys - 43.0,
            "F_broad_phys": Fb_phys,
            "x_own": logL_own - 43.0,
            "F_broad_own": Fb_own,
        })

    # Precompute 5-point noiseless exponential -> single-Gaussian width calibration per template
    exp_grid = np.linspace(2.65, 3.45, 5)
    for t in templates:
        d0 = dict(t["d"])
        pars0 = t["base_pars"].copy()
        pars0["F_broad1"] = t["F_broad_phys"]
        names_e = c.layout("Halpha", "exponential", False)
        y_g_equiv = []
        for ye in exp_grid:
            pars0["W_broad1"] = 10.0 ** float(ye)
            pa_e = [pars0[n] for n in names_e]
            d0["y"] = c.components(pa_e, d0, "exponential", False)["total"]
            f_g = c.fit_model(d0, "single_gaussian", False, seeds=1)
            y_g_equiv.append(math.log10(float(f_g["params"][f_g["names"].index("W_broad1")])))
        t["exp_calib_x"] = exp_grid
        t["exp_calib_y"] = np.array(y_g_equiv, dtype=float)

    assert len(templates) == 49
    _WORKER_TEMPLATES = templates
    return _WORKER_TEMPLATES


def _fit_and_evaluate_spectrum(d, t, rep, true_x, true_y, delta_y_inj, kind, Fb_inj):
    """Run exact T2/M4 null+broad (+/- absorption) fitting and Section 2.3 QC rules."""
    nulls = [fit.null_model(d, a) for a in (False, True)]
    null = nulls[1] if nulls[0]["BIC"] - nulls[1]["BIC"] > 6.0 else nulls[0]
    candidates = [c.fit_model(d, "single_gaussian", a, seeds=3) for a in (False, True)]
    no, ab = candidates
    f = ab if no["BIC"] - ab["BIC"] > 6.0 else no
    f = fit.metrics(f, d, null)

    bounds = set(f["bound_parameters"])
    snr_b = float(f["broad_SNR_conditional"]) if f["broad_SNR_conditional"] is not None else 0.0
    checks = dict(
        broadSNR=(f["broad_SNR_conditional"] is not None and snr_b >= 5.0),
        deltaBIC=(float(f["delta_BIC_null"]) >= 6.0),
        widthLSF=(float(f["FWHM_intrinsic_kms"]) > float(f["LSF_FWHM_kms"])),
        widthBoundary=("W_broad1" not in bounds),
        chi2=(float(f["reduced_chi2"]) < 2.0),
        converged=bool(f["success"]),
        absBoundary=not (f["absorption"] and bounds.intersection({"tau_abs", "W_abs"})),
    )

    fnames = f["names"]
    pp = np.array(f["params"], dtype=float)
    ii = [fnames.index("F_broad1"), fnames.index("W_broad1")]
    fb_fit_flux = float(f["broad_flux_erg_s_cm2"])
    w_fit_kms = float(f["FWHM_intrinsic_kms"])
    wn_fit_kms = float(pp[fnames.index("W_narrow")])
    rec_y = math.log10(w_fit_kms)

    if pp[ii[0]] > 1e-12 and pp[ii[1]] > 1e-12 and fb_fit_flux > 0.0:
        cov = np.array(f["covariance_conditional"], dtype=float)[np.ix_(ii, ii)]
        J = np.diag(1.0 / (math.log(10.0) * pp[ii]))
        logcov = J @ cov @ J.T
        cov_ok = bool(np.isfinite(logcov).all() and np.linalg.eigvalsh(logcov).min() > 0.0)
        rec_x = float(true_x) + math.log10(fb_fit_flux / (Fb_inj * c.FUNIT))
        logcov_list = logcov.tolist()
    else:
        cov_ok = False
        rec_x = None
        logcov_list = None

    selected = bool(all(checks.values()))
    if selected:
        assert cov_ok, f"Selected fit has non-PD covariance: tpl={t['template_id']}, rep={rep}"

    if kind == "single_gaussian":
        true_y_g_equiv = float(true_y)
    else:
        # Linear interpolation/extrapolation on template's noiseless exponential->Gaussian calibration
        slope_end = (t["exp_calib_y"][-1] - t["exp_calib_y"][0]) / (t["exp_calib_x"][-1] - t["exp_calib_x"][0])
        if true_y < t["exp_calib_x"][0]:
            true_y_g_equiv = float(t["exp_calib_y"][0] + slope_end * (true_y - t["exp_calib_x"][0]))
        elif true_y > t["exp_calib_x"][-1]:
            true_y_g_equiv = float(t["exp_calib_y"][-1] + slope_end * (true_y - t["exp_calib_x"][-1]))
        else:
            true_y_g_equiv = float(np.interp(true_y, t["exp_calib_x"], t["exp_calib_y"]))

    bias_equiv = rec_y - true_y_g_equiv
    resid_eq3 = (rec_y - (A_GAUSS + B_GAUSS * rec_x)) if rec_x is not None else (rec_y - (A_GAUSS + B_GAUSS * true_x))

    return {
        "replicate": int(rep),
        "template_id": int(t["template_id"]),
        "template_group": t["template_group"],
        "seed_idx": int(t["seed_idx"]),
        "source_id": t["source_id"],
        "true_x": float(true_x),
        "true_y": float(true_y),
        "true_y_gauss_equiv": float(true_y_g_equiv),
        "delta_y_inj": float(delta_y_inj),
        "selected": selected,
        "failed_checks": [k for k, ok in checks.items() if not ok],
        "recovered_SNR_broad": snr_b,
        "recovered_delta_BIC": float(f["delta_BIC_null"]),
        "recovered_W_narrow_kms": wn_fit_kms,
        "recovered_x": rec_x,
        "recovered_y": rec_y,
        "bias_vs_inj_equiv_dex": float(bias_equiv),
        "resid_vs_Eq3_dex": float(resid_eq3),
        "cov_xy_log10": logcov_list,
    }


def worker_random_draw(args):
    """Simulate one random-scatter spectrum (matching stage2 RNG SeedSequence on frozen29)."""
    sc_id, kind, a0, b0, scatter, lum_mode, rep, template_id = args
    templates = build_all_49_templates()
    t = templates[template_id]
    d = dict(t["d"])
    rng = np.random.default_rng(np.random.SeedSequence([20261113, sc_id, rep, t["seed_idx"]]))

    true_x = t["x_phys"] if lum_mode == "phys" else t["x_own"]
    Fb_inj = t["F_broad_phys"] if lum_mode == "phys" else t["F_broad_own"]

    eps_y = scatter * float(rng.normal())
    yy = a0 + b0 * true_x + eps_y
    width = float(10.0 ** yy)

    pars = t["base_pars"].copy()
    pars["F_broad1"] = Fb_inj
    pars["W_broad1"] = width
    names = c.layout("Halpha", kind, False)
    pa = [pars[n] for n in names]
    noiseless = c.components(pa, d, kind, False)["total"]
    d["y"] = noiseless + rng.normal(0.0, d["err"])

    return _fit_and_evaluate_spectrum(d, t, rep, true_x, math.log10(width), eps_y, kind, Fb_inj)


def worker_paired_offset(args):
    """Simulate paired common-random-number spectra at exact offsets delta_y in {-0.10, 0.00, +0.10}."""
    grid_seed_id, kind, a0, b0, rep, template_id, dy_offset = args
    templates = build_all_49_templates()
    t = templates[template_id]
    d = dict(t["d"])
    # Common random numbers across dy_offset: RNG seed depends ONLY on (grid_seed_id, rep, seed_idx)!
    rng = np.random.default_rng(np.random.SeedSequence([20261115, grid_seed_id, rep, t["seed_idx"]]))
    noise_vec = rng.normal(0.0, d["err"])

    true_x = t["x_phys"]
    Fb_inj = t["F_broad_phys"]
    yy = a0 + b0 * true_x + float(dy_offset)
    width = float(10.0 ** yy)

    pars = t["base_pars"].copy()
    pars["F_broad1"] = Fb_inj
    pars["W_broad1"] = width
    names = c.layout("Halpha", kind, False)
    pa = [pars[n] for n in names]
    noiseless = c.components(pa, d, kind, False)["total"]
    d["y"] = noiseless + noise_vec

    res = _fit_and_evaluate_spectrum(d, t, rep, true_x, math.log10(width), float(dy_offset), kind, Fb_inj)
    res["grid_offset_dex"] = float(dy_offset)
    return res


def fit_replicate_posterior(sel_rows):
    """Bayesian 2D-covariance posterior (lrdstats.post2) and OLS residual SD on passing rows."""
    x = np.array([r["recovered_x"] for r in sel_rows], dtype=float)
    y = np.array([r["recovered_y"] for r in sel_rows], dtype=float)
    C = np.array([r["cov_xy_log10"] for r in sel_rows], dtype=float)
    b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
    ps = p.sum(0)
    sig_q = [float(v) for v in S.quant(s, ps)]
    sig_p95 = float(S.quant(s, ps, qs=(0.95,))[0])
    b_q = [float(v) for v in S.quant(b, p.sum(1))]
    a_hat = float((p * ah).sum())
    co = np.polyfit(x, y, 1)
    ols_sd = float(np.std(y - np.polyval(co, x), ddof=1))
    return {
        "N_selected": len(sel_rows),
        "sigma_q": sig_q,
        "sigma_p95_onesided": sig_p95,
        "slope_q": b_q,
        "a_hat": a_hat,
        "ols_sd": ols_sd,
        "ols_slope": float(co[0]),
    }


def summarize_suite(name, kind, a0, b0, scatter, rows, n_reps=100):
    """Summarize both the 49-parent sample and the 29-frozen subset, plus failed-fit S/N bias."""
    by_rep_all = collections.defaultdict(list)
    by_rep_29 = collections.defaultdict(list)
    rej_all = collections.Counter()
    rej_29 = collections.Counter()
    pass_by_tpl = collections.Counter()
    tot_by_tpl = collections.Counter()

    for r in rows:
        by_rep_all[r["replicate"]].append(r)
        tot_by_tpl[r["template_id"]] += 1
        if r["selected"]:
            pass_by_tpl[r["template_id"]] += 1
        for fc in r["failed_checks"]:
            rej_all[fc] += 1
        if r["template_group"] == "frozen29":
            by_rep_29[r["replicate"]].append(r)
            for fc in r["failed_checks"]:
                rej_29[fc] += 1

    assert len(by_rep_all) == n_reps
    reps_all, reps_29 = [], []
    for rep in range(n_reps):
        rr_all = by_rep_all[rep]
        rr_29 = by_rep_29[rep]
        assert len(rr_all) == 49 and len(rr_29) == 29
        sel_all = [x for x in rr_all if x["selected"]]
        sel_29 = [x for x in rr_29 if x["selected"]]
        assert len(sel_all) >= 4 and len(sel_29) >= 4
        reps_all.append(fit_replicate_posterior(sel_all))
        reps_29.append(fit_replicate_posterior(sel_29))

    def _pack_reps(reps, n_templates, rej_cnt):
        ols_vals = np.array([r["ols_sd"] for r in reps])
        sig_meds = np.array([r["sigma_q"][2] for r in reps])
        sig_p95s = np.array([r["sigma_p95_onesided"] for r in reps])
        b_meds = np.array([r["slope_q"][2] for r in reps])
        a_hats = np.array([r["a_hat"] for r in reps])
        ns = np.array([r["N_selected"] for r in reps])
        return {
            "N_replicates": n_reps,
            "total_attempts": n_reps * n_templates,
            "total_selected": int(ns.sum()),
            "N_selected_min_p16_p50_p84_max": [float(v) for v in np.quantile(ns, [0.0, 0.16, 0.50, 0.84, 1.0])],
            "quality_rejections": dict(rej_cnt),
            "bayes_sigma_med_quantiles": [float(v) for v in np.quantile(sig_meds, [0.025, 0.16, 0.50, 0.84, 0.975])],
            "bayes_sigma_p95_median": float(np.median(sig_p95s)),
            "bayes_hits_le_obs_Fsig": int((sig_meds <= OBS_SIG_MED).sum()),
            "bayes_hits_le_00540301": int((sig_meds <= 0.0540301).sum()),
            "ols_sd_quantiles": [float(v) for v in np.quantile(ols_vals, [0.025, 0.16, 0.50, 0.84, 0.975])],
            "ols_hits_le_obs_Frms": int((ols_vals <= OBS_OLS_SD).sum()),
            "bayes_slope_med_quantiles": [float(v) for v in np.quantile(b_meds, [0.025, 0.16, 0.50, 0.84, 0.975])],
            "recovered_a_hat_p16_p50_p84": [float(v) for v in np.quantile(a_hats, [0.16, 0.50, 0.84])],
        }

    # S/N-binned width bias among QC-failing fits (both all 49 templates and failed20 templates specifically)
    snr_bins = [
        ("SNR_lt_2", 0.0, 2.0),
        ("SNR_2_to_3", 2.0, 3.0),
        ("SNR_3_to_4", 3.0, 4.0),
        ("SNR_4_to_7", 4.0, 7.0),
        ("SNR_ge_7", 7.0, 1e9),
    ]

    def _bin_failed_bias(subset_rows):
        out_bins = {}
        failed_rows = [r for r in subset_rows if not r["selected"]]
        for bname, lo, hi in snr_bins:
            sub = [r for r in failed_rows if lo <= r["recovered_SNR_broad"] < hi]
            if not sub:
                out_bins[bname] = {"N": 0}
                continue
            b_eq = np.array([r["bias_vs_inj_equiv_dex"] for r in sub], dtype=float)
            b_rel = np.array([r["resid_vs_Eq3_dex"] for r in sub], dtype=float)
            wn600 = sum(1 for r in sub if r["recovered_W_narrow_kms"] >= 599.0)
            out_bins[bname] = {
                "N": len(sub),
                "bias_vs_inj_equiv_p16_p50_p84": [float(v) for v in np.quantile(b_eq, [0.16, 0.50, 0.84])],
                "frac_bias_vs_inj_equiv_gt_0p09": float((b_eq > 0.09).mean()),
                "resid_vs_Eq3_p16_p50_p84": [float(v) for v in np.quantile(b_rel, [0.16, 0.50, 0.84])],
                "frac_resid_vs_Eq3_gt_0p09": float((b_rel > 0.09).mean()),
                "frac_W_narrow_hit_600": float(wn600 / len(sub)),
            }
        return out_bins

    # Binned pass rate vs injected offset from continuous draws
    offset_bins_cont = {}
    for oname, olo, ohi in [("minus_0p10_bin", -0.15, -0.05), ("zero_0p00_bin", -0.05, 0.05), ("plus_0p10_bin", 0.05, 0.15)]:
        sub_all = [r for r in rows if olo <= r["delta_y_inj"] < ohi]
        sub_f20 = [r for r in sub_all if r["template_group"] == "failed20"]
        offset_bins_cont[oname] = {
            "N_all49": len(sub_all),
            "pass_rate_all49": float(np.mean([r["selected"] for r in sub_all])) if sub_all else None,
            "N_failed20": len(sub_f20),
            "pass_rate_failed20": float(np.mean([r["selected"] for r in sub_f20])) if sub_f20 else None,
        }

    templates = build_all_49_templates()
    per_template_pass = []
    for t in templates:
        tid = t["template_id"]
        per_template_pass.append({
            "template_id": tid,
            "source_id": t["source_id"],
            "template_group": t["template_group"],
            "grating": t["grating"],
            "x_logL43": t["x_phys"],
            "pass_rate": float(pass_by_tpl[tid] / max(tot_by_tpl[tid], 1)),
        })

    return {
        "name": name,
        "profile": kind,
        "injected_a0": a0,
        "injected_b0": b0,
        "injected_sigma_parent_dex": scatter,
        "parent49": _pack_reps(reps_all, 49, rej_all),
        "subset_frozen29": _pack_reps(reps_29, 29, rej_29),
        "failed_fits_SNR_bins_all49": _bin_failed_bias(rows),
        "failed_fits_SNR_bins_failed20_only": _bin_failed_bias([r for r in rows if r["template_group"] == "failed20"]),
        "continuous_offset_pass_rates": offset_bins_cont,
        "per_template_pass_rates": per_template_pass,
    }


def summarize_paired_offset_grid(rows_grid, n_reps_grid):
    """Summarize exact common-random-number pass rates at delta_y in {-0.10, 0.00, +0.10} dex."""
    out = {}
    templates = build_all_49_templates()
    for dy in (-0.10, 0.00, 0.10):
        key = f"dy_{dy:+.2f}"
        sub = [r for r in rows_grid if abs(r["grid_offset_dex"] - dy) < 1e-6]
        assert len(sub) == 49 * n_reps_grid
        sub29 = [r for r in sub if r["template_group"] == "frozen29"]
        sub20 = [r for r in sub if r["template_group"] == "failed20"]
        by_tpl = collections.defaultdict(list)
        for r in sub:
            by_tpl[r["template_id"]].append(r["selected"])
        out[key] = {
            "pass_rate_all49": float(np.mean([r["selected"] for r in sub])),
            "pass_rate_frozen29": float(np.mean([r["selected"] for r in sub29])),
            "pass_rate_failed20": float(np.mean([r["selected"] for r in sub20])),
            "per_template_pass_rate": {
                templates[tid]["source_id"]: float(np.mean(by_tpl[tid])) for tid in range(49)
            },
        }
    out["ratio_plus010_over_zero_all49"] = out["dy_+0.10"]["pass_rate_all49"] / max(out["dy_+0.00"]["pass_rate_all49"], 1e-9)
    out["ratio_plus010_over_zero_failed20"] = out["dy_+0.10"]["pass_rate_failed20"] / max(out["dy_+0.00"]["pass_rate_failed20"], 1e-9)
    return out


def run_t3b_selection_corrected_posterior(offset_grid_summary):
    """selection-corrected upper limit (executed only if C1 truncation criterion triggers): selection-weighted 29-LRD posterior."""
    w1_in = json.loads((WIDTHS / "regression_inputs.json").read_text())
    rows29 = [r for r in w1_in["single_gaussian"] if r["cls"] == "LRD"]
    assert len(rows29) == 29
    x = np.array([r["x"] for r in rows29], dtype=float)
    y = np.array([r["y"] for r in rows29], dtype=float)
    C = np.array([r["cov"] for r in rows29], dtype=float)

    # Unweighted baseline posterior on (b, sigma_int) grid
    b_grid, s_grid, p_base, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])

    # Build smooth empirical selection curve S_i(dy) for each of the 29 passing objects by piecewise-linear
    # log-odds interpolation across dy in {-0.10, 0.00, +0.10}
    dy_nodes = np.array([-0.10, 0.00, 0.10])
    gh_z, gh_w = np.polynomial.hermite.hermgauss(25)
    gh_z = gh_z * math.sqrt(2.0)
    gh_w = gh_w / math.sqrt(math.pi)

    p_mat = np.zeros((29, 3), dtype=float)
    for i, r in enumerate(rows29):
        sid = r["id"] if "id" in r else FROZEN_ALL[LRD_INDICES[i]]["source_id"]
        for k, dkey in enumerate(("dy_-0.10", "dy_+0.00", "dy_+0.10")):
            p_mat[i, k] = np.clip(offset_grid_summary[dkey]["per_template_pass_rate"][sid], 0.02, 0.995)
    logit_mat = np.log(p_mat / (1.0 - p_mat))
    slope_hi = (logit_mat[:, 2] - logit_mat[:, 1]) / 0.10
    slope_lo = (logit_mat[:, 1] - logit_mat[:, 0]) / 0.10

    # Marginal selection normalization Z(sigma_int) = prod_i E_{dy ~ N(0, sigma_int^2)}[S_i(dy)]
    log_Z_s = np.zeros_like(s_grid)
    for js, sig in enumerate(s_grid):
        dy_pts = sig * gh_z[None, :]  # (1, 25)
        lg = np.where(
            dy_pts >= 0.0,
            logit_mat[:, 1:2] + slope_hi[:, None] * dy_pts,
            logit_mat[:, 1:2] + slope_lo[:, None] * dy_pts,
        )
        s_prob = 1.0 / (1.0 + np.exp(-np.clip(lg, -20.0, 20.0)))
        marg_pass_i = np.sum(s_prob * gh_w[None, :], axis=1)
        log_Z_s[js] = float(np.sum(np.log(np.clip(marg_pass_i, 1e-12, 1.0))))

    # Truncation-corrected posterior divides likelihood by prod_i P(pass_i | sigma_int)
    log_w = -log_Z_s
    log_w -= np.max(log_w)
    w_s = np.exp(log_w)
    p_corr = p_base * w_s[None, :]
    p_corr /= p_corr.sum()
    ps_corr = p_corr.sum(0)
    return {
        "sigma_q": [float(v) for v in S.quant(s_grid, ps_corr)],
        "sigma_p95_onesided": float(S.quant(s_grid, ps_corr, qs=(0.95,))[0]),
        "slope_q": [float(v) for v in S.quant(b_grid, p_corr.sum(1))],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke-test", action="store_true", help="Run 1 replicate of 49 templates to verify pipeline.")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--work", type=Path, default=ROOT / "work" / "parent49_injection", help="checkpoint directory")
    args = parser.parse_args()
    require_spectra()

    t_start = time.monotonic()
    # Warm up and verify all 49 templates + noiseless calibration in main process first
    t0_cal = time.monotonic()
    tpls = build_all_49_templates()
    print(f"[Init] Loaded {len(tpls)} valid parent templates + noiseless calibrations in {time.monotonic() - t0_cal:.2f}s", flush=True)

    if args.smoke_test:
        t0 = time.monotonic()
        tasks = [(11, "single_gaussian", A_GAUSS, B_GAUSS, 0.10, "phys", 0, tid) for tid in range(49)]
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
            rows = list(pool.map(worker_random_draw, tasks, chunksize=4))
        assert len(rows) == 49
        n_sel_49 = sum(r["selected"] for r in rows)
        n_sel_29 = sum(r["selected"] for r in rows if r["template_group"] == "frozen29")
        print(f"[Smoke Test PASSED] 49 fits in {time.monotonic() - t0:.2f}s | selected: {n_sel_49}/49 (frozen29: {n_sel_29}/29)", flush=True)
        return

    ckpt_dir = args.work / "checkpoints"
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    primary_scenarios = [
        (11, "gauss_rel_s010", "single_gaussian", A_GAUSS, B_GAUSS, 0.10),
        (13, "exp_matched_s010", "exponential", A_EXP, B_EXP, 0.10),
        (20, "gauss_rel_s005", "single_gaussian", A_GAUSS, B_GAUSS, 0.05),
        (21, "exp_matched_s005", "exponential", A_EXP, B_EXP, 0.05),
    ]

    suite_summaries = {}
    cached_rows_s010 = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        # 1. Run 4 primary 49-template x 100-trial scenarios
        for sc_id, name, kind, a0, b0, scatter in primary_scenarios:
            ckpt_path = ckpt_dir / f"{name}_parent49.json"
            if ckpt_path.exists():
                saved = json.loads(ckpt_path.read_text())
                rows = saved["rows"]
                summary = saved["summary"]
                print(f"[Cached] {name} ({len(rows)} rows)", flush=True)
            else:
                t0 = time.monotonic()
                tasks = [
                    (sc_id, kind, a0, b0, scatter, "phys", rep, tid)
                    for rep in range(100)
                    for tid in range(49)
                ]
                rows = []
                for idx_done, res in enumerate(pool.map(worker_random_draw, tasks, chunksize=8), 1):
                    rows.append(res)
                    if idx_done % 980 == 0 or idx_done == len(tasks):
                        print(f"[{name}] {idx_done}/{len(tasks)} fits ({time.monotonic() - t0:.1f}s)", flush=True)
                summary = summarize_suite(name, kind, a0, b0, scatter, rows, n_reps=100)
                tmp = ckpt_path.with_suffix(".tmp")
                tmp.write_text(json.dumps({"summary": summary, "rows": rows}, ensure_ascii=False))
                tmp.rename(ckpt_path)

            # HARD SELF-CHECK against data/widths/sloped_injection_recovery.json on the 29-frozen subset for sc_id 11 & 13!
            if name in ("gauss_rel_s010", "exp_matched_s010"):
                ref_s2 = stage2_ref["sloped_simulations"][name]
                sub29 = summary["subset_frozen29"]
                assert sub29["total_selected"] == ref_s2["total_selected"], (
                    f"Stage2 mismatch in {name} total_selected: {sub29['total_selected']} vs {ref_s2['total_selected']}"
                )
                assert sub29["bayes_hits_le_obs_Fsig"] == ref_s2["bayes_hits_le_obs_Fsig"], (
                    f"Stage2 mismatch in {name} hits: {sub29['bayes_hits_le_obs_Fsig']} vs {ref_s2['bayes_hits_le_obs_Fsig']}"
                )
                assert abs(sub29["bayes_sigma_med_quantiles"][2] - ref_s2["bayes_sigma_med_quantiles"][2]) < 1e-10, (
                    f"Stage2 mismatch in {name} median sigma_int: {sub29['bayes_sigma_med_quantiles'][2]} vs {ref_s2['bayes_sigma_med_quantiles'][2]}"
                )
                print(f"  [ASSERT PASSED] {name} 29-frozen subset matches sloped_injection_recovery.json to < 1e-10!", flush=True)
                cached_rows_s010[name] = rows

            suite_summaries[name] = summary
            p49 = summary["parent49"]
            f47 = summary["failed_fits_SNR_bins_all49"]["SNR_4_to_7"]
            print(
                f"[DONE {name:18s}] Parent49 N_sel_med={p49['N_selected_min_p16_p50_p84_max'][2]:4.1f} | "
                f"sigma_int med={p49['bayes_sigma_med_quantiles'][2]:.4f} "
                f"(hits<=0.05403: {p49['bayes_hits_le_00540301']:2d}/100) | "
                f"b_med={p49['bayes_slope_med_quantiles'][2]:.3f} | "
                f"Failed SNR4-7 (N={f47.get('N', 0)}): >+0.09 equiv={f47.get('frac_bias_vs_inj_equiv_gt_0p09', 0.0):.3f}, "
                f">+0.09 Eq3={f47.get('frac_resid_vs_Eq3_gt_0p09', 0.0):.3f}",
                flush=True,
            )

        # 2. Sensitivity check: own-pipeline L_Ha on the 20 failed templates at sigma_parent = 0.10 dex
        ownL_summaries = {}
        for sc_id, base_name, kind, a0, b0 in [
            (11, "gauss_rel_s010", "single_gaussian", A_GAUSS, B_GAUSS),
            (13, "exp_matched_s010", "exponential", A_EXP, B_EXP),
        ]:
            own_name = f"{base_name}_ownL"
            ckpt_own = ckpt_dir / f"{own_name}_failed20.json"
            if ckpt_own.exists():
                rows_f20_own = json.loads(ckpt_own.read_text())["rows_failed20"]
                print(f"[Cached] {own_name} ({len(rows_f20_own)} failed20 rows)", flush=True)
            else:
                t0 = time.monotonic()
                tasks_own = [
                    (sc_id, kind, a0, b0, 0.10, "own", rep, tid)
                    for rep in range(100)
                    for tid in range(29, 49)
                ]
                rows_f20_own = []
                for idx_done, res in enumerate(pool.map(worker_random_draw, tasks_own, chunksize=8), 1):
                    rows_f20_own.append(res)
                    if idx_done % 1000 == 0 or idx_done == len(tasks_own):
                        print(f"[{own_name}] {idx_done}/{len(tasks_own)} fits ({time.monotonic() - t0:.1f}s)", flush=True)
                tmp = ckpt_own.with_suffix(".tmp")
                tmp.write_text(json.dumps({"rows_failed20": rows_f20_own}, ensure_ascii=False))
                tmp.rename(ckpt_own)

            rows_f29 = [r for r in cached_rows_s010[base_name] if r["template_group"] == "frozen29"]
            merged_own = rows_f29 + rows_f20_own
            ownL_summaries[own_name] = summarize_suite(own_name, kind, a0, b0, 0.10, merged_own, n_reps=100)
            p49_o = ownL_summaries[own_name]["parent49"]
            print(
                f"[DONE {own_name:22s}] Parent49 N_sel_med={p49_o['N_selected_min_p16_p50_p84_max'][2]:4.1f} | "
                f"sigma_int med={p49_o['bayes_sigma_med_quantiles'][2]:.4f} "
                f"(hits<=0.05403: {p49_o['bayes_hits_le_00540301']:2d}/100)",
                flush=True,
            )

        # 3. Paired common-random-number offset grid at delta_y in {-0.10, 0.00, +0.10} dex (N_grid = 12 reps)
        n_reps_grid = 12
        offset_grids = {}
        for g_id, g_name, kind, a0, b0 in [
            (31, "offset_grid_gauss", "single_gaussian", A_GAUSS, B_GAUSS),
            (32, "offset_grid_exp", "exponential", A_EXP, B_EXP),
        ]:
            ckpt_g = ckpt_dir / f"{g_name}.json"
            if ckpt_g.exists():
                rows_g = json.loads(ckpt_g.read_text())["rows"]
                print(f"[Cached] {g_name} ({len(rows_g)} rows)", flush=True)
            else:
                t0 = time.monotonic()
                tasks_g = [
                    (g_id, kind, a0, b0, rep, tid, dy)
                    for rep in range(n_reps_grid)
                    for tid in range(49)
                    for dy in (-0.10, 0.00, 0.10)
                ]
                rows_g = []
                for idx_done, res in enumerate(pool.map(worker_paired_offset, tasks_g, chunksize=8), 1):
                    rows_g.append(res)
                    if idx_done % 588 == 0 or idx_done == len(tasks_g):
                        print(f"[{g_name}] {idx_done}/{len(tasks_g)} fits ({time.monotonic() - t0:.1f}s)", flush=True)
                tmp = ckpt_g.with_suffix(".tmp")
                tmp.write_text(json.dumps({"rows": rows_g}, ensure_ascii=False))
                tmp.rename(ckpt_g)
            offset_grids[g_name] = summarize_paired_offset_grid(rows_g, n_reps_grid)
            og = offset_grids[g_name]
            print(
                f"[DONE {g_name:18s}] Pass rate all49 (-0.10 / 0.00 / +0.10): "
                f"{og['dy_-0.10']['pass_rate_all49']:.3f} / {og['dy_+0.00']['pass_rate_all49']:.3f} / {og['dy_+0.10']['pass_rate_all49']:.3f} "
                f"(ratio +0.10/0.00 = {og['ratio_plus010_over_zero_all49']:.3f})",
                flush=True,
            )

    # 4. Evaluate pre-registered criteria C1 and C2
    c2_hits = {
        "gauss_s010_physL_hits": suite_summaries["gauss_rel_s010"]["parent49"]["bayes_hits_le_00540301"],
        "exp_s010_physL_hits": suite_summaries["exp_matched_s010"]["parent49"]["bayes_hits_le_00540301"],
        "gauss_s010_ownL_hits": ownL_summaries["gauss_rel_s010_ownL"]["parent49"]["bayes_hits_le_00540301"],
        "exp_s010_ownL_hits": ownL_summaries["exp_matched_s010_ownL"]["parent49"]["bayes_hits_le_00540301"],
    }
    c2_pass = all(v <= 5 for v in c2_hits.values())

    # Check C1 across all 4 primary scenarios
    c1_eval = {}
    trigger_t3b = False
    for name in ("gauss_rel_s005", "gauss_rel_s010", "exp_matched_s005", "exp_matched_s010"):
        s_bin = suite_summaries[name]["failed_fits_SNR_bins_all49"]["SNR_4_to_7"]
        frac_eq = s_bin.get("frac_bias_vs_inj_equiv_gt_0p09", 0.0)
        frac_rel = s_bin.get("frac_resid_vs_Eq3_gt_0p09", 0.0)
        og_key = "offset_grid_gauss" if "gauss" in name else "offset_grid_exp"
        ratio_all = offset_grids[og_key]["ratio_plus010_over_zero_all49"]
        ratio_f20 = offset_grids[og_key]["ratio_plus010_over_zero_failed20"]
        if frac_eq >= 0.40:
            regime = "artifact_regime (>=40%)"
        elif frac_eq <= 0.20 and (ratio_all < 0.50 or ratio_f20 < 0.50):
            regime = "truncation_regime (<=20% and >50% pass-rate drop at +0.10 dex)"
            trigger_t3b = True
        else:
            regime = "intermediate_regime (report exact numbers without single-cause claim)"
        c1_eval[name] = {
            "N_failed_SNR4to7_all49": s_bin.get("N", 0),
            "frac_bias_vs_inj_equiv_gt_0p09": frac_eq,
            "frac_resid_vs_Eq3_gt_0p09": frac_rel,
            "pass_rate_ratio_plus010_over_zero_all49": ratio_all,
            "pass_rate_ratio_plus010_over_zero_failed20": ratio_f20,
            "regime": regime,
        }

    t3b_result = None
    if trigger_t3b:
        t3b_result = {
            "triggered": True,
            "gauss_grid_corrected_29LRD": run_t3b_selection_corrected_posterior(offset_grids["offset_grid_gauss"]),
            "exp_grid_corrected_29LRD": run_t3b_selection_corrected_posterior(offset_grids["offset_grid_exp"]),
        }
        print(f"[selection-corrected upper limit TRIGGERED] Selection-corrected 95% upper limits: {t3b_result}", flush=True)
    else:
        t3b_result = {"triggered": False, "reason": "Pre-registered C1 truncation condition not met."}
        print("[selection-corrected upper limit NOT triggered] Pre-registered C1 truncation condition not met.", flush=True)

    final_out = {
        "thresholds": {
            "observed_bayes_sigma_int_median_dex": OBS_SIG_MED,
            "observed_ols_residual_sd_dex": OBS_OLS_SD,
        },
        "primary_parent49_suites": suite_summaries,
        "sensitivity_ownL_parent49_suites": ownL_summaries,
        "paired_offset_pass_rate_grids": offset_grids,
        "preregistered_evaluation": {
            "C1_SNR4to7_evaluation": c1_eval,
            "C2_scatter_compression_hits_le_0054": c2_hits,
            "C2_passed_all_le_5_of_100": c2_pass,
        },
        "todo_3b_selection_corrected_posterior": t3b_result,
        "elapsed_total_s": round(time.monotonic() - t_start, 2),
    }
    out_path = RESULTS / "parent49_injection.json"
    out_path.write_text(json.dumps(final_out, indent=2, ensure_ascii=False))
    print(f"=== ALL TODO 3 SIMULATIONS & SELF-CHECKS PASSED ({final_out['elapsed_total_s']}s) -> {out_path} ===", flush=True)


if __name__ == "__main__":
    main()
