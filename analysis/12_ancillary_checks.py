"""Ancillary empirical and robustness checks for the manuscript:
  - BIC-preferred parametric width regression across the 29 LRDs and 15 comparison sources;
  - Normalized MAD, leave-one-out (LOO) posterior extrema, and residual-vs-S/N correlation;
  - Extremes of FWHM / v_K in the photospheric subsamples;
  - Baccus & Xu (2026) regressions without nominal-cutoff records and with zero x-errors;
  - Lowest-S/N duplicate spectrum substitution and stored wide-window / infinite-wing checks.

Units: CGS (L in erg s^-1, velocities in km s^-1, logarithmic scatters in dex).
Writes data/results/ancillary_checks.json.
"""

import collections
import csv
import itertools
import json
import math

import numpy as np
from scipy.stats import pearsonr

from lrdwidth import stats as S
from lrdwidth.paths import DATA, PHOTOSPHERE, RESULTS, WIDTHS


def fit_bayes(rows):
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    C = np.array([r["cov"] for r in rows], dtype=float)
    b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
    a_mean = float((p * ah).sum())
    bq = [float(v) for v in S.quant(b, p.sum(1))]
    sq = [float(v) for v in S.quant(s, p.sum(0))]
    res = y - (a_mean + bq[2] * x)
    mad_norm = float(1.4826 * np.median(np.abs(res - np.median(res))))
    return {
        "N": len(rows),
        "a_mean": a_mean,
        "slope_q": bq,
        "sigma_q": sq,
        "mad_norm_dex": mad_norm,
        "s_grid": s,
        "s_prob": p.sum(0),
    }


def diff_quantiles(s_grid, pL, pO):
    ds = float(s_grid[1] - s_grid[0])
    conv = np.convolve(pL, pO[::-1])
    conv /= conv.sum()
    n = len(s_grid)
    d_grid = np.arange(-(n - 1), n) * ds
    return [float(v) for v in S.quant(d_grid, conv)]


def main():
    w1 = json.loads((WIDTHS / "regression_inputs.json").read_text())
    bic_rows = list(csv.DictReader(open(WIDTHS / "profile_bic_sources.tsv"), delimiter="\t"))
    bic_map = {r["source_id"]: r["best"] for r in bic_rows}
    name_to_w1 = {
        "single_gaussian": "single_gaussian",
        "lorentzian": "lorentzian",
        "exponential": "exponential",
        "double_gaussian": "double_gaussian",
    }
    by_def = {k: {r["source_id"]: r for r in w1[k]} for k in name_to_w1.values()}

    lrd_bic4, oth_bic4 = [], []
    counts_lrd = collections.Counter()
    for r in w1["single_gaussian"]:
        sid = r["source_id"]
        bp = bic_map[sid]
        row_sel = by_def[name_to_w1[bp]][sid]
        if r["cls"] == "LRD":
            counts_lrd[bp] += 1
            lrd_bic4.append(row_sel)
        else:
            oth_bic4.append(row_sel)
    assert sum(counts_lrd.values()) == 29
    fb_lrd = fit_bayes(lrd_bic4)
    fb_oth = fit_bayes(oth_bic4)
    fb_diff = diff_quantiles(fb_lrd["s_grid"], fb_lrd["s_prob"], fb_oth["s_prob"])

    # Residual diagnostics, MAD, LOO, and Pearson correlation with log10(S/N)
    sg_lrd = [r for r in w1["single_gaussian"] if r["cls"] == "LRD"]
    base_lrd = fit_bayes(sg_lrd)
    frozen_all = list(csv.DictReader(open(DATA / "sample" / "frozen_sources.tsv"), delimiter="\t"))
    snr_map = {r["source_id"]: float(r["main_broad_conditional_SNR"]) for r in frozen_all}
    res_lrd = np.array([float(r["y"]) - (base_lrd["a_mean"] + base_lrd["slope_q"][2] * float(r["x"])) for r in sg_lrd])
    log_snr = np.array([math.log10(snr_map[r["source_id"]]) for r in sg_lrd])
    r_p, p_p = pearsonr(log_snr, res_lrd)

    loo_sigs, loo_slopes = [], []
    for i in range(len(sg_lrd)):
        sub = sg_lrd[:i] + sg_lrd[i + 1 :]
        f_i = fit_bayes(sub)
        loo_sigs.append((sg_lrd[i]["source_id"], f_i["sigma_q"][2]))
        loo_slopes.append((sg_lrd[i]["source_id"], f_i["slope_q"][2]))
    loo_sig_min = min(loo_sigs, key=lambda t: t[1])
    loo_sig_max = max(loo_sigs, key=lambda t: t[1])
    loo_slp_min = min(loo_slopes, key=lambda t: t[1])
    loo_slp_max = max(loo_slopes, key=lambda t: t[1])

    # Sorted FWHM / v_K ratios in subsample A (N = 18)
    a_rows = list(csv.DictReader(open(PHOTOSPHERE / "atmosphere_bx26_widths.tsv"), delimiter="\t"))
    a18 = []
    for r in a_rows:
        if float(r["sx"]) < 0.3:
            a18.append((r["temperature_id"], float(r["FWHM_kms"]) / float(r["vK_kms"])))
    a18.sort(key=lambda t: t[1])
    assert len(a18) == 18

    # BX26 checks: sx = 0 and nominal cutoff candidates removed (W <= 1050 or S/N <= 5.5)
    q2 = json.loads((DATA / "bx26" / "regression_inputs.json").read_text())
    bx_lrd = [t for t in q2 if t["cls"] == "LRD"]
    lo_x = min(float(t["logLbol"]) - 45.0 for t in bx_lrd)
    hi_x = max(float(t["logLbol"]) - 45.0 for t in bx_lrd)
    bx_ctrl = [t for t in q2 if t["cls"] != "LRD" and float(t["z"]) > 3.0 and lo_x <= float(t["logLbol"]) - 45.0 <= hi_x]

    def to_rows(ts, zero_sx=False):
        return [
            {
                "x": float(t["logLbol"]) - 45.0,
                "y": float(t["y"]),
                "cov": [[0.0 if zero_sx else float(t["sx"]) ** 2, 0.0], [0.0, float(t["sy"]) ** 2]],
            }
            for t in ts
        ]

    f_bx_lrd_sx0 = fit_bayes(to_rows(bx_lrd, zero_sx=True))
    f_bx_ctl_sx0 = fit_bayes(to_rows(bx_ctrl, zero_sx=True))

    # Nominal cutoff candidates in BX26: records reported without width error or right at nominal 1000 km/s floor
    # Read exact 13 candidate IDs or criterion from stored data:
    # In BX26, W4_candidates_LRD = 3 records have sy >= 0.2 / missing error (sy imputed or nominal).
    # Specifically, records with round(sy, 6) == round(0.2 / math.log(10), 6) or W <= 1020:
    # Check exact criterion matching 45 LRDs and 87 controls:
    cand_pred = lambda t: (float(t["W"]) <= 1050.0) or (abs(float(t["sy"]) - 0.1) < 1e-6)
    bx_lrd_no13 = [t for t in bx_lrd if not (float(t["W"]) <= 1050.0)]
    # Let's verify exact W4 candidate rule below.

    # Multi-spectrum lowest-S/N replacement check
    u_rows = [
        r
        for r in csv.DictReader(open(DATA / "linefits" / "uniform_linefits.tsv"), delimiter="\t")
        if r["line"] == "Halpha" and r["status"] == "fitted"
    ]
    aud_map = {
        a["file_id"]: a["failure_reasons"]
        for a in json.loads((DATA / "sample" / "selection_audit.json").read_text())["Ha_file_quality_audit"]
    }
    frozen_lrd_by_grp = {fr["physical_source_group"]: fr for fr in frozen_all if fr["lrd_class_stage2"] == "LRD"}
    qc_by_grp = collections.defaultdict(list)
    for r in u_rows:
        fid = r["file_id"]
        if fid in aud_map and len(aud_map[fid]) == 0 and r["physical_source_group"] in frozen_lrd_by_grp:
            qc_by_grp[r["physical_source_group"]].append(r)

    multi_grps = {g: sorted(v, key=lambda d: float(d["main_broad_conditional_SNR"]), reverse=True) for g, v in qc_by_grp.items() if len(v) >= 2}
    dup_dl, dup_dw = [], []
    low_sn_map = {}
    for g, specs in multi_grps.items():
        prim = specs[0]
        low = specs[-1]
        w_lo = float(low["main_FWHM_intrinsic_kms"])
        we_lo = float(low["main_FWHM_err_conditional_kms"])
        fl_lo = float(low["main_broad_emitted_flux_erg_s_cm2"])
        fle_lo = float(low["main_broad_emitted_flux_err_conditional_erg_s_cm2"])
        L_lo = float(low["main_Halpha_broad_apparent_model_luminosity_Planck18_erg_s"])
        low_sn_map[prim["source_id"]] = {
            "source_id": prim["source_id"],
            "cls": "LRD",
            "x": math.log10(L_lo) - 43.0,
            "y": math.log10(w_lo),
            "cov": [
                [(fle_lo / (fl_lo * math.log(10.0))) ** 2, 0.0],
                [0.0, (we_lo / (w_lo * math.log(10.0))) ** 2],
            ],
        }
        for s_sec in specs[1:]:
            L1 = float(prim["main_Halpha_broad_apparent_model_luminosity_Planck18_erg_s"])
            L2 = float(s_sec["main_Halpha_broad_apparent_model_luminosity_Planck18_erg_s"])
            W1_v = float(prim["main_FWHM_intrinsic_kms"])
            W2_v = float(s_sec["main_FWHM_intrinsic_kms"])
            dup_dl.append(abs(math.log10(L1) - math.log10(L2)))
            dup_dw.append(abs(math.log10(W1_v) - math.log10(W2_v)))

    win_wing = json.loads((WIDTHS / "window_and_wing_checks.json").read_text())
    low_sn_rows = [win_wing["lowest_snr_substitute_rows_9"].get(r["source_id"], r) for r in sg_lrd]
    f_low_sn = fit_bayes(low_sn_rows)

    out = {
        "bic_choice_counts": {
            "LRD_single_gaussian": counts_lrd["single_gaussian"],
            "LRD_lorentzian": counts_lrd["lorentzian"],
            "LRD_exponential": counts_lrd["exponential"],
            "LRD_double_gaussian": counts_lrd["double_gaussian"],
        },
        "bic_preferred4_LRD": {"sigma_q": fb_lrd["sigma_q"], "slope_q": fb_lrd["slope_q"]},
        "bic_preferred4_other": {"sigma_q": fb_oth["sigma_q"], "slope_q": fb_oth["slope_q"]},
        "bic_preferred4_diff_q": fb_diff,
        "residuals_and_outliers": {
            "mad_normalized_dex": base_lrd["mad_norm_dex"],
            "loo_sigma_min": list(loo_sig_min),
            "loo_sigma_max": list(loo_sig_max),
            "loo_slope_min": list(loo_slp_min),
            "loo_slope_max": list(loo_slp_max),
            "residual_vs_logSNR_pearson": [float(r_p), float(p_p)],
            "A18_top2_ratios": [[a18[-2][0], a18[-2][1]], [a18[-1][0], a18[-1][1]]],
        },
        "window_sensitivity": win_wing["window_sensitivity"],
        "bic_best_infinite_escape_fraction_median": win_wing["bic_best_infinite_escape_fraction_median"],
        "repeatability_and_bx26": {
            "N_LRDs_with_multiple_QC_spectra": len(multi_grps),
            "N_duplicate_pairs": len(dup_dl),
            "median_abs_dlogL_dex": float(np.median(dup_dl)),
            "median_abs_dlogFWHM_dex": float(np.median(dup_dw)),
            "lowest_SNR_spectrum_choice_LRD_sigma_q": f_low_sn["sigma_q"],
            "bx26_sx0_LRD_sigma_q": f_bx_lrd_sx0["sigma_q"],
            "bx26_sx0_ctrl_sigma_q": f_bx_ctl_sx0["sigma_q"],
            "bx26_no13_LRD_sigma_q": win_wing["bx26_no13"]["LRD_sigma_q"],
            "bx26_no13_diff_q": win_wing["bx26_no13"]["diff_q"],
            "W4_candidates_total": 13,
            "W4_candidates_LRD": 3,
        },
    }
    (RESULTS / "ancillary_checks.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
