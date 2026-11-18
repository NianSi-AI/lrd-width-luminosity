"""Calibrate unmodeled per-spectrum width and flux variance from the 9 multi-spectrum LRDs,
re-evaluate the Bayesian intrinsic-scatter posterior with empirical error floors, and compute
the high-S/N residual dispersion, Baccus & Xu (2026) sky deduplication, refitted SDSS width
offset, and photospheric-to-line luminosity dispersion.

Units: CGS throughout (luminosities in erg s^-1, velocities in km s^-1, dispersions in dex).
Writes data/results/repeat_calibration.json and data/widths/repeat_pairs_22.tsv.
"""

import collections
import csv
import itertools
import json
import math

import numpy as np
from scipy.optimize import brentq

from lrdwidth import stats as S
from lrdwidth.paths import DATA, PHOTOSPHERE, RESULTS, WIDTHS


def fit_bayes(rows, extra_sy_dex: float = 0.0, extra_sx_dex: float = 0.0):
    """Evaluate exact 2D posterior over (b, sigma_int) with intercept marginalized analytically."""
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    C = np.array([r["cov"] for r in rows], dtype=float)
    cxx = C[:, 0, 0] + extra_sx_dex**2
    cyy = C[:, 1, 1] + extra_sy_dex**2
    cxy = C[:, 0, 1]
    b, s, p, ah = S.post2(x, y, cxx, cyy, cxy)
    a_mean = float((p * ah).sum())
    ps = p.sum(0)
    bq = [float(v) for v in S.quant(b, p.sum(1))]
    sq = [float(v) for v in S.quant(s, ps)]
    s_p95_onesided = float(S.quant(s, ps, qs=(0.95,))[0])
    co = np.polyfit(x, y, 1)
    ols_res = y - np.polyval(co, x)
    return {
        "N": len(rows),
        "a_mean": a_mean,
        "slope_q": bq,
        "sigma_q": sq,
        "sigma_p95_onesided": s_p95_onesided,
        "ols_sd_ddof1": float(np.std(ols_res, ddof=1)),
        "ols_sd_ddof2": float(np.std(ols_res, ddof=2)),
        "s_grid": s,
        "s_prob": ps,
    }


def diff_quantiles(s_grid, pL, pO):
    """Posterior quantiles of Delta sigma_int = sigma_int(LRD) - sigma_int(other)."""
    ds = float(s_grid[1] - s_grid[0])
    conv = np.convolve(pL, pO[::-1])
    conv /= conv.sum()
    n = len(s_grid)
    d_grid = np.arange(-(n - 1), n) * ds
    return [float(v) for v in S.quant(d_grid, conv)]


def solve_extra_floor_per_spec(pairs):
    """Solve mean(Delta^2 / (s1^2 + s2^2 + 2 * s_floor^2)) = 1 for per-spectrum s_floor (dex)."""
    def obj(sf):
        return float(np.mean([(d**2) / (s1**2 + s2**2 + 2.0 * sf**2) for d, s1, s2 in pairs])) - 1.0

    if obj(0.0) <= 0.0:
        return 0.0
    return float(brentq(obj, 0.0, 1.0))


def summarize_pairs(plist):
    zs = [d / math.sqrt(s1**2 + s2**2) for d, s1, s2 in plist]
    ds = [abs(d) for d, _, _ in plist]
    return {
        "N_pairs": len(plist),
        "norm_diff_rms": float(np.sqrt(np.mean(np.square(zs)))),
        "extra_floor_per_spec_dex": solve_extra_floor_per_spec(plist),
        "raw_pair_rms_over_sqrt2_dex": float(np.sqrt(np.mean([d**2 for d, _, _ in plist])) / math.sqrt(2.0)),
        "median_abs_diff_dex": float(np.median(ds)),
    }


def main():
    frozen_all = list(csv.DictReader(open(DATA / "sample" / "frozen_sources.tsv"), delimiter="\t"))
    w1 = json.loads((WIDTHS / "regression_inputs.json").read_text())
    sg_lrd = [dict(r) for r in w1["single_gaussian"] if r["cls"] == "LRD"]
    sg_oth = [dict(r) for r in w1["single_gaussian"] if r["cls"] != "LRD"]
    assert len(sg_lrd) == 29 and len(sg_oth) == 15

    snr_map = {r["source_id"]: float(r["main_broad_conditional_SNR"]) for r in frozen_all}
    for r in sg_lrd + sg_oth:
        r["snr"] = snr_map[r["source_id"]]

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

    def parse_spec(r):
        w = float(r["main_FWHM_intrinsic_kms"])
        we = float(r["main_FWHM_err_conditional_kms"])
        fl = float(r["main_broad_emitted_flux_erg_s_cm2"])
        fle = float(r["main_broad_emitted_flux_err_conditional_erg_s_cm2"])
        L = float(r["main_Halpha_broad_apparent_model_luminosity_Planck18_erg_s"])
        return {
            "file_id": r["file_id"],
            "file": r["file"],
            "source_id": r["source_id"],
            "group": r["physical_source_group"],
            "grating": r["grating"],
            "snr": float(r["main_broad_conditional_SNR"]),
            "y": math.log10(w),
            "sy": we / (w * math.log(10.0)),
            "x": math.log10(L),
            "sx": fle / (fl * math.log(10.0)),
            "FWHM_kms": w,
            "FWHM_err_kms": we,
        }

    qc_by_grp = collections.defaultdict(list)
    for r in u_rows:
        fid = r["file_id"]
        if fid in aud_map and len(aud_map[fid]) == 0 and r["physical_source_group"] in frozen_lrd_by_grp:
            qc_by_grp[r["physical_source_group"]].append(parse_spec(r))

    multi_lrd = {g: sorted(v, key=lambda d: d["snr"], reverse=True) for g, v in qc_by_grp.items() if len(v) >= 2}
    assert len(multi_lrd) == 9

    best_rest_w, all_pairs_w, same_grat_w, cross_grat_w = [], [], [], []
    best_rest_x, all_pairs_x, same_grat_x = [], [], []
    tsv_rows = []

    for g, specs in sorted(multi_lrd.items()):
        prim = specs[0]
        assert prim["file_id"] == frozen_lrd_by_grp[g]["file_id"]
        for sec in specs[1:]:
            best_rest_w.append((prim["y"] - sec["y"], prim["sy"], sec["sy"]))
            best_rest_x.append((prim["x"] - sec["x"], prim["sx"], sec["sx"]))
        for idx_a, idx_b in itertools.combinations(range(len(specs)), 2):
            a, b = specs[idx_a], specs[idx_b]
            dy = a["y"] - b["y"]
            dx = a["x"] - b["x"]
            same_g = a["grating"] == b["grating"]
            all_pairs_w.append((dy, a["sy"], b["sy"]))
            all_pairs_x.append((dx, a["sx"], b["sx"]))
            if same_g:
                same_grat_w.append((dy, a["sy"], b["sy"]))
                same_grat_x.append((dx, a["sx"], b["sx"]))
            else:
                cross_grat_w.append((dy, a["sy"], b["sy"]))
            tsv_rows.append(
                {
                    "source_id": prim["source_id"],
                    "is_primary_pair": idx_a == 0,
                    "same_grating": same_g,
                    "grating_a": a["grating"],
                    "grating_b": b["grating"],
                    "snr_a": f"{a['snr']:.2f}",
                    "snr_b": f"{b['snr']:.2f}",
                    "FWHM_a_kms": f"{a['FWHM_kms']:.1f}",
                    "FWHM_b_kms": f"{b['FWHM_kms']:.1f}",
                    "dlogFWHM_dex": f"{dy:.5f}",
                    "formal_pair_sigma_dex": f"{math.sqrt(a['sy']**2 + b['sy']**2):.5f}",
                    "norm_diff_z": f"{dy / math.sqrt(a['sy']**2 + b['sy']**2):.3f}",
                    "dlogLHa_dex": f"{dx:.5f}",
                }
            )

    sw_best = summarize_pairs(best_rest_w)
    sw_all = summarize_pairs(all_pairs_w)
    sw_same = summarize_pairs(same_grat_w)
    sw_cross = summarize_pairs(cross_grat_w)
    sx_best = summarize_pairs(best_rest_x)
    sx_all = summarize_pairs(all_pairs_x)
    sx_same = summarize_pairs(same_grat_x)

    assert sw_best["N_pairs"] == 14 and abs(sw_best["norm_diff_rms"] - 2.81) < 0.01
    assert abs(sw_all["extra_floor_per_spec_dex"] - 0.060) < 0.001
    assert abs(sw_same["extra_floor_per_spec_dex"] - 0.044) < 0.001

    with open(WIDTHS / "repeat_pairs_22.tsv", "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(tsv_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(tsv_rows)

    fl_same = sw_same["extra_floor_per_spec_dex"]
    fl_all = sw_all["extra_floor_per_spec_dex"]
    fl_best = sw_best["extra_floor_per_spec_dex"]

    post_formal = fit_bayes(sg_lrd, 0.0, 0.0)

    def run_floor_pair(sy_f, sx_f=0.0):
        fL = fit_bayes(sg_lrd, sy_f, sx_f)
        fO = fit_bayes(sg_oth, sy_f, sx_f)
        dq = diff_quantiles(fL["s_grid"], fL["s_prob"], fO["s_prob"])
        return {
            "sy_floor_dex": sy_f,
            "sx_floor_dex": sx_f,
            "LRD_sigma_q": fL["sigma_q"],
            "LRD_p95_onesided": fL["sigma_p95_onesided"],
            "LRD_slope_q": fL["slope_q"],
            "other_sigma_q": fO["sigma_q"],
            "diff_q": dq,
        }

    refit_formal = run_floor_pair(0.0)
    refit_same = run_floor_pair(0.044)
    refit_all = run_floor_pair(0.060)
    refit_same_unrounded = run_floor_pair(fl_same)
    refit_all_unrounded = run_floor_pair(fl_all)
    refit_best = run_floor_pair(fl_best)
    refit_same_xy = run_floor_pair(0.044, sx_same["extra_floor_per_spec_dex"])
    refit_all_xy = run_floor_pair(0.060, sx_all["extra_floor_per_spec_dex"])

    assert abs(refit_same["LRD_sigma_q"][2] - 0.027) < 0.0005
    assert abs(refit_all["LRD_sigma_q"][2] - 0.018) < 0.0005

    lrd_sorted_L = sorted(sg_lrd, key=lambda r: float(r["x"]))
    post_noTop3 = fit_bayes(lrd_sorted_L[:-3], 0.0, 0.0)

    sn15_rows = [r for r in sg_lrd if r["snr"] >= 15.0]
    assert len(sn15_rows) == 12
    a_full = post_formal["a_mean"]
    b_full = post_formal["slope_q"][2]
    res_sn15_full = np.array([float(r["y"]) - (a_full + b_full * float(r["x"])) for r in sn15_rows])
    post_sn15 = fit_bayes(sn15_rows)
    sn15_check = {
        "N": len(sn15_rows),
        "raw_ols_sd_ddof1_dex": post_sn15["ols_sd_ddof1"],
        "raw_ols_sd_ddof2_dex": post_sn15["ols_sd_ddof2"],
        "raw_rms_about_full_relation_ddof2_dex": float(np.sqrt(np.sum(res_sn15_full**2) / (len(sn15_rows) - 2))),
        "sigma_q": post_sn15["sigma_q"],
    }
    assert abs(sn15_check["raw_ols_sd_ddof2_dex"] - 0.056) < 0.001

    q2_in = json.loads((DATA / "bx26" / "regression_inputs.json").read_text())
    bx_lrd_raw = [t for t in q2_in if t["cls"] == "LRD"]
    clusters = []
    for t in bx_lrd_raw:
        ra, dec = float(t["ra"]), float(t["dec"])
        matched = False
        for cl in clusters:
            ra0, dec0 = float(cl[0]["ra"]), float(cl[0]["dec"])
            sep = math.hypot((ra - ra0) * math.cos(math.radians(dec0)), dec - dec0) * 3600.0
            if sep < 1.0:
                cl.append(t)
                matched = True
                break
        if not matched:
            clusters.append([t])
    assert len(bx_lrd_raw) == 48 and len(clusters) == 45
    bx_lrd_dedup = [
        {
            "id": min(cl, key=lambda r: float(r["sy"]))["id"],
            "x": float(min(cl, key=lambda r: float(r["sy"]))["logLbol"]) - 45.0,
            "y": float(min(cl, key=lambda r: float(r["sy"]))["y"]),
            "cov": [
                [float(min(cl, key=lambda r: float(r["sy"]))["sx"]) ** 2, 0.0],
                [0.0, float(min(cl, key=lambda r: float(r["sy"]))["sy"]) ** 2],
            ],
        }
        for cl in clusters
    ]
    post_bx_dedup = fit_bayes(bx_lrd_dedup)
    assert abs(post_bx_dedup["sigma_q"][2] - 0.060) < 0.001

    q2_ov = list(csv.DictReader(open(DATA / "bx26" / "overlap_with_this_work.tsv"), delimiter="\t"))
    w_ratios_all = [float(r["log10_width_ratio"]) for r in q2_ov]
    w_ratios_lrd = [
        float(r["log10_width_ratio"])
        for r in q2_ov
        if r["published_class"] == "LRD" and r["frozen_class"] == "LRD"
    ]

    sdss_in = json.loads((DATA / "sdss" / "regression_inputs.json").read_text())["single_gaussian_QC"]
    post_sdss = fit_bayes(sdss_in)
    x_med_lrd = float(np.median([float(r["x"]) for r in sg_lrd]))
    logW_sdss_at_med = post_sdss["a_mean"] + post_sdss["slope_q"][2] * x_med_lrd
    logW_lrd_at_med = post_formal["a_mean"] + post_formal["slope_q"][2] * x_med_lrd
    dlogW_sdss_lrd = logW_sdss_at_med - logW_lrd_at_med
    assert abs(dlogW_sdss_lrd - 0.264) < 0.005

    ds_samples = json.loads((PHOTOSPHERE / "samples.json").read_text())
    lph_ha = {}
    for key in ("Curtis_A_own", "Curtis_B_own"):
        diffs = [(float(r["x"]) + 45.0) - (float(r["line_x"]) + 43.0) for r in ds_samples[key]]
        lph_ha[key] = {
            "N": len(diffs),
            "mean_dex": float(np.mean(diffs)),
            "sd_dex": float(np.std(diffs, ddof=1)),
        }

    out = {
        "duplicate_spectra": {
            "N_multi_LRDs": 9,
            "width_best_vs_rest_14": sw_best,
            "width_all_pairs_22": sw_all,
            "width_same_grating_7": sw_same,
            "width_cross_grating_15": sw_cross,
            "flux_best_vs_rest_14": sx_best,
            "flux_all_pairs_22": sx_all,
            "flux_same_grating_7": sx_same,
        },
        "posterior_with_floors": {
            "formal_zero_floor": refit_formal,
            "same_grating_sy_floor": refit_same,
            "all_pairs_sy_floor": refit_all,
            "same_grating_unrounded_sy_floor": refit_same_unrounded,
            "all_pairs_unrounded_sy_floor": refit_all_unrounded,
            "best_vs_rest_sy_floor": refit_best,
            "same_grating_xy_floor": refit_same_xy,
            "all_pairs_xy_floor": refit_all_xy,
            "noTop3_formal_p95_onesided": post_noTop3["sigma_p95_onesided"],
        },
        "high_snr_LRDs_ge15": sn15_check,
        "bx26_checks": {
            "N_records_LRD": 48,
            "N_unique_sky_LRD": len(clusters),
            "duplicate_clusters": [[x["id"] for x in cl] for cl in clusters if len(cl) > 1],
            "dedup_LRD_sigma_q": post_bx_dedup["sigma_q"],
            "shared_29_width_ratio_sd_dex": float(np.std(w_ratios_all, ddof=1)),
            "shared_29_width_ratio_mean_dex": float(np.mean(w_ratios_all)),
            "shared_16_LRD_width_ratio_sd_dex": float(np.std(w_ratios_lrd, ddof=1)),
        },
        "sdss_224_offset": {
            "N": len(sdss_in),
            "LRD_median_logLHa": x_med_lrd + 43.0,
            "SDSS_FWHM_at_LRD_med_kms": 10.0**logW_sdss_at_med,
            "LRD_FWHM_at_LRD_med_kms": 10.0**logW_lrd_at_med,
            "dlogFWHM_SDSS_minus_LRD_dex": dlogW_sdss_lrd,
        },
        "photosphere_log_Lph_over_LHa": lph_ha,
    }
    (RESULTS / "repeat_calibration.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
