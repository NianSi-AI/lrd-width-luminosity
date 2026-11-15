"""Parent-sample covariance regressions, repeat-spectrum random-effects models,
and exterior-envelope gravitational potential bounds:
  1. Full 2x2 fitting-covariance regressions across the clean-label parent LRD sample
     (N = 29 primary, N = 31 adding high-S/N absorption-bound sources, N = 35 with S/N >= 4,
     N = 38 with S/N >= 3, and N = 39 with S/N >= 3 or Delta BIC > 0);
  2. Empirical per-spectrum width error floor across 6 repeat-pair slicings, leave-one-object-out
     floor sensitivity across the 9 multi-spectrum LRDs, and bivariate random-effects hierarchical
     posteriors on all 43 QC-passing spectra (and 39 spectra when RUBIES-EGS-42232 is reduced to
     its primary spectrum);
  3. Dimensionless exterior-envelope gravitational potential ratio delta = |Phi_ext(R_ph)| / (G M_BH / R_ph),
     fractional escape-speed shifts, and relaxed bound-orbit radius ceilings across the published
     photospheric temperature samples (A, B, C, E).

Units: CGS luminosities (erg s^-1, Planck18), velocities in km s^-1, logarithmic scatters in dex.
Writes:
  - data/results/parent_covariance_audit.json
  - data/results/parent_failed21_audit.tsv
  - data/results/repeat_random_effects.json
  - data/results/repeat_floor_slicings.tsv
  - data/results/exterior_potential_check.json
"""

import csv
import itertools
import json
import math
import os

for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "1")

import numpy as np
from scipy.optimize import brentq

from lrdwidth import stats as S
from lrdwidth.paths import PHOTOSPHERE, RESULTS, WIDTHS

G_CGS = 6.67430e-8
C_CGS = 2.99792458e10
SIGMA_SB_CGS = 5.670374419e-5
KAPPA_ES_MS = 0.398
KAPPA_ES_EXACT = 6.6524587321e-25 / 1.67262192369e-24


def run_bayes_post2(rows):
    """Evaluate exact 2D grid posterior over (b, sigma_int) with full 2x2 covariance."""
    x = np.array([float(r["x"]) for r in rows], dtype=float)
    y = np.array([float(r["y"]) for r in rows], dtype=float)
    C = np.array([r["cov"] for r in rows], dtype=float)
    cxx, cyy, cxy = C[:, 0, 0], C[:, 1, 1], C[:, 0, 1]
    b, s, p, ah = S.post2(x, y, cxx, cyy, cxy)
    a_mean = float((p * ah).sum())
    ps = p.sum(0)
    bq = [float(v) for v in S.quant(b, p.sum(1))]
    sq = [float(v) for v in S.quant(s, ps)]
    s_p95 = float(S.quant(s, ps, qs=(0.95,))[0])
    b0, s0, p0, _ = S.post2(x, y, cxx, cyy, np.zeros_like(cxy))
    ps0 = p0.sum(0)
    sq0 = [float(v) for v in S.quant(s0, ps0)]
    s0_p95 = float(S.quant(s0, ps0, qs=(0.95,))[0])
    bq0 = [float(v) for v in S.quant(b0, p0.sum(1))]
    co = np.polyfit(x, y, 1)
    ols_res = y - np.polyval(co, x)
    return {
        "N": len(rows),
        "a_mean": a_mean,
        "slope_q": bq,
        "sigma_q": sq,
        "sigma_p95_onesided": s_p95,
        "raw_ols_sd_ddof1": float(np.std(ols_res, ddof=1)),
        "zero_cov_control": {
            "sigma_q": sq0,
            "sigma_p95_onesided": s0_p95,
            "slope_med": bq0[2],
        },
    }


def solve_weighted_floor(pairs, weights=None):
    w = np.ones(len(pairs), dtype=float) if weights is None else np.asarray(weights, dtype=float)
    w = w / w.sum()
    d2 = np.array([p[0] ** 2 for p in pairs], dtype=float)
    v0 = np.array([p[1] ** 2 + p[2] ** 2 for p in pairs], dtype=float)

    def obj(f):
        return float(np.sum(w * d2 / (v0 + 2.0 * f * f))) - 1.0

    if obj(0.0) <= 0.0:
        return 0.0
    return float(brentq(obj, 0.0, 1.0, xtol=1e-12))


def fit_29_with_floor(sg_lrd, sy_floor_dex: float):
    x = np.array([float(r["x"]) for r in sg_lrd], dtype=float)
    y = np.array([float(r["y"]) for r in sg_lrd], dtype=float)
    C = np.array([r["cov"] for r in sg_lrd], dtype=float)
    b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1] + sy_floor_dex ** 2, C[:, 0, 1])
    ps = p.sum(0)
    sq = [float(v) for v in S.quant(s, ps)]
    bq = [float(v) for v in S.quant(b, p.sum(1))]
    p95 = float(S.quant(s, ps, qs=(0.95,))[0])
    return {
        "sy_floor_dex": float(sy_floor_dex),
        "sigma_q": sq,
        "sigma_med": sq[2],
        "sigma_p95_onesided": p95,
        "slope_q": bq,
        "a_mean": float((p * ah).sum()),
    }


def ranef_fixed_primary_exact(groups, use_cxy: bool = True, nb=161, ns=200, nf=200, f_fixed=None):
    b_grid = -0.2 + (np.arange(nb) + 0.5) * (0.8 / nb)
    s_grid = (np.arange(ns) + 0.5) * (0.20 / ns)
    f_grid = np.array([float(f_fixed)], dtype=float) if f_fixed is not None else (np.arange(nf) + 0.5) * (0.20 / nf)
    n_f = len(f_grid)
    B = b_grid[:, None, None]
    Sg = s_grid[None, :, None]
    Fg = f_grid[None, None, :]

    sw_a = np.zeros((nb, ns, n_f), dtype=float)
    swr_a = np.zeros((nb, ns, n_f), dtype=float)
    quad_0 = np.zeros((nb, ns, n_f), dtype=float)
    logdet = np.zeros((nb, ns, n_f), dtype=float)

    for _, specs in groups.items():
        m = len(specs)
        sp0 = specs[0]
        cxx0 = float(sp0["cov"][0][0])
        cxy0 = float(sp0["cov"][0][1]) if use_cxy else 0.0
        cyy0 = float(sp0["cov"][1][1])
        x0 = float(sp0["x"])
        if m == 1:
            V = Sg * Sg + B * B * cxx0 - 2.0 * B * cxy0 + cyy0 + Fg * Fg
            d = float(sp0["y"]) - B * x0
            w = 1.0 / V
            sw_a += w
            swr_a += w * d
            quad_0 += w * d * d
            logdet += np.log(V)
        else:
            y_vec = np.array([float(sp["y"]) for sp in specs], dtype=float)
            cyy_vec = np.array([float(sp["cov"][1][1]) for sp in specs], dtype=float)
            step = 25
            for i0 in range(0, nb, step):
                i1 = min(nb, i0 + step)
                bs = b_grid[i0:i1]
                B_sub = bs[:, None, None]
                nb_sub = i1 - i0
                base_var = Sg * Sg + B_sub * B_sub * cxx0
                Sigma = np.broadcast_to(base_var[..., None, None], (nb_sub, ns, n_f, m, m)).copy()
                Sigma[..., 0, :] -= (B_sub * cxy0)[..., None]
                Sigma[..., :, 0] -= (B_sub * cxy0)[..., None]
                Sigma[..., 0, 0] += (B_sub * cxy0)
                for k in range(m):
                    Sigma[..., k, k] += cyy_vec[k] + Fg * Fg
                _, ldet = np.linalg.slogdet(Sigma)
                d_vec = np.broadcast_to((y_vec[None, :] - bs[:, None] * x0)[:, None, None, :, None], (nb_sub, ns, n_f, m, 1))
                ones = np.ones((nb_sub, ns, n_f, m, 1), dtype=float)
                sol = np.linalg.solve(Sigma, np.concatenate([ones, d_vec], axis=-1))
                u_ones, u_d = sol[..., 0], sol[..., 1]
                sw_a[i0:i1] += np.sum(u_ones, axis=-1)
                swr_a[i0:i1] += np.sum(u_d, axis=-1)
                quad_0[i0:i1] += np.sum(d_vec[..., 0] * u_d, axis=-1)
                logdet[i0:i1] += ldet

    ah = swr_a / sw_a
    asd = 1.0 / np.sqrt(sw_a)
    norm = S.ndtr((4.5 - ah) / asd) - S.ndtr((2.0 - ah) / asd)
    lp = -0.5 * (logdet + quad_0 - swr_a * swr_a / sw_a) - 0.5 * np.log(sw_a) + np.log(np.maximum(norm, 1e-300))
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    ps = p.sum((0, 2))
    pb = p.sum((1, 2))
    out = {
        "N_objects": len(groups),
        "N_spectra": int(sum(len(v) for v in groups.values())),
        "sigma_q": [float(v) for v in S.quant(s_grid, ps)],
        "sigma_med": float(S.quant(s_grid, ps, qs=(0.5,))[0]),
        "sigma_p95_onesided": float(S.quant(s_grid, ps, qs=(0.95,))[0]),
        "slope_q": [float(v) for v in S.quant(b_grid, pb)],
    }
    if f_fixed is None:
        pf = p.sum((0, 1))
        fq = [float(v) for v in S.quant(f_grid, pf)]
        out["f_q"] = fq
        out["f_med"] = fq[2]
        out["f_95_interval"] = [fq[0], fq[4]]
    return out


def ranef_per_spectrum_residual_2dcov(groups, nb=161, ns=200, nf=200):
    b_grid = -0.2 + (np.arange(nb) + 0.5) * (0.8 / nb)
    s_grid = (np.arange(ns) + 0.5) * (0.20 / ns)
    f_grid = (np.arange(nf) + 0.5) * (0.20 / nf)
    B, Sg, Fg = b_grid[:, None, None], s_grid[None, :, None], f_grid[None, None, :]
    sw_a = swr_a = quad_0 = logdet = 0.0
    for _, specs in groups.items():
        W_i = D_i = Q_i = ldet_i = 0.0
        for sp in specs:
            cxx, cxy, cyy = float(sp["cov"][0][0]), float(sp["cov"][0][1]), float(sp["cov"][1][1])
            V = cyy + B * B * cxx - 2.0 * B * cxy + Fg * Fg
            w = 1.0 / V
            d = float(sp["y"]) - B * float(sp["x"])
            W_i = W_i + w
            D_i = D_i + w * d
            Q_i = Q_i + w * d * d
            ldet_i = ldet_i + np.log(V)
        denom = 1.0 + Sg * Sg * W_i
        sw_a = sw_a + W_i / denom
        swr_a = swr_a + D_i / denom
        quad_0 = quad_0 + (Q_i - Sg * Sg * D_i * D_i / denom)
        logdet = logdet + ldet_i + np.log(denom)
    ah = swr_a / sw_a
    asd = 1.0 / np.sqrt(sw_a)
    norm = S.ndtr((4.5 - ah) / asd) - S.ndtr((2.0 - ah) / asd)
    lp = -0.5 * (logdet + quad_0 - swr_a * swr_a / sw_a) - 0.5 * np.log(sw_a) + np.log(np.maximum(norm, 1e-300))
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    ps, pf, pb = p.sum((0, 2)), p.sum((0, 1)), p.sum((1, 2))
    fq = [float(v) for v in S.quant(f_grid, pf)]
    return {
        "N_objects": len(groups),
        "N_spectra": int(sum(len(v) for v in groups.values())),
        "sigma_q": [float(v) for v in S.quant(s_grid, ps)],
        "sigma_med": float(S.quant(s_grid, ps, qs=(0.5,))[0]),
        "sigma_p95_onesided": float(S.quant(s_grid, ps, qs=(0.95,))[0]),
        "f_q": fq,
        "f_med": fq[2],
        "f_95_interval": [fq[0], fq[4]],
        "slope_q": [float(v) for v in S.quant(b_grid, pb)],
    }


def delta_m_over_medd(T_eff_K: float, kappa_cgs: float, r_out_over_Rph: float = 100.0, kappa_es: float = KAPPA_ES_MS, exact_finite_rout: bool = False) -> float:
    pref = (4.0 * math.pi * G_CGS * C_CGS) / (SIGMA_SB_CGS * (T_eff_K ** 4) * kappa_es * kappa_cgs)
    val = pref * math.log(r_out_over_Rph)
    if exact_finite_rout:
        val /= (1.0 - 1.0 / r_out_over_Rph)
    return val


def delta_vesc_frac(dm_over_medd: float, lam: float = 1.0) -> float:
    return math.sqrt(1.0 + lam * dm_over_medd) - 1.0


def delta_m_general_powerlaw(T_eff_K: float, kappa_cgs: float, alpha: float, r_out_over_Rph: float = 100.0, kappa_es: float = KAPPA_ES_MS) -> float:
    u = r_out_over_Rph
    if abs(alpha - 1.0) < 1e-12:
        num, den = u - 1.0, math.log(u)
    elif abs(alpha - 2.0) < 1e-12:
        num, den = math.log(u), 1.0 - u ** (-1.0)
    else:
        num = (u ** (2.0 - alpha) - 1.0) / (2.0 - alpha)
        den = (u ** (1.0 - alpha) - 1.0) / (1.0 - alpha)
    pref = (4.0 * math.pi * G_CGS * C_CGS) / (SIGMA_SB_CGS * (T_eff_K ** 4) * kappa_es * kappa_cgs)
    return pref * (num / den)


def main():
    w1 = json.loads((WIDTHS / "regression_inputs.json").read_text())
    sg_lrd = [r for r in w1["single_gaussian"] if r["cls"] == "LRD"]
    assert len(sg_lrd) == 29

    bundle = json.loads((WIDTHS / "parent_and_repeat_covariances.json").read_text())
    spec_by_grp = bundle["qc_spectra_by_group_29_43"]
    failed_audit = bundle["failed_clean_lrd_21"]

    # --- Part 1: Parent-sample covariance audit ---
    base_29 = [v[0] for v in spec_by_grp.values()]
    post_29 = run_bayes_post2(base_29)
    assert abs(post_29["sigma_q"][2] - 0.0540268) < 1e-5
    assert abs(post_29["sigma_p95_onesided"] - 0.0747209) < 1e-5
    assert abs(post_29["slope_q"][2] - 0.1700390) < 1e-5

    out_of_band = [r for r in failed_audit if not r["in_passband"]]
    valid_failed = [r for r in failed_audit if r["in_passband"]]
    assert len(out_of_band) == 1 and out_of_band[0]["source_id"] == "RUBIES-UDS-40579"
    assert len(valid_failed) == 20

    abs_only_2 = [r for r in valid_failed if r["M1_failure_reasons"] == ["adopted_absorption_parameters_not_at_bounds"]]
    add_sn4 = [r for r in valid_failed if r["SNR_broad"] >= 4.0 and "W_broad1" not in r["bound_parameters"]]
    add_sn3 = [r for r in valid_failed if r["SNR_broad"] >= 3.0 and "W_broad1" not in r["bound_parameters"]]
    add_sn3_or_bic = [r for r in valid_failed if (r["SNR_broad"] >= 3.0 or r["delta_BIC_null"] > 0.0) and "W_broad1" not in r["bound_parameters"]]
    assert len(abs_only_2) == 2 and len(add_sn4) == 6 and len(add_sn3) == 9 and len(add_sn3_or_bic) == 10

    reg_results = {
        "row1_frozen_29": post_29,
        "row2_plus_2_abs_bound_31": run_bayes_post2(base_29 + abs_only_2),
        "row3_plus_valid_SNge4_35": run_bayes_post2(base_29 + add_sn4),
        "row4_plus_valid_SNge3_38": run_bayes_post2(base_29 + add_sn3),
        "row5_ancillary_SNge3_or_dBICgt0_39": run_bayes_post2(base_29 + add_sn3_or_bic),
    }
    tsv_fields = [
        "source_id", "file_id", "grating", "z_fit", "lambda_Ha_obs_um", "in_passband",
        "err_Ha_over_inband_median", "SNR_broad", "delta_BIC_null", "FWHM_kms",
        "log10_LHa_erg_s", "residual_Eq3_dex", "W_narrow_kms", "W_narrow_bound_side",
        "bound_parameters", "M1_failure_reasons", "literature_FWHM_kms", "literature_log10_LHa_erg_s",
    ]
    with (RESULTS / "parent_failed21_audit.tsv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=tsv_fields, delimiter="\t")
        w.writeheader()
        for r in failed_audit:
            row_out = {k: r[k] for k in tsv_fields}
            row_out["SNR_broad"] = f"{r['SNR_broad']:.2f}"
            row_out["delta_BIC_null"] = f"{r['delta_BIC_null']:.2f}"
            row_out["FWHM_kms"] = f"{r['FWHM_kms']:.1f}"
            row_out["log10_LHa_erg_s"] = f"{r['log10_LHa_erg_s']:.3f}"
            row_out["residual_Eq3_dex"] = f"{r['residual_Eq3_dex']:+.4f}"
            row_out["W_narrow_kms"] = f"{r['W_narrow_kms']:.1f}"
            row_out["bound_parameters"] = ",".join(r["bound_parameters"])
            row_out["M1_failure_reasons"] = ",".join(r["M1_failure_reasons"])
            row_out["literature_FWHM_kms"] = f"{r['literature_FWHM_kms']:.1f}" if r["literature_FWHM_kms"] else ""
            row_out["literature_log10_LHa_erg_s"] = f"{r['literature_log10_LHa_erg_s']:.3f}" if r["literature_log10_LHa_erg_s"] else ""
            w.writerow(row_out)

    parent_summary = {
        "N_parent_clean_total": 50,
        "N_frozen_LRD": 29,
        "N_failed_clean_LRD": 21,
        "N_out_of_band_excluded": len(out_of_band),
        "out_of_band_sources": [r["source_id"] for r in out_of_band],
        "N_parent_valid_inband": 49,
        "regressions_full_2x2_covariance": reg_results,
        "failed_21_audit": failed_audit,
    }
    (RESULTS / "parent_covariance_audit.json").write_text(json.dumps(parent_summary, indent=2) + "\n")

    # --- Part 2: Repeat-spectrum reweighting and 43-spectrum random-effects models ---
    prim_only_grps = {g: [v[0]] for g, v in spec_by_grp.items()}
    re_zero_check = ranef_fixed_primary_exact(prim_only_grps, use_cxy=True, nb=401, ns=600, f_fixed=0.0)
    base_check = fit_29_with_floor(sg_lrd, 0.0)
    assert abs(re_zero_check["sigma_med"] - base_check["sigma_med"]) < 1e-4

    multi_by_grp = {g: v for g, v in sorted(spec_by_grp.items()) if len(v) >= 2}
    all_pair_records = []
    for g, specs in multi_by_grp.items():
        sid = specs[0]["source_id"]
        for i, j in itertools.combinations(range(len(specs)), 2):
            a, b = specs[i], specs[j]
            all_pair_records.append({
                "source_id": sid,
                "group": g,
                "idx_a": i,
                "idx_b": j,
                "same_grating": a["grating"] == b["grating"],
                "dy": a["y"] - b["y"],
                "s1": a["sy_cond"],
                "s2": b["sy_cond"],
                "n_pairs_in_obj": len(specs) * (len(specs) - 1) // 2,
            })
    assert len(all_pair_records) == 22

    def eval_slicing(recs, weights=None):
        plist = [(r["dy"], r["s1"], r["s2"]) for r in recs]
        f_val = solve_weighted_floor(plist, weights=weights)
        res = fit_29_with_floor(sg_lrd, f_val)
        return {"N_pairs": len(recs), "N_objects": len({r["source_id"] for r in recs}), **res}

    rem_12 = [r for r in all_pair_records if r["source_id"] != "RUBIES-EGS-42232"]
    s6_base = eval_slicing(rem_12)
    s6 = {
        "dropped_source_id": "RUBIES-EGS-42232",
        "dropped_N_pairs": len(all_pair_records) - len(rem_12),
        "remaining_N_pairs": len(rem_12),
        **s6_base,
    }
    slicings = {
        "1_all_pairs_22": eval_slicing(all_pair_records),
        "2_same_grating_7": eval_slicing([r for r in all_pair_records if r["same_grating"]]),
        "3_cross_grating_15": eval_slicing([r for r in all_pair_records if not r["same_grating"]]),
        "4_top2_one_pair_per_source_9": eval_slicing([r for r in all_pair_records if r["idx_a"] == 0 and r["idx_b"] == 1]),
        "5_equal_object_weight_9": eval_slicing(all_pair_records, weights=[1.0 / (9.0 * r["n_pairs_in_obj"]) for r in all_pair_records]),
        "6_drop_RUBIES_EGS_42232_12": s6,
    }
    assert abs(slicings["1_all_pairs_22"]["sy_floor_dex"] - 0.06010) < 5e-5
    assert abs(slicings["2_same_grating_7"]["sy_floor_dex"] - 0.04449) < 5e-5

    loo_slicings = {}
    for sid_drop in sorted({r["source_id"] for r in all_pair_records}):
        sub = [r for r in all_pair_records if r["source_id"] != sid_drop]
        loo_slicings[sid_drop] = eval_slicing(sub)

    re_43_b1 = ranef_fixed_primary_exact(spec_by_grp, use_cxy=True)
    re_43_b2 = ranef_per_spectrum_residual_2dcov(spec_by_grp)
    re_43_b3 = ranef_fixed_primary_exact(spec_by_grp, use_cxy=False)

    grps_39 = {g: ([v[0]] if v[0]["source_id"] == "RUBIES-EGS-42232" else v) for g, v in spec_by_grp.items()}
    re_39_b1 = ranef_fixed_primary_exact(grps_39, use_cxy=True)
    re_39_b2 = ranef_per_spectrum_residual_2dcov(grps_39)
    re_39_b3 = ranef_fixed_primary_exact(grps_39, use_cxy=False)

    re_loo_primary = {}
    for g_target, v_target in multi_by_grp.items():
        sid_t = v_target[0]["source_id"]
        grps_mod = {g: ([v[0]] if g == g_target else v) for g, v in spec_by_grp.items()}
        r_b1 = ranef_fixed_primary_exact(grps_mod, use_cxy=True)
        re_loo_primary[sid_t] = {
            "N_spectra": r_b1["N_spectra"],
            "f_med": r_b1["f_med"],
            "sigma_med": r_b1["sigma_med"],
            "sigma_p95_onesided": r_b1["sigma_p95_onesided"],
            "slope_med": r_b1["slope_q"][2],
        }

    all_p95 = (
        [v["sigma_p95_onesided"] for v in slicings.values()]
        + [v["sigma_p95_onesided"] for v in loo_slicings.values()]
        + [re_43_b1["sigma_p95_onesided"], re_43_b2["sigma_p95_onesided"], re_43_b3["sigma_p95_onesided"]]
        + [re_39_b1["sigma_p95_onesided"], re_39_b2["sigma_p95_onesided"], re_39_b3["sigma_p95_onesided"]]
        + [v["sigma_p95_onesided"] for v in re_loo_primary.values()]
    )

    with (RESULTS / "repeat_floor_slicings.tsv").open("w", newline="") as fp:
        fields = ["slicing", "N_objects", "N_pairs", "sy_floor_dex", "sigma_med_dex", "sigma_16_dex", "sigma_84_dex", "sigma_p95_onesided_dex", "slope_med"]
        w = csv.DictWriter(fp, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for k, v in slicings.items():
            w.writerow({
                "slicing": k,
                "N_objects": v["N_objects"],
                "N_pairs": v["N_pairs"],
                "sy_floor_dex": f"{v['sy_floor_dex']:.5f}",
                "sigma_med_dex": f"{v['sigma_med']:.4f}",
                "sigma_16_dex": f"{v['sigma_q'][1]:.4f}",
                "sigma_84_dex": f"{v['sigma_q'][3]:.4f}",
                "sigma_p95_onesided_dex": f"{v['sigma_p95_onesided']:.4f}",
                "slope_med": f"{v['slope_q'][2]:.4f}",
            })

    repeat_summary = {
        "pairwise_slicings_6": slicings,
        "pairwise_leave_one_object_out_9": loo_slicings,
        "random_effects_models": {
            "exact_2dcov_fixed_primary_x": {
                "all_43_spectra": re_43_b1,
                "42232_primary_only_39_spectra": re_39_b1,
            },
            "exact_2dcov_per_spectrum_residual": {
                "all_43_spectra": re_43_b2,
                "42232_primary_only_39_spectra": re_39_b2,
            },
            "diagonal_fixed_primary_x_control": {
                "all_43_spectra": re_43_b3,
                "42232_primary_only_39_spectra": re_39_b3,
            },
        },
        "random_effects_reduce_one_object_to_primary_9": re_loo_primary,
        "criterion_evaluation": {
            "min_p95_onesided_dex": float(min(all_p95)),
            "max_p95_onesided_dex": float(max(all_p95)),
            "all_le_0p065": bool(max(all_p95) <= 0.065),
            "in_0p045_to_0p065_band": bool(0.045 <= min(all_p95) and max(all_p95) <= 0.065),
        },
    }
    (RESULTS / "repeat_random_effects.json").write_text(json.dumps(repeat_summary, indent=2) + "\n")

    # --- Part 3: Exterior-envelope gravitational potential & escape-speed check ---
    sxmax = 0.3
    temps = {}
    for tag, fn in [
        ("Sample_A_Curtis_BX26", "atmosphere_bx26_widths.tsv"),
        ("Sample_B_Curtis_own", "atmosphere_our_widths.tsv"),
        ("Sample_C_deGraaff_all", "blackbody_sources.tsv"),
    ]:
        rows = [r for r in csv.DictReader((PHOTOSPHERE / fn).open(), delimiter="\t") if float(r["sx"]) < sxmax]
        temps[tag] = [float(r["T_K"]) for r in rows]
        if tag == "Sample_C_deGraaff_all":
            temps["Sample_E_deGraaff_own"] = [float(r["T_K"]) for r in rows if r["width_source"] in ("frozen", "T2_quality")]
    assert len(temps["Sample_A_Curtis_BX26"]) == 18 and len(temps["Sample_B_Curtis_own"]) == 15
    assert len(temps["Sample_C_deGraaff_all"]) == 9 and len(temps["Sample_E_deGraaff_own"]) == 5

    all_T = temps["Sample_A_Curtis_BX26"] + temps["Sample_B_Curtis_own"] + temps["Sample_C_deGraaff_all"]
    T_min_obs, T_med_obs, T_max_obs = float(min(all_T)), float(np.median(all_T)), float(max(all_T))
    c_5000 = (4.0 * math.pi * G_CGS * C_CGS) / (SIGMA_SB_CGS * (5000.0 ** 4) * KAPPA_ES_MS)
    dm_bench = delta_m_over_medd(5000.0, 1e-3, 100.0)
    dv_bench = delta_vesc_frac(dm_bench, lam=1.0)
    assert abs(dm_bench - 8.2093186e-3) < 1e-8

    kappas = [0.1, 0.02, 8.7e-3, 1e-3, 3e-4, 1e-4]
    T_points = [
        ("T_nominal_min_3600K", 3600.0),
        ("T_sample_min_3637K", T_min_obs),
        ("T_Curtis_AppA_4109K", 4109.0),
        ("T_sample_median_4604K", T_med_obs),
        ("T_fiducial_5000K", 5000.0),
        ("T_sample_max_6381K", T_max_obs),
        ("T_nominal_max_6400K", 6400.0),
    ]
    grid_rows = []
    for k_val in kappas:
        row_T = {}
        for t_label, t_val in T_points:
            dm_asym = delta_m_over_medd(t_val, k_val, 100.0, exact_finite_rout=False)
            dm_fin = delta_m_over_medd(t_val, k_val, 100.0, exact_finite_rout=True)
            row_T[t_label] = {
                "T_K": t_val,
                "DeltaM_over_MEdd": dm_asym,
                "DeltaM_over_MEdd_finite_rout100": dm_fin,
                "Delta_vesc_pct_at_lambda1": 100.0 * delta_vesc_frac(dm_asym, lam=1.0),
                "Delta_vesc_pct_at_lambda1_finite_rout100": 100.0 * delta_vesc_frac(dm_fin, lam=1.0),
                "Delta_vesc_pct_at_lambda0p29": 100.0 * delta_vesc_frac(dm_asym, lam=0.29),
                "relaxed_lambda_ceiling_from_1": 1.0 / (1.0 - dm_asym),
                "relaxed_lambda_ceiling_finite_rout100": 1.0 / (1.0 - dm_fin),
            }
        all_dv_pct = [100.0 * delta_vesc_frac(delta_m_over_medd(tv, k_val, 100.0), lam=1.0) for tv in all_T]
        grid_rows.append({
            "kappa_cm2_g": k_val,
            "by_temperature": row_T,
            "sample_ABC_42fits_dv_pct_p16_p50_p84_max": [
                float(np.percentile(all_dv_pct, 16)),
                float(np.percentile(all_dv_pct, 50)),
                float(np.percentile(all_dv_pct, 84)),
                float(np.max(all_dv_pct)),
            ],
        })

    alpha_checks = {}
    for alpha in (1.5, 2.0, 2.5, 3.0):
        dm_5000 = delta_m_general_powerlaw(5000.0, 1e-3, alpha, 100.0)
        dm_tmin = delta_m_general_powerlaw(T_min_obs, 1e-3, alpha, 100.0)
        alpha_checks[f"alpha_{alpha:.1f}"] = {
            "T5000_DeltaM_over_MEdd": dm_5000,
            "T5000_Delta_vesc_pct": 100.0 * delta_vesc_frac(dm_5000),
            "Tmin_DeltaM_over_MEdd": dm_tmin,
            "Tmin_Delta_vesc_pct": 100.0 * delta_vesc_frac(dm_tmin),
        }

    pot_out = {
        "constants_cgs": {
            "G": G_CGS,
            "c": C_CGS,
            "sigma_SB": SIGMA_SB_CGS,
            "kappa_es_ms": KAPPA_ES_MS,
            "kappa_es_exact": KAPPA_ES_EXACT,
            "prefactor_5000K_cm2_g": c_5000,
        },
        "self_check_benchmark": {
            "T_eff_K": 5000.0,
            "kappa_cm2_g": 1e-3,
            "r_out_over_Rph": 100.0,
            "DeltaM_over_MEdd": dm_bench,
            "Delta_vesc_pct": 100.0 * dv_bench,
        },
        "sample_temperature_summary_K": {
            "N_total_kept_fits_ABC": len(all_T),
            "T_min_obs_K": T_min_obs,
            "T_median_obs_K": T_med_obs,
            "T_max_obs_K": T_max_obs,
        },
        "opacity_temperature_grid": grid_rows,
        "density_slope_sensitivity_kappa1e3_rout100": alpha_checks,
    }
    (RESULTS / "exterior_potential_check.json").write_text(json.dumps(pot_out, indent=2) + "\n")


if __name__ == "__main__":
    main()
