"""Sensitivity of broad-line widths to electron-scattering optical depth for Gaussian and top-hat cores convolved
with the exponential scattering kernel of Rusakov et al. (2026), with and without an unscattered core component.

Units: velocities in km s^-1, optical depth tau_e dimensionless, scatters in dex (base-10 logarithms).
Writes data/results/scattering_kernel.json.
"""

import json
import math

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erfcx

from lrdwidth.paths import RESULTS

RUS_A = 428.0  # km s^-1 per unit tau_e at T_e = 10^4 K (Rusakov et al. 2026)
RUS_B = 370.0  # km s^-1 intercept at T_e = 10^4 K (Rusakov et al. 2026)
FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def W_rusakov(tau_e: float, T4: float = 1.0) -> float:
    """Wing e-folding scale W in km s^-1 (Rusakov et al. 2026)."""
    return (RUS_A * tau_e + RUS_B) * math.sqrt(T4)


def D_wing(tau_e: float) -> float:
    """Analytical logarithmic derivative d log W / d log tau_e."""
    return (RUS_A * tau_e) / (RUS_A * tau_e + RUS_B)


def gauss_profile(v: np.ndarray, fwhm_in: float) -> np.ndarray:
    """Unit-area Gaussian profile with given FWHM (km s^-1)."""
    sig = fwhm_in * FWHM_TO_SIGMA
    return np.exp(-0.5 * (v / sig) ** 2) / (math.sqrt(2.0 * math.pi) * sig)


def gauss_exp_conv_exact(v: np.ndarray, fwhm_in: float, W: float) -> np.ndarray:
    """Exact analytical convolution of a unit-area Gaussian with E(u) = exp(-|u|/W) / (2W)."""
    sig = fwhm_in * FWHM_TO_SIGMA
    av = np.abs(v)
    inv_rt2 = 1.0 / math.sqrt(2.0)
    zp = (sig / W + av / sig) * inv_rt2
    zm = (sig / W - av / sig) * inv_rt2
    base = np.exp(-0.5 * (av / sig) ** 2)
    term_p = base * erfcx(zp)
    term_m = np.where(
        zm >= 0.0,
        base * erfcx(np.maximum(zm, 0.0)),
        2.0 * np.exp(0.5 * (sig / W) ** 2 - av / W) - base * erfcx(np.maximum(-zm, 0.0)),
    )
    return (0.25 / W) * (term_p + term_m)


def tophat_exp_conv_exact(v: np.ndarray, v_term: float, W: float) -> np.ndarray:
    """Exact analytical convolution of a unit-area top-hat on [-v_term, v_term] with E(u) = exp(-|u|/W) / (2W)
    (Liu, Feng & Ho 2026, Eqs. 7-8)."""
    av = np.abs(v)
    out = np.empty_like(av, dtype=float)
    inside = av <= v_term
    out[inside] = (1.0 - np.exp(-v_term / W) * np.cosh(av[inside] / W)) / (2.0 * v_term)
    out[~inside] = (np.exp(-av[~inside] / W) * math.sinh(v_term / W)) / (2.0 * v_term)
    return out


def measure_direct_fwhm(v_pos: np.ndarray, prof_pos: np.ndarray) -> float:
    """Direct FWHM on a symmetric monotone-decreasing non-negative profile defined on v_pos >= 0."""
    half_peak = 0.5 * prof_pos[0]
    idx = int(np.searchsorted(-prof_pos, -half_peak))
    assert 0 < idx < len(v_pos), f"Half-maximum bracket failed: idx={idx}"
    v1, v2 = float(v_pos[idx - 1]), float(v_pos[idx])
    y1, y2 = float(prof_pos[idx - 1] - half_peak), float(prof_pos[idx] - half_peak)
    v_half = v1 - y1 * (v2 - v1) / (y2 - y1)
    return 2.0 * v_half


def measure_w50(v_pos: np.ndarray, prof_pos: np.ndarray) -> float:
    """Interquartile velocity width W_50 = v_75 - v_25."""
    dv = float(v_pos[1] - v_pos[0])
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (prof_pos[:-1] + prof_pos[1:]) * dv)))
    total_half = cum[-1]
    assert abs(total_half - 0.5) < 1e-4, f"Half-integral {total_half} != 0.5"
    target = 0.5 * total_half
    idx = int(np.searchsorted(cum, target))
    v1, v2 = float(v_pos[idx - 1]), float(v_pos[idx])
    c1, c2 = float(cum[idx - 1]), float(cum[idx])
    v_q = v1 + (target - c1) * (v2 - v1) / (c2 - c1)
    return 2.0 * v_q


def measure_single_gaussian_fwhm(v_pos: np.ndarray, prof_pos: np.ndarray) -> float:
    """Unweighted least-squares fit of a single Gaussian over [-V_MAX, +V_MAX]."""
    w_sqrt = np.sqrt(np.where(v_pos == 0.0, 1.0, 2.0))
    w50_init = measure_w50(v_pos, prof_pos)
    sig0 = max(w50_init / 1.3489795, 50.0)
    amp0 = float(prof_pos[0])

    def resid(p):
        amp, sig = p
        model = amp * np.exp(-0.5 * (v_pos / sig) ** 2)
        return w_sqrt * (model - prof_pos)

    res = least_squares(
        resid,
        x0=[amp0, sig0],
        bounds=([0.0, 1.0], [np.inf, 5e4]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    assert res.success, f"Gaussian fit failed: {res.message}"
    return float(res.x[1] / FWHM_TO_SIGMA)


def eval_widths(tau_e: float, mode: str, fwhm_in: float, v_pos: np.ndarray) -> dict:
    W = W_rusakov(tau_e)
    if mode == "pure_conv_gauss":
        prof = gauss_exp_conv_exact(v_pos, fwhm_in, W)
    elif mode == "pure_conv_tophat":
        prof = tophat_exp_conv_exact(v_pos, 0.5 * fwhm_in, W)
    elif mode == "rt_unscat_plus_scat_gauss":
        f0 = gauss_profile(v_pos, fwhm_in)
        fc = gauss_exp_conv_exact(v_pos, fwhm_in, W)
        et = math.exp(-tau_e)
        prof = et * f0 + (1.0 - et) * fc
    else:
        raise ValueError(mode)
    return {
        "gauss_fit_fwhm": measure_single_gaussian_fwhm(v_pos, prof),
        "direct_fwhm": measure_direct_fwhm(v_pos, prof),
        "w50": measure_w50(v_pos, prof),
    }


def log_deriv(
    tau_e: float, mode: str, fwhm_in: float, v_pos: np.ndarray, rel_step: float = 1e-3
) -> dict:
    """Centered finite-difference d log(width) / d log(tau_e)."""
    t_hi = tau_e * math.exp(rel_step)
    t_lo = tau_e * math.exp(-rel_step)
    w_mid = eval_widths(tau_e, mode, fwhm_in, v_pos)
    w_hi = eval_widths(t_hi, mode, fwhm_in, v_pos)
    w_lo = eval_widths(t_lo, mode, fwhm_in, v_pos)
    out = {}
    for k in ("gauss_fit_fwhm", "direct_fwhm", "w50"):
        d_eff = math.log(w_hi[k] / w_lo[k]) / (2.0 * rel_step)
        out[k] = {
            "width_kms": w_mid[k],
            "D_eff": d_eff,
            "ratio_to_D_wing": d_eff / D_wing(tau_e),
        }
    return out


def main():
    v_full = np.linspace(-25000.0, 25000.0, 200001)
    v_pos = v_full[v_full >= 0.0]

    # Self-check 1: direct numerical convolution vs exact analytical formulas
    v_small = np.linspace(-15000.0, 15000.0, 60001)
    dv_s = float(v_small[1] - v_small[0])
    W_test = W_rusakov(1.5)
    kern = np.exp(-np.abs(v_small) / W_test) / (2.0 * W_test)
    g0 = gauss_profile(v_small, 800.0)
    num_conv_g = np.convolve(g0, kern, mode="same") * dv_s
    ana_conv_g = gauss_exp_conv_exact(v_small, 800.0, W_test)
    mask_mid = np.abs(v_small) < 8000.0
    assert np.max(np.abs(num_conv_g[mask_mid] - ana_conv_g[mask_mid])) < 1e-7

    th0 = np.where(np.abs(v_small) <= 600.0, 1.0 / 1200.0, 0.0)
    num_conv_th = np.convolve(th0, kern, mode="same") * dv_s
    ana_conv_th = tophat_exp_conv_exact(v_small, 600.0, W_test)
    assert np.max(np.abs(num_conv_th[mask_mid] - ana_conv_th[mask_mid])) < 5e-7

    # Self-check 2: Liu, Feng & Ho (2026) Table 1 row 1 (v_term=622, sigma_v=848 -> FWHM_mod=1910 km/s)
    prof_egg = tophat_exp_conv_exact(v_pos, 622.0, 848.0)
    fwhm_egg = measure_direct_fwhm(v_pos, prof_egg)
    assert abs(fwhm_egg - 1910.0) < 1.5

    # Self-check 3: narrow-core limit recovers D_wing
    d_lim = log_deriv(1.0, "pure_conv_gauss", 0.05, v_pos)
    for k in ("gauss_fit_fwhm", "direct_fwhm", "w50"):
        assert abs(d_lim[k]["D_eff"] - D_wing(1.0)) < 5e-5

    r1 = json.load(open(RESULTS / "width_luminosity.json"))
    r6 = json.load(open(RESULTS / "photosphere_constrained.json"))
    r11 = json.load(open(RESULTS / "photosphere_simulations.json"))
    s0 = json.load(open(RESULTS / "repeat_calibration.json"))["posterior_with_floors"]["formal_zero_floor"]["LRD_p95_onesided"]
    rats = {
        "PA_Curtis_BX": r6["A"]["ratio_median"],
        "PB_Curtis_own": r6["B"]["ratio_median"],
        "PE_deGraaff_own": r11["main_E"]["ratio_p16_50_84"][1],
    }
    grav_max_frac = 1.0 / min(rats.values())
    assert round(100.0 * grav_max_frac) == 40

    taus = [0.5, 1.0, 2.0, 3.0]
    fwhm_ins = [500.0, 750.0, 1000.0]
    modes = ["pure_conv_gauss", "pure_conv_tophat", "rt_unscat_plus_scat_gauss"]

    grid = {}
    for mode in modes:
        grid[mode] = {}
        for fin in fwhm_ins:
            key_fin = f"fwhm_in_{int(fin)}"
            grid[mode][key_fin] = {}
            for t in taus:
                res = log_deriv(t, mode, fin, v_pos)
                d_w = D_wing(t)
                for m_v in res.values():
                    m_v["sigma_log_tau_budget_dex"] = s0 / m_v["D_eff"]
                    m_v["sigma_log_R_budget_dex"] = 2.0 * s0 / m_v["D_eff"]
                grid[mode][key_fin][f"tau_{t:g}"] = {
                    "tau_e": t,
                    "W_wing_kms": W_rusakov(t),
                    "D_wing": d_w,
                    "sigma_log_tau_wing_dex": s0 / d_w,
                    "sigma_log_R_wing_dex": 2.0 * s0 / d_w,
                    "measures": res,
                }

    fins_dense = np.linspace(500.0, 1000.0, 21)
    taus_23 = np.linspace(2.0, 3.0, 21)

    conv_half_D = [
        log_deriv(0.5, "pure_conv_gauss", float(fi), v_pos)["gauss_fit_fwhm"]["D_eff"]
        for fi in fins_dense
    ]
    conv_23_D = [
        log_deriv(float(ta), "pure_conv_gauss", float(fi), v_pos)["gauss_fit_fwhm"]["D_eff"]
        for fi in fins_dense
        for ta in taus_23
    ]
    rt_half_boost = [
        eval_widths(0.5, "rt_unscat_plus_scat_gauss", float(fi), v_pos)["gauss_fit_fwhm"]
        / float(fi)
        - 1.0
        for fi in fins_dense
    ]
    rt_23_D = [
        log_deriv(float(ta), "rt_unscat_plus_scat_gauss", float(fi), v_pos)["gauss_fit_fwhm"][
            "D_eff"
        ]
        for fi in fins_dense
        for ta in taus_23
    ]

    summary_ranges = {
        "fwhm_in_range_kms": [500, 1000],
        "pure_conv_gauss_tau0p5": {
            "D_eff_min": min(conv_half_D),
            "D_eff_max": max(conv_half_D),
            "sigma_log_tau_min_dex": s0 / max(conv_half_D),
            "sigma_log_tau_max_dex": s0 / min(conv_half_D),
        },
        "pure_conv_gauss_tau2to3": {
            "D_eff_min": min(conv_23_D),
            "D_eff_max": max(conv_23_D),
            "sigma_log_tau_min_dex": s0 / max(conv_23_D),
            "sigma_log_tau_max_dex": s0 / min(conv_23_D),
            "sigma_log_R_min_dex": 2.0 * s0 / max(conv_23_D),
            "sigma_log_R_max_dex": 2.0 * s0 / min(conv_23_D),
        },
        "rt_unscat_plus_scat_gauss_tau0p5": {
            "max_fractional_boost_over_core": max(rt_half_boost),
            "min_fractional_boost_over_core": min(rt_half_boost),
        },
        "rt_unscat_plus_scat_gauss_tau2to3": {
            "D_eff_min": min(rt_23_D),
            "D_eff_max": max(rt_23_D),
            "sigma_log_tau_min_dex": s0 / max(rt_23_D),
            "sigma_log_tau_max_dex": s0 / min(rt_23_D),
            "sigma_log_R_min_dex": 2.0 * s0 / max(rt_23_D),
            "sigma_log_R_max_dex": 2.0 * s0 / min(rt_23_D),
        },
    }
    assert summary_ranges["rt_unscat_plus_scat_gauss_tau2to3"]["D_eff_min"] > D_wing(3.0)

    out = {
        "super_eddington_gravity_share": {
            "median_FWHM_over_vK": rats,
            "min_median_ratio": min(rats.values()),
            "max_median_ratio": max(rats.values()),
            "gravity_share_max_at_lambda1_k1": grav_max_frac,
            "gravity_share_max_pct_rounded": int(round(100.0 * grav_max_frac)),
        },
        "scattering_sensitivity_grid": grid,
        "extrema_summary": summary_ranges,
    }
    with open(RESULTS / "scattering_kernel.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
