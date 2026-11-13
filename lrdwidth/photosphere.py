"""Width against line luminosity and photospheric luminosity, with simulations under the two hypotheses.

Models (flat priors on the true predictors, intercept integrated analytically):
  partial regression   y = a + b1 x_line + b2 x_ph + N(0, s^2),
                       V_i = s^2 + C_yy + b1^2 C_ll + b2^2 s_x^2 - 2 b1 C_ly
  fixed scaling        y - 0.25 x_ph - t = a + N(0, s^2),  V_i = s^2 + C_yy + (0.25 s_x)^2 + s_t^2
Priors: a ~ U(2, 4.5), b ~ U(-1, 1), s ~ U(0, 0.6).
"""

import json
import math

import numpy as np

from .stats import ndtr, quant

G = 6.6743e-08
CL = 29979245800.0
SB = 5.670374419e-05
KAPPA = 6.6524587321e-25 / 1.67262192369e-24


def log_vk(logL, T):
    """log10 of the Keplerian speed (km/s) at the photosphere of a source at its Eddington luminosity."""
    L = 10**logL
    M = L * KAPPA / (4 * math.pi * G * CL)
    R = math.sqrt(L / (4 * math.pi * SB * T**4))
    return math.log10(math.sqrt(G * M / R) / 100000.0)


def partial(xl, xp, y, cll, sx, cyy, cly, nb1=101, nb2=101, ns=150, a_lo=2.0, a_hi=4.5):
    """y = a + b1 x_line + b2 x_ph + N(0, s^2) with correlated errors on (x_line, y); a uniform on [a_lo, a_hi] or flat."""
    b1 = -1 + (np.arange(nb1) + 0.5) * 2 / nb1
    b2 = -1 + (np.arange(nb2) + 0.5) * 2 / nb2
    s = (np.arange(ns) + 0.5) * 0.6 / ns
    B1 = b1[:, None, None]
    B2 = b2[None, :, None]
    S = s[None, None, :]
    sw = 0
    swr = 0
    swr2 = 0
    slv = 0
    for i in range(len(y)):
        V = S * S + cyy[i] + B1 * B1 * cll[i] + B2 * B2 * sx[i] ** 2 - 2 * B1 * cly[i]
        R = y[i] - B1 * xl[i] - B2 * xp[i]
        w = 1 / V
        sw = sw + w
        swr = swr + w * R
        swr2 = swr2 + w * R * R
        slv = slv + np.log(V)
    ah = swr / sw
    asd = 1 / np.sqrt(sw)
    norm = ndtr((a_hi - ah) / asd) - ndtr((a_lo - ah) / asd) if a_lo is not None else 1.0
    lp = -0.5 * (slv + swr2 - swr * swr / sw) - 0.5 * np.log(sw) + np.log(np.maximum(norm, 1e-300))
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    return (b1, b2, s, p)


def partial_fast(xl, xp, y, cll, sx, cyy, cly, nb1=61, nb2=81, ns=100):
    """Median of b2 on a coarser grid with a flat intercept prior; used inside the simulations."""
    b1, b2, s, p = partial(xl, xp, y, cll, sx, cyy, cly, nb1, nb2, ns, a_lo=None)
    return quant(b2, p.sum((0, 2)), (0.5,))[0]


def one_dim(r, v, ns=2400, a_lo=None, a_hi=None):
    """r_i = a + N(0, s^2 + v_i): posterior quantiles of s and the conditional estimate of a."""
    s = (np.arange(ns) + 0.5) * 0.6 / ns
    S = s[:, None]
    V = S * S + v[None, :]
    w = 1 / V
    sw = w.sum(1)
    swr = (w * r).sum(1)
    swr2 = (w * r * r).sum(1)
    lp = -0.5 * (np.log(V).sum(1) + swr2 - swr * swr / sw) - 0.5 * np.log(sw)
    if a_lo is not None:
        ah = swr / sw
        asd = 1 / np.sqrt(sw)
        lp = lp + np.log(np.maximum(ndtr((a_hi - ah) / asd) - ndtr((a_lo - ah) / asd), 1e-300))
    p = np.exp(lp - lp.max())
    p /= p.sum()
    q = quant(s, p)
    i = int(np.argmin(abs(s - q[2])))
    return (q, float(swr[i] / sw[i]))


def line_only(xl, y, cll, cyy, cly, nb=801, ns=1200, a_lo=2.0, a_hi=4.5):
    """y = a + b x_line + N(0, s^2) with correlated errors: quantiles of b and s, and the intercept at their medians."""
    b = -1 + (np.arange(nb) + 0.5) * 2 / nb
    s = (np.arange(ns) + 0.5) * 0.6 / ns
    B = b[:, None]
    S = s[None, :]
    sw = 0
    swr = 0
    swr2 = 0
    slv = 0
    for i in range(len(y)):
        V = S * S + cyy[i] + B * B * cll[i] - 2 * B * cly[i]
        R = y[i] - B * xl[i]
        w = 1 / V
        sw = sw + w
        swr = swr + w * R
        swr2 = swr2 + w * R * R
        slv = slv + np.log(V)
    ah = swr / sw
    asd = 1 / np.sqrt(sw)
    norm = ndtr((a_hi - ah) / asd) - ndtr((a_lo - ah) / asd)
    lp = -0.5 * (slv + swr2 - swr * swr / sw) - 0.5 * np.log(sw) + np.log(np.maximum(norm, 1e-300))
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    qb = quant(b, p.sum(1))
    qs = quant(s, p.sum(0))
    i = int(np.argmin(abs(b - qb[2])))
    j = int(np.argmin(abs(s - qs[2])))
    return (qb, qs, float(ah[i, j]))


def load_main(path):
    """Main sample: build the regression variables from the per-source table."""
    J = json.load(open(path))
    rows = []
    for s in J:
        w = s["widths"]["single_gaussian"]
        c = w["cov_log_luminosity_width"]
        T16, T50, T84 = s["T_p16_p50_p84_K"]
        L16, L50, L84 = s["logL_MBB_p16_p50_p84"]
        rows.append(
            dict(
                id=s["source_id"],
                srcid=s["dg_srcid"],
                xl=w["logL_demagnified_Planck18_erg_s"] - 43,
                cll=c[0][0],
                cly=c[0][1],
                cyy=c[1][1],
                y=math.log10(w["width_kms"]),
                xp=L50 - 45,
                sx=(L84 - L16) / 2,
                t=math.log10(T50 / 5000),
                st=(math.log10(T84) - math.log10(T16)) / 2,
                logvk=log_vk(L50, T50),
            )
        )
    return rows


def load_sample(path, key):
    """A photospheric sample from the stored per-source regression inputs."""
    J = json.load(open(path))[key]
    rows = []
    for r in J:
        rows.append(
            dict(
                id=r["temperature_id"],
                xl=r["line_x"],
                cll=r["line_sx"] ** 2,
                cly=r.get("line_width_cov") or 0.0,
                cyy=r["sy"] ** 2,
                y=r["y"],
                xp=r["x"],
                sx=r["sx"],
                t=r["t"],
                st=r["st"],
                logvk=r["y"] - r["q_dex"],
            )
        )
    return rows


def arr(rows, k):
    return np.array([r[k] for r in rows], float)


def analyse(rows, nsim, seed, base=None):
    """Observed regressions for one sample, and simulations under the photospheric and the line-luminosity hypotheses."""
    xl, cll, cly, cyy, y, xp, sx, t, st, lvk = [
        arr(rows, k) for k in ("xl", "cll", "cly", "cyy", "y", "xp", "sx", "t", "st", "logvk")
    ]
    q = y - lvk
    out = dict(
        N=len(y),
        ratio_p16_50_84=[float(10**v) for v in np.quantile(q, [0.16, 0.5, 0.84])],
        q_median=float(np.median(q)),
        q_sd=float(np.std(q, ddof=1)),
    )
    qf, a_ph = one_dim(y - 0.25 * xp - t, cyy + (0.25 * sx) ** 2 + st**2, a_lo=2.0, a_hi=4.5)
    out["sigma_fixed_photosphere"] = qf
    qb, qs, a_line = line_only(xl, y, cll, cyy, cly)
    out["line_only"] = dict(b=qb, sigma=qs)
    b1, b2, s, p = partial(xl, xp, y, cll, sx, cyy, cly)
    out["partial"] = dict(
        b1=quant(b1, p.sum((1, 2))), b2=quant(b2, p.sum((0, 2))), sigma=quant(s, p.sum((0, 1)))
    )
    # effect of the bounded intercept prior, and of the coarser grid used in the simulations
    b1u, b2u, su, pu = partial(xl, xp, y, cll, sx, cyy, cly, a_lo=None)
    out["partial"]["b2_median_unbounded_intercept"] = quant(b2u, pu.sum((0, 2)), (0.5,))[0]
    ob = partial_fast(xl, xp, y, cll, sx, cyy, cly)
    out["partial"]["b2_median_sim_grid"] = ob
    # generating parameters: the flat-prior medians of this fit and, if given, the hierarchical baseline
    gens = {
        "own_flat_baseline": dict(a_ph=a_ph, s_ph=qf[2], a_line=a_line, b_line=qb[2], s_line=qs[2])
    }
    if base:
        gens["hierarchical_baseline"] = dict(
            a_ph=base["ph"]["a"],
            s_ph=base["ph"]["sigma"],
            a_line=base["line"]["a"],
            b_line=base["line"]["b_line"],
            s_line=base["line"]["sigma"],
        )
    out["sim"] = {}
    n = len(y)
    for gname, g in gens.items():
        rng = np.random.default_rng([20261113, seed, len(gname)])
        recP = np.empty(nsim)
        recL = np.empty(nsim)
        for it in range(nsim):
            # errors on line luminosity and width are correlated: draw them from the 2x2 covariance of each source
            z1 = rng.normal(0, 1, n)
            z2 = rng.normal(0, 1, n)
            sl_ = np.sqrt(cll)
            rho = cly / np.sqrt(cll * cyy)
            el = sl_ * z1
            ey = np.sqrt(cyy) * (rho * z1 + np.sqrt(1 - rho**2) * z2)
            xlo = xl + el
            xpo = xp + rng.normal(0, 1, n) * sx
            yP = g["a_ph"] + 0.25 * xp + t + rng.normal(0, g["s_ph"], n) + ey
            yL = g["a_line"] + g["b_line"] * xl + rng.normal(0, g["s_line"], n) + ey
            recP[it] = partial_fast(xlo, xpo, yP, cll, sx, cyy, cly)
            recL[it] = partial_fast(xlo, xpo, yL, cll, sx, cyy, cly)
        crit = float(np.quantile(recL, 0.95))
        out["sim"][gname] = dict(
            gen=g,
            N_sim=nsim,
            observed_b2=ob,
            Hph_p16_50_84=[float(v) for v in np.quantile(recP, [0.16, 0.5, 0.84])],
            Hline_p16_50_84=[float(v) for v in np.quantile(recL, [0.16, 0.5, 0.84])],
            n_ph_le_obs=int((recP <= ob).sum()),
            n_line_ge_obs=int((recL >= ob).sum()),
            p_ph=float((recP <= ob).mean()),
            p_line=float((recL >= ob).mean()),
            power=float((recP > crit).mean()),
        )
    return out
