"""Regressions of line width on line luminosity and on photospheric parameters."""

import csv
import json
import math

import numpy as np

from . import stats
from .stats import diffq, hpd_level, quant  # noqa: F401  (re-exported)


def post2(x, y, cxx, cyy, cxy, **grid):
    """Posterior of slope and intrinsic scatter for y = a + b x; see ``stats.post2``. Returns the two grids and the posterior."""
    return stats.post2(x, y, cxx, cyy, cxy, **grid)[:3]


def post3g(x1, x2, y, s1, s2, sy, lim1=(-1, 1), lim2=(-1, 1), n1=161, n2=161, ns=241):
    """Two-predictor posterior on a grid; see ``stats.post3``."""
    return stats.post3(x1, x2, y, s1, s2, sy, lim1, lim2, n1, n2, ns)


def ols_sd(x, y):
    """Standard deviation of the residuals about an ordinary least-squares line, and its slope."""
    c = np.polyfit(x, y, 1)
    r = y - np.polyval(c, x)
    return (float(np.std(r, ddof=1)), float(c[0]))


def post_lum_temp(x, t, y, sx, st, sy, nb=121, nc=181, ns=241):
    """y = a + b x + c t + N(0, s^2) with independent errors on x, t and y; flat prior on a, integrated analytically.

    Returns the grids of b, c and s, the joint posterior, and the conditional mean of a."""
    b = -1 + (np.arange(nb) + 0.5) * 2 / nb
    cc = -3 + (np.arange(nc) + 0.5) * 6 / nc
    s = (np.arange(ns) + 0.5) * 0.6 / ns
    B = b[:, None, None]
    Cg = cc[None, :, None]
    S = s[None, None, :]
    sw = 0
    swr = 0
    swr2 = 0
    slv = 0
    for i in range(len(x)):
        V = S * S + sy[i] ** 2 + B * B * sx[i] ** 2 + Cg * Cg * st[i] ** 2
        R = y[i] - B * x[i] - Cg * t[i]
        w = 1 / V
        sw = sw + w
        swr = swr + w * R
        swr2 = swr2 + w * R * R
        slv = slv + np.log(V)
    lp = -0.5 * (slv + swr2 - swr * swr / sw) - 0.5 * np.log(sw)
    ah = swr / sw
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    return (b, cc, s, p, ah)


def fit(rows):
    """Width-luminosity fit of one sample: posterior quantiles of scatter and slope, the posterior-mean intercept, OLS summary."""
    x = np.array([r["x"] for r in rows])
    y = np.array([r["y"] for r in rows])
    C = np.array([r["cov"] for r in rows])
    b, s, p = post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
    B = b[:, None]
    S = s[None, :]
    sw = 0
    swr = 0
    for xi, yi, ci in zip(x, y, C):
        V = S * S + ci[1, 1] + B * B * ci[0, 0] - 2 * B * ci[0, 1]
        sw = sw + 1 / V
        swr = swr + (yi - B * xi) / V
    a = float((p * swr / sw).sum())
    sd, sl = ols_sd(x, y)
    return (
        dict(
            N=len(rows),
            sigma=quant(s, p.sum(0)),
            slope=quant(b, p.sum(1)),
            a=a,
            ols_sd=sd,
            ols_slope=sl,
        ),
        (s, p.sum(0)),
    )


def pair(L, O):
    """Fit LRDs and comparison sources and return both, with the quantiles of the difference in intrinsic scatter."""
    rl, ml = fit(L)
    ro, mo = fit(O)
    return dict(LRD=rl, other=ro, diff=diffq(ml[0], ml[1], mo[1]))


def sig1d(resid, var, n=3000, smax=0.6):
    """Intrinsic scatter of residuals about a fixed relation (only a constant and the scatter are fitted)."""
    s = (np.arange(n) + 0.5) * smax / n
    V = s[None, :] ** 2 + var[:, None]
    w = 1 / V
    sw = w.sum(0)
    swr = (w * resid[:, None]).sum(0)
    swr2 = (w * resid[:, None] ** 2).sum(0)
    lp = -0.5 * (np.log(V).sum(0) + swr2 - swr**2 / sw) - 0.5 * np.log(sw)
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    return quant(s, p)


def ori(i1, i2, h, n=200001):
    """Scatter (dex) of the projected speed for v^2 proportional to sin^2 i + h^2, with cos i uniform between i1 and i2 (degrees)."""
    mu = np.linspace(math.cos(math.radians(i2)), math.cos(math.radians(i1)), n)
    v = 0.5 * np.log10(1 - mu**2 + h * h + 1e-300)
    v = v[np.isfinite(v)]
    return float(np.std(v))


def load_line_inputs(data):
    """Line luminosities and widths of our sample and of the BX26 table, indexed by source identifier."""
    ours = json.load(open(data / "widths" / "regression_inputs.json"))["single_gaussian"]
    bx = json.load(open(data / "bx26" / "regression_inputs.json"))
    return (ours, {r["source_id"]: r for r in ours}, bx, {t["id"]: t for t in bx})


def load_photosphere_sample(path, tag, byf, byx):
    """One photospheric sample with the line luminosity of each source (tag "A": BX26 widths; otherwise our widths)."""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    xph = np.array([float(r["x"]) for r in rows])
    t = np.array([float(r["t"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    sph = np.array([float(r["sx"]) for r in rows])
    st = np.array([float(r["st"]) for r in rows])
    sy = np.array([float(r["sy"]) for r in rows])
    if tag == "A":
        xl = np.array([byx[r["width_id"]]["logLbol"] - 45 for r in rows])
        sl = np.array([byx[r["width_id"]]["sx"] for r in rows])
    else:
        xl = np.array([byf[r["width_id"]]["x"] for r in rows])
        sl = np.array([math.sqrt(byf[r["width_id"]]["cov"][0][0]) for r in rows])
    return (rows, xl, sl, xph, sph, t, st, y, sy)
