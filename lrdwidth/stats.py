"""Grid posteriors and summary statistics.

All posteriors are evaluated on regular grids of cell centres. Quantiles are read from the cumulative distribution
evaluated at the cell edges (see ``quant``).
"""

import math

import numpy as np

erf = np.vectorize(math.erf)


def ndtr(z):
    return 0.5 * (1 + erf(z / math.sqrt(2)))


def quant(grid, pm, qs=(0.025, 0.16, 0.5, 0.84, 0.975)):
    c = np.cumsum(pm) / pm.sum()
    h = grid[1] - grid[0]
    return [
        float(
            np.interp(
                v, np.concatenate([[0], c]), np.concatenate([[grid[0] - h / 2], grid + h / 2])
            )
        )
        for v in qs
    ]


def post2(x, y, cxx, cyy, cxy, b_lo=-1, b_hi=1, nb=801, s_lo=0, s_hi=0.6, ns=601, a_lo=2, a_hi=4.5):
    """y = a + b x + N(0, s^2); per-source covariance; a uniform on [a_lo, a_hi], integrated analytically."""
    b = b_lo + (np.arange(nb) + 0.5) * (b_hi - b_lo) / nb
    s = s_lo + (np.arange(ns) + 0.5) * (s_hi - s_lo) / ns
    B = b[:, None]
    S = s[None, :]
    sw = 0
    swr = 0
    swr2 = 0
    slv = 0
    for xi, yi, a, c, d in zip(x, y, cxx, cyy, cxy):
        V = S * S + c + B * B * a - 2 * B * d
        R = yi - B * xi
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
    return (b, s, p, ah)


def diffq(s, p1, p2):
    n = len(s)
    d = np.convolve(p1, p2[::-1])
    g = (np.arange(2 * n - 1) - (n - 1)) * (s[1] - s[0])
    return quant(g, d)


def post3(x1, x2, y, s1, s2, sy, lim1=(-1, 1), lim2=(-3, 3), n1=161, n2=241, ns=301, smax=0.6):
    """y = a + b1 x1 + b2 x2 + N(0, s^2); independent errors on all three axes; flat prior on a, integrated analytically."""
    b1 = lim1[0] + (np.arange(n1) + 0.5) * (lim1[1] - lim1[0]) / n1
    b2 = lim2[0] + (np.arange(n2) + 0.5) * (lim2[1] - lim2[0]) / n2
    s = (np.arange(ns) + 0.5) * smax / ns
    B1 = b1[:, None, None]
    B2 = b2[None, :, None]
    S = s[None, None, :]
    sw = 0
    swr = 0
    swr2 = 0
    slv = 0
    for i in range(len(y)):
        V = S * S + sy[i] ** 2 + B1 * B1 * s1[i] ** 2 + B2 * B2 * s2[i] ** 2
        R = y[i] - B1 * x1[i] - B2 * x2[i]
        w = 1 / V
        sw = sw + w
        swr = swr + w * R
        swr2 = swr2 + w * R * R
        slv = slv + np.log(V)
    lp = -0.5 * (slv + swr2 - swr * swr / sw) - 0.5 * np.log(sw)
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    return (b1, b2, s, p)


def hpd_level(g1, g2, p12, v1, v2):
    """Posterior mass inside the iso-density contour through (v1, v2)."""
    i = np.interp(v1, g1, np.arange(len(g1)))
    j = np.interp(v2, g2, np.arange(len(g2)))
    i0 = int(i)
    j0 = int(j)
    f = i - i0
    g = j - j0
    d = (
        p12[i0, j0] * (1 - f) * (1 - g)
        + p12[i0 + 1, j0] * f * (1 - g)
        + p12[i0, j0 + 1] * (1 - f) * g
        + p12[i0 + 1, j0 + 1] * f * g
    )
    return float(p12[p12 > d].sum())


def contour_levels(p12, masses=(0.68, 0.95)):
    v = np.sort(p12.ravel())[::-1]
    c = np.cumsum(v)
    return [float(v[np.searchsorted(c, m)]) for m in masses]


# Keplerian speed at the photosphere of a source radiating at lambda times the Eddington luminosity of the enclosed mass
G = 6.6743e-08
CLIGHT = 29979245800.0
SIGSB = 5.670374419e-05
KAPPA = 6.6524587321e-25 / 1.67262192369e-24
VK_COEF = math.sqrt(KAPPA / (4 * math.pi * CLIGHT)) * (4 * math.pi * SIGSB) ** 0.25
VK0_KMS = VK_COEF * 1e45**0.25 * 5000 / 100000.0
