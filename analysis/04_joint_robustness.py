"""Robustness of the joint regression: leave-one-out, the partial slope on temperature, a three-predictor fit, and the
slope recovered from simulated photospheric widths.

Writes data/results/joint_robustness.json.
"""

import json

import numpy as np

from lrdwidth.paths import DATA, PHOTOSPHERE, RESULTS
from lrdwidth.regression import load_line_inputs, load_photosphere_sample, post3g, quant

W1, byf, Q, byx = load_line_inputs(DATA)
out = {}
for tag, fn in (("A", "atmosphere_bx26_widths.tsv"), ("B", "atmosphere_our_widths.tsv")):
    rows, xl, sl, xph, sph, t, st, y, sy = load_photosphere_sample(PHOTOSPHERE / fn, tag, byf, byx)
    # leave one out
    meds = []
    his = []
    los = []
    for i in range(len(y)):
        k = np.arange(len(y)) != i
        b1, b2, s, p = post3g(xl[k], xph[k], y[k], sl[k], sph[k], sy[k], n1=81, n2=81, ns=121)
        q = quant(b2, p.sum((0, 2)))
        meds.append(q[2])
        los.append(q[0])
        his.append(q[4])
    # partial slope on temperature
    b1, c, s, p = post3g(xl, t, y, sl, st, sy, lim2=(-3, 3), n2=241)
    qc = quant(c, p.sum((0, 2)))
    # three predictors: line luminosity, L_ph, T (coarse grid)
    n1 = 61
    n2 = 61
    n3 = 81
    ns = 81
    g1 = -0.3 + (np.arange(n1) + 0.5) * 1.0 / n1
    g2 = -0.5 + (np.arange(n2) + 0.5) * 1.0 / n2
    g3 = -3 + (np.arange(n3) + 0.5) * 6 / n3
    gs = (np.arange(ns) + 0.5) * 0.3 / ns
    B1 = g1[:, None, None, None]
    B2 = g2[None, :, None, None]
    C3 = g3[None, None, :, None]
    S = gs[None, None, None, :]
    sw = 0
    swr = 0
    swr2 = 0
    slv = 0
    for i in range(len(y)):
        V = S * S + sy[i] ** 2 + B1 * B1 * sl[i] ** 2 + B2 * B2 * sph[i] ** 2 + C3 * C3 * st[i] ** 2
        R = y[i] - B1 * xl[i] - B2 * xph[i] - C3 * t[i]
        w = 1 / V
        sw = sw + w
        swr = swr + w * R
        swr2 = swr2 + w * R * R
        slv = slv + np.log(V)
    lp = -0.5 * (slv + swr2 - swr * swr / sw) - 0.5 * np.log(sw)
    lp -= lp.max()
    p4 = np.exp(lp)
    p4 /= p4.sum()
    q1 = quant(g1, p4.sum((1, 2, 3)))
    q2 = quant(g2, p4.sum((0, 2, 3)))
    q3 = quant(g3, p4.sum((0, 1, 3)))
    out[tag] = dict(
        loo_bph_median_range=[min(meds), max(meds)],
        loo_bph_hi_max=max(his),
        loo_bph_lo_min=min(los),
        c_partial=qc,
        three=dict(b_line=q1, b_ph=q2, c=q3),
    )
    # if the widths followed L_ph^0.25 T with 0.05 dex of scatter: the partial slope recovered from the same predictors
    rng = np.random.default_rng(20261112)
    rec = []
    for it in range(40):
        ysim = 3.3 + 0.25 * xph + t + rng.normal(0, 0.05, len(y)) + rng.normal(0, 1, len(y)) * sy
        b1, b2, s, p = post3g(xl, xph, ysim, sl, sph, sy, n1=81, n2=81, ns=121)
        rec.append(quant(b2, p.sum((0, 2)))[2])
    out[tag]["power_bph_if_photospheric"] = [
        float(np.percentile(rec, 16)),
        float(np.median(rec)),
        float(np.percentile(rec, 84)),
    ]
json.dump(out, open(RESULTS / "joint_robustness.json", "w"), indent=1)
