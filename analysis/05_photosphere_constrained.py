"""Photospheric comparison restricted to sources whose photospheric luminosity is constrained
(sigma of log L_ph below SXMAX).

The likelihood V = s^2 + sy^2 + b^2 sx^2 + c^2 st^2 assumes a flat prior on the true x. When sx far exceeds the spread
of x in the sample that prior fails and the -0.5 ln V term biases |b| low, so such sources are excluded.

Writes data/results/photosphere_constrained.json.
"""

import csv
import json
import math

import numpy as np

from lrdwidth.paths import DATA, PHOTOSPHERE, RESULTS
from lrdwidth.regression import (
    hpd_level,
    load_line_inputs,
    load_photosphere_sample,
    post2,
    post3g,
    post_lum_temp,
    quant,
    sig1d,
)

SXMAX = 0.3
W1, byf, Q, byx = load_line_inputs(DATA)
out = {"SXMAX": SXMAX}
for tag, fn in (
    ("A", "atmosphere_bx26_widths.tsv"),
    ("B", "atmosphere_our_widths.tsv"),
    ("C", "blackbody_sources.tsv"),
):
    rows = list(csv.DictReader(open(PHOTOSPHERE / fn), delimiter="\t"))
    keep = np.array([float(r["sx"]) < SXMAX for r in rows])
    x = np.array([float(r["x"]) for r in rows])[keep]
    t = np.array([float(r["t"]) for r in rows])[keep]
    y = np.array([float(r["y"]) for r in rows])[keep]
    sx = np.array([float(r["sx"]) for r in rows])[keep]
    st = np.array([float(r["st"]) for r in rows])[keep]
    sy = np.array([float(r["sy"]) for r in rows])[keep]
    q = y - (0.25 * x + t + math.log10(1492.630535))
    b, cc, s, p, ah = post_lum_temp(x, t, y, sx, st, sy, 161, 241, 301)
    pbc = p.sum(2)
    r = dict(
        N=int(keep.sum()),
        dropped=[rw["temperature_id"] for rw, k in zip(rows, keep) if not k],
        ratio_median=float(10 ** np.median(q)),
        ratio_p16=float(10 ** np.percentile(q, 16)),
        ratio_p84=float(10 ** np.percentile(q, 84)),
        q_sd=float(np.std(q, ddof=1)),
        b=quant(b, p.sum((1, 2))),
        c=quant(cc, p.sum((0, 2))),
        sigma=quant(s, p.sum((0, 1))),
        hpd=hpd_level(b, cc, pbc, 0.25, 1.0),
        hpd_b025_c0=hpd_level(b, cc, pbc, 0.25, 0.0),
        P_c_lt_1=float(p.sum((0, 2))[cc < 1].sum()),
        P_b_lt_025=float(p.sum((1, 2))[b < 0.25].sum()),
        sigma_fixed=sig1d(q, sy**2 + (0.25 * sx) ** 2 + st**2),
        sd_t=float(np.std(t, ddof=1)),
        sd_x=float(np.std(x, ddof=1)),
        logL_range=[float(x.min() + 45), float(x.max() + 45)],
        T_range=[float(5000 * 10 ** t.min()), float(5000 * 10 ** t.max())],
    )
    # normalization at fixed exponents (weighted mean of the ratio) and its error
    var = sy**2 + (0.25 * sx) ** 2 + st**2 + r["sigma_fixed"][2] ** 2
    w = 1 / var
    r["q_mean"] = float((w * q).sum() / w.sum())
    r["q_mean_err"] = float(1 / math.sqrt(w.sum()))
    r["ratio_mean"] = 10 ** r["q_mean"]
    # c fixed at 1 or at 0, fitting b only
    for cfix in (0, 1):
        bb, ss, pp = post2(
            x, y - cfix * t, sx**2, sy**2 + (cfix * st) ** 2, np.zeros(len(x)), a_lo=-50, a_hi=50
        )
        r[f"b_cfix{cfix}"] = quant(bb, pp.sum(1))
        r[f"sigma_cfix{cfix}"] = quant(ss, pp.sum(0))
    out[tag] = r
    if tag in "AB":
        rows2, xl, sl, xph, sph, t2, st2, y2, sy2 = load_photosphere_sample(
            PHOTOSPHERE / fn, tag, byf, byx
        )
        xl = xl[keep]
        sl = sl[keep]
        b1, b2, s, p = post3g(xl, x, y, sl, sx, sy)
        pb = p.sum(2)
        r["joint"] = dict(
            b_line=quant(b1, p.sum((1, 2))),
            b_ph=quant(b2, p.sum((0, 2))),
            sigma=quant(s, p.sum((0, 1))),
            hpd_0_025=hpd_level(b1, b2, pb, 0.0, 0.25),
        )
        bb, ss, pp = post2(xl, y, sl**2, sy**2, np.zeros(len(y)), a_lo=-50, a_hi=50)
        r["vs_line"] = dict(sigma=quant(ss, pp.sum(0)), slope=quant(bb, pp.sum(1)))
        # regress on the line luminosity and on the photospheric prediction x_pred = 0.25 x_ph + t: y = a + b1 x_line + g x_pred (photospheric scaling: g = 1, b1 = 0)
        xp = 0.25 * x + t
        sp = np.sqrt((0.25 * sx) ** 2 + st**2)
        b1, g, s, p = post3g(xl, xp, y, sl, sp, sy, lim2=(-1, 3), n2=241)
        pb = p.sum(2)
        r["joint_pred"] = dict(
            b_line=quant(b1, p.sum((1, 2))),
            g=quant(g, p.sum((0, 2))),
            sigma=quant(s, p.sum((0, 1))),
            hpd_photospheric=hpd_level(b1, g, pb, 0.0, 1.0),
            hpd_line_only=float(
                pb[
                    pb > pb[int(np.argmax(pb[:, int(np.argmin(abs(g)))])), int(np.argmin(abs(g)))]
                ].sum()
            ),
        )
        r["corr_line_ph"] = float(np.corrcoef(xl, x)[0, 1])
        c2 = np.polyfit(xl, x, 1)
        r["sd_Lph_at_fixed_line"] = float(np.std(x - np.polyval(c2, xl), ddof=2))
json.dump(out, open(RESULTS / "photosphere_constrained.json", "w"), indent=1)
