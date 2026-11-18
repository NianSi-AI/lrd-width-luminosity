"""Width-luminosity relation for six measures of width, and the widths against the photospheric Keplerian speed.

Writes data/results/width_luminosity.json.
"""

import csv
import json

import numpy as np

from lrdwidth import stats
from lrdwidth.paths import PHOTOSPHERE, RESULTS, WIDTHS
from lrdwidth.regression import diffq, hpd_level, ols_sd, post2, post_lum_temp, quant

kap, coef, vK0 = stats.KAPPA, stats.VK_COEF, stats.VK0_KMS
out = {}
# six measures of width
W1 = json.load(open(WIDTHS / "regression_inputs.json"))
for kind, rows in W1.items():
    res = {}
    marg = {}
    for cls, sel in (("LRD", lambda r: r["cls"] == "LRD"), ("other", lambda r: r["cls"] != "LRD")):
        rr = [r for r in rows if sel(r)]
        x = np.array([r["x"] for r in rr])
        y = np.array([r["y"] for r in rr])
        C = np.array([r["cov"] for r in rr])
        b, s, p = post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
        marg[cls] = p.sum(0)
        sd, sl = ols_sd(x, y)
        res[cls] = dict(
            N=len(rr), sigma=quant(s, p.sum(0)), slope=quant(b, p.sum(1)), ols_sd=sd, ols_slope=sl
        )
    res["diff"] = diffq(s, marg["LRD"], marg["other"])
    out["width_" + kind] = res
# Keplerian speed at the photosphere
out["photosphere_coef"] = dict(kappa=kap, coef_cgs=coef, vK_1e45_5000K_kms=vK0)
for tag, fn in (
    ("A", "atmosphere_bx26_widths.tsv"),
    ("B", "atmosphere_our_widths.tsv"),
    ("C", "blackbody_sources.tsv"),
):
    rows = list(csv.DictReader(open(PHOTOSPHERE / fn), delimiter="\t"))
    T = np.array([float(r["T_K"]) for r in rows])
    logL = np.array([float(r["logL_ph_erg_s"]) for r in rows])
    F = np.array([float(r["FWHM_kms"]) for r in rows])
    sx = np.array([float(r["sx"]) for r in rows])
    st = np.array([float(r["st"]) for r in rows])
    sy = np.array([float(r["sy"]) for r in rows])
    x = logL - 45
    t = np.log10(T / 5000)
    y = np.log10(F)
    vK = coef * 10 ** (logL / 4) * T / 100000.0
    q = y - np.log10(vK)
    assert np.allclose(x, [float(r["x"]) for r in rows], atol=1e-09) and np.allclose(
        t, [float(r["t"]) for r in rows], atol=1e-09
    )
    dq = float(np.max(np.abs(q - np.array([float(r["q_dex"]) for r in rows]))))
    r = {}
    for gname, (nb, nc, ns) in (("g1", (121, 181, 241)), ("g2", (161, 241, 301))):
        b, cc, s, p, ah = post_lum_temp(x, t, y, sx, st, sy, nb, nc, ns)
        pbc = p.sum(2)
        r[gname] = dict(
            b=quant(b, p.sum((1, 2))),
            c=quant(cc, p.sum((0, 2))),
            sigma=quant(s, p.sum((0, 1))),
            a=float((p * ah).sum()),
            hpd=hpd_level(b, cc, pbc, 0.25, 1.0),
            hpd_c0=hpd_level(b, cc, pbc, 0.25, 0.0),
            P_c_lt_1=float(p.sum((0, 2))[cc < 1].sum()),
            P_b_lt_025=float(p.sum((1, 2))[b < 0.25].sum()),
        )
    # exponents fixed at (0.25, 1): only the normalization and the scatter remain
    sq2 = sy**2 + (0.25 * sx) ** 2 + st**2
    s = (np.arange(2000) + 0.5) * 0.6 / 2000
    S = s[None, :]
    V = S * S + sq2[:, None]
    w = 1 / V
    sw = w.sum(0)
    swr = (w * q[:, None]).sum(0)
    swr2 = (w * q[:, None] ** 2).sum(0)
    lp = -0.5 * (np.log(V).sum(0) + swr2 - swr**2 / sw) - 0.5 * np.log(sw)
    lp -= lp.max()
    pp = np.exp(lp)
    pp /= pp.sum()
    qmean = float((pp * swr / sw).sum())
    r.update(
        N=len(rows),
        q_median=float(np.median(q)),
        q_p16=float(np.percentile(q, 16)),
        q_p84=float(np.percentile(q, 84)),
        q_sd=float(np.std(q, ddof=1)),
        ratio_median=float(10 ** np.median(q)),
        max_abs_dq_vs_table=dq,
        sigma_q_fixed=quant(s, pp),
        q_mean_fixed=qmean,
        ratio_mean_fixed=float(10**qmean),
        med_sx=float(np.median(sx)),
        med_st=float(np.median(st)),
        med_sy=float(np.median(sy)),
        T_range=[float(T.min()), float(T.max())],
        logL_range=[float(logL.min()), float(logL.max())],
        sd_t=float(np.std(t, ddof=1)),
        sd_x=float(np.std(x, ddof=1)),
        corr_xt=float(np.corrcoef(x, t)[0, 1]),
    )
    out["photosphere_" + tag] = r
json.dump(out, open(RESULTS / "width_luminosity.json", "w"), indent=1)
