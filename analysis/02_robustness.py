"""Robustness of the relation, the BX26 table, virial-mass residuals, same-source comparisons with photospheric
parameters, and the orientation scatter of disk-like velocity fields.

Writes data/results/robustness.json.
"""

import csv
import json

import numpy as np

from lrdwidth.paths import DATA, PHOTOSPHERE, RESULTS, WIDTHS
from lrdwidth.regression import ols_sd, ori, pair, post2, quant, sig1d

out = {}
# main sample and robustness
W1 = json.load(open(WIDTHS / "regression_inputs.json"))["single_gaussian"]
L = [r for r in W1 if r["cls"] == "LRD"]
O = [r for r in W1 if r["cls"] != "LRD"]
out["frozen_main"] = pair(L, O)
top = sorted(L, key=lambda r: -r["x"])[:3]
L3 = [r for r in L if r not in top]
out["frozen_noTop3"] = pair(L3, O)
out["frozen_noTop3"]["dropped"] = [
    (t["source_id"], round(t["x"] + 43, 2), bool(t.get("absorption"))) for t in top
]
wid = sorted(O, key=lambda r: -r["y"])[:2]
O2 = [r for r in O if r not in wid]
out["frozen_noWide2"] = pair(L, O2)
out["frozen_noWide2"]["dropped"] = [(t["source_id"], round(10 ** t["y"])) for t in wid]
noabs = [r for r in L if not r.get("absorption")]
out["frozen_noAbs"] = pair(noabs, O)
out["N_abs"] = [
    sum((bool(r.get("absorption")) for r in L)),
    sum((bool(r.get("absorption")) for r in O)),
]
# descriptive ranges
desc = lambda rr: dict(
    L=[
        float(min((r["x"] for r in rr)) + 43),
        float(np.median([r["x"] for r in rr]) + 43),
        float(max((r["x"] for r in rr)) + 43),
    ],
    W=[
        float(10 ** min((r["y"] for r in rr))),
        float(10 ** np.median([r["y"] for r in rr])),
        float(10 ** max((r["y"] for r in rr))),
    ],
)
out["frozen_desc"] = dict(LRD=desc(L), other=desc(O))
# instrumental width halved
W3 = json.load(open(WIDTHS / "regression_inputs_lsf_half.json"))
L5 = [r for r in W3 if r["cls"] == "LRD"]
O5 = [r for r in W3 if r["cls"] != "LRD"]
out["frozen_LSFhalf"] = pair(L5, O5)
dl = np.array([abs(r["delta_logwidth_dex"]) for r in W3])
out["LSFhalf_absdelta"] = dict(median=float(np.median(dl)), max=float(dl.max()))
# the BX26 table
Q = json.load(open(DATA / "bx26" / "regression_inputs.json"))


def bxrow(t):
    return dict(
        id=t["id"],
        x=t["logLbol"] - 45,
        y=t["y"],
        cov=[[t["sx"] ** 2, 0], [0, t["sy"] ** 2]],
        z=t["z"],
        cls=t["cls"],
        logM=t["logM"],
        g=t["grating"],
    )


R = [bxrow(t) for t in Q]
XL = [r for r in R if r["cls"] == "LRD"]
lo = min((r["x"] for r in XL))
hi = max((r["x"] for r in XL))
XC = [r for r in R if r["cls"] != "LRD" and r["z"] > 3 and (lo <= r["x"] <= hi)]
out["bx_main"] = pair(XL, XC)
out["bx_desc"] = dict(
    Lrange=[lo + 45, hi + 45],
    zmin_LRD=min((r["z"] for r in XL)),
    NLRD_zlt3=sum((r["z"] < 3 for r in XL)),
)
ov = {
    r["published_id"]
    for r in csv.DictReader(open(DATA / "bx26" / "overlap_with_this_work.tsv"), delimiter="\t")
}
XLn = [r for r in XL if r["id"] not in ov]
lo2 = min((r["x"] for r in XLn))
hi2 = max((r["x"] for r in XLn))
XCn = [
    r
    for r in R
    if r["cls"] != "LRD" and r["z"] > 3 and (lo <= r["x"] <= hi) and (r["id"] not in ov)
]
out["bx_noOverlap"] = pair(XLn, XCn)
# virial masses: residuals of the tabulated log M about a line in log L_bol
for nm, rr in (("LRD", XL), ("ctrl", XC)):
    x = np.array([r["x"] for r in rr])
    m = np.array([r["logM"] for r in rr])
    sd, sl = ols_sd(x, m)
    lam = x + 45 - m
    out["bx_Mvir_" + nm] = dict(resid_sd=sd, slope=sl, sd_loglam=float(np.std(lam, ddof=1)))
# same sources: scatter of the widths about the line luminosity and about photospheric quantities
byid_f = {r["source_id"]: r for r in W1}
byid_x = {r["id"]: r for r in R}
for tag, fn in (
    ("A", "atmosphere_bx26_widths.tsv"),
    ("B", "atmosphere_our_widths.tsv"),
    ("C", "blackbody_sources.tsv"),
):
    rows = list(csv.DictReader(open(PHOTOSPHERE / fn), delimiter="\t"))
    x = np.array([float(r["x"]) for r in rows])
    t = np.array([float(r["t"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    sx = np.array([float(r["sx"]) for r in rows])
    st = np.array([float(r["st"]) for r in rows])
    sy = np.array([float(r["sy"]) for r in rows])
    res = dict(
        N=len(rows),
        sources=[r["width_id"] for r in rows],
        width_source=sorted(set((r["width_source"] for r in rows))),
    )
    z0 = np.zeros(len(x))
    # (i) photospheric luminosity only (c = 0)
    b, s, p = post2(x, y, sx**2, sy**2, z0, a_lo=-50, a_hi=50)
    res["vs_Lph"] = dict(sigma=quant(s, p.sum(0)), slope=quant(b, p.sum(1)))
    # (ii) exponents fixed at (0.25, 1)
    q = y - 0.25 * x - t
    res["fixed_pred"] = dict(sigma=sig1d(q, sy**2 + (0.25 * sx) ** 2 + st**2))
    # (iii) temperature only
    b, s, p = post2(t, y, st**2, sy**2, z0, b_lo=-3, b_hi=3, nb=901, a_lo=-50, a_hi=50)
    res["vs_T"] = dict(sigma=quant(s, p.sum(0)), slope=quant(b, p.sum(1)))
    # (iv) Halpha luminosity, same sources
    xa = []
    ca = []
    ok = []
    for r in rows:
        wid_ = r["width_id"]
        src_ = r["width_source"]
        if wid_ in byid_f and src_ in ("frozen",):
            xa.append(byid_f[wid_]["x"])
            ca.append(byid_f[wid_]["cov"])
            ok.append(True)
        elif wid_ in byid_x:
            xa.append(byid_x[wid_]["x"])
            ca.append(byid_x[wid_]["cov"])
            ok.append(True)
        else:
            xa.append(np.nan)
            ca.append([[0, 0], [0, 0]])
            ok.append(False)
    ok = np.array(ok)
    xa = np.array(xa)
    ca = np.array(ca)
    res["N_with_line_lum"] = int(ok.sum())
    if tag in "AB" and ok.all():
        b, s, p = post2(xa, y, ca[:, 0, 0], sy**2, ca[:, 0, 1], a_lo=-50, a_hi=50)
        res["vs_Lline"] = dict(sigma=quant(s, p.sum(0)), slope=quant(b, p.sum(1)))
        res["corr_Lline_Lph"] = float(np.corrcoef(xa, x)[0, 1])
        res["ols_Lph_on_Lline"] = ols_sd(xa, x)
    out["pred_" + tag] = res
out["orient"] = {
    "thin_0_45": ori(1e-06, 45, 0),
    "thin_0_60": ori(1e-06, 60, 0),
    "thin_25_45": ori(25, 45, 0),
    "h0.5_0_60": ori(0, 60, 0.5),
    "h0.7_0_60": ori(0, 60, 0.7),
    "h1_0_60": ori(0, 60, 1.0),
    "h0.3_0_60": ori(0, 60, 0.3),
    "h0.5_0_90": ori(0, 90, 0.5),
    "h1_0_90": ori(0, 90, 1.0),
    "thin_45_90": ori(45, 90, 0),
    "thin_60_90": ori(60, 90, 0),
}
json.dump(out, open(RESULTS / "robustness.json", "w"), indent=1, default=str)
