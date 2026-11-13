"""Draw Figures 1-5 of the manuscript (PDF and PNG) into paper/figures/.

The regressions shown are recomputed here from the per-source inputs with lrdwidth.stats.

Usage: python scripts/make_figures.py
"""

import csv
import json
import math

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lrdwidth import stats as S
from lrdwidth.paths import DATA, PAPER, PHOTOSPHERE, RESULTS, WIDTHS

mpl.use("Agg")
OUT = PAPER / "figures"
RED = "#B8323C"
BLUE = "#2F6FA7"
INK = "#1a1a1a"
MUTE = "#6b6b6b"
GOLD = "#B07A12"
mpl.rcParams.update(
    {
        "font.family": "STIXGeneral",
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.3,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.frameon": False,
        "axes.edgecolor": INK,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "pdf.fonttype": 42,
        "savefig.dpi": 300,
        "legend.handletextpad": 0.3,
    }
)
W1 = 3.4
W2 = 7.1


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf")
    fig.savefig(f"{OUT}/{name}.png")
    plt.close(fig)


def kms_axis(ax, lo, hi, ticks=(1000, 1500, 2000, 3000, 4000, 6000), axis="y"):
    tk = [t for t in ticks if lo <= t <= hi]
    mn = [np.log10(t) for t in range(500, 9000, 100) if lo <= t <= hi]
    if axis == "y":
        ax.set_ylim(np.log10(lo), np.log10(hi))
        ax.set_yticks([np.log10(t) for t in tk])
        ax.set_yticklabels([str(t) for t in tk])
        ax.set_yticks(mn, minor=True)
    else:
        ax.set_xlim(np.log10(lo), np.log10(hi))
        ax.set_xticks([np.log10(t) for t in tk])
        ax.set_xticklabels([str(t) for t in tk])
        ax.set_xticks(mn, minor=True)


def fit2(rows):
    x = np.array([r["x"] for r in rows])
    y = np.array([r["y"] for r in rows])
    C = np.array([r["cov"] for r in rows])
    b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
    return (float((p * ah).sum()), S.quant(b, p.sum(1))[2], S.quant(s, p.sum(0)))


# Figure 1: the relation in both data sets
W1in = json.load(open(WIDTHS / "regression_inputs.json"))
G = W1in["single_gaussian"]
Q = json.load(open(DATA / "bx26" / "regression_inputs.json"))
qso = json.load(open(DATA / "qso1" / "qso1_fit.json"))
ov = {
    r["published_id"]
    for r in csv.DictReader(open(DATA / "bx26" / "overlap_with_this_work.tsv"), delimiter="\t")
}
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 3.2))
for cls, col, mk, z in (("non-LRD", BLUE, "^", 2), ("LRD", RED, "o", 4)):
    r = [t for t in G if (t["cls"] == "LRD") == (cls == "LRD")]
    x = np.array([t["x"] for t in r]) + 43
    y = np.array([t["y"] for t in r])
    C = np.array([t["cov"] for t in r])
    a, b, sq = fit2(r)
    s = sq[2]
    xx = np.linspace(x.min() - 0.06, x.max() + 0.06, 50)
    yy = a + b * (xx - 43)
    a1.fill_between(
        xx, yy - s, yy + s, color=col, alpha=0.14 if cls == "LRD" else 0.09, lw=0, zorder=z - 2
    )
    a1.plot(xx, yy, color=col, lw=1.3, ls="-" if cls == "LRD" else (0, (4, 2)), zorder=z)
    a1.errorbar(
        x,
        y,
        xerr=np.sqrt(C[:, 0, 0]),
        yerr=np.sqrt(C[:, 1, 1]),
        fmt="none",
        ecolor=col,
        elinewidth=0.6,
        alpha=0.55,
        zorder=z,
    )
    a1.scatter(
        x,
        y,
        s=24 if mk == "o" else 28,
        marker=mk,
        facecolor=col,
        edgecolor="white",
        linewidth=0.6,
        zorder=z + 1,
    )
xq = math.log10(qso["L_delensed_erg_s"])
yq = math.log10(qso["FWHM_kms"])
eq = qso["FWHM_err_conditional_kms"] / qso["FWHM_kms"] / math.log(10)
exq = math.hypot(math.sqrt(qso["log_cov_xy"][0][0]), qso["mu_error"] / qso["mu"] / math.log(10))
a1.errorbar(xq, yq, xerr=exq, yerr=eq, fmt="none", ecolor=INK, elinewidth=0.7, zorder=6)
a1.scatter(xq, yq, s=34, marker="D", facecolor="white", edgecolor=INK, linewidth=0.9, zorder=7)
a1.annotate(
    "A2744-QSO1", (xq, yq), xytext=(6, 7), textcoords="offset points", fontsize=8, color=INK
)
a1.legend(
    handles=[
        Line2D([], [], marker="o", color=RED, mec="white", ms=5.5, lw=1.3, label="LRDs (29)"),
        Line2D(
            [],
            [],
            marker="^",
            color=BLUE,
            mec="white",
            ms=6,
            lw=1.3,
            ls=(0, (4, 2)),
            label="other broad-line sources (15)",
        ),
    ],
    loc="lower right",
    handlelength=2.6,
    borderaxespad=0.6,
)
a1.set_xlim(41.45, 43.72)
kms_axis(a1, 700, 4600)
a1.set_xlabel("$\\log\\,(L_{\\rm H\\alpha}\\,/\\,{\\rm erg\\,s^{-1}})$")
a1.set_ylabel("broad H$\\alpha$ FWHM  (km s$^{-1}$)")
a1.text(0.04, 0.95, "(a) this work", transform=a1.transAxes, va="top", fontsize=9.5)


def bxrow(t):
    return dict(
        id=t["id"],
        x=t["logLbol"] - 45,
        y=t["y"],
        cov=[[t["sx"] ** 2, 0], [0, t["sy"] ** 2]],
        z=t["z"],
        cls=t["cls"],
    )


R = [bxrow(t) for t in Q]
XL = [r for r in R if r["cls"] == "LRD"]
lo = min((r["x"] for r in XL))
hi = max((r["x"] for r in XL))
XC = [r for r in R if r["cls"] != "LRD" and r["z"] > 3 and (lo <= r["x"] <= hi)]
assert len(XL) == 48 and len(XC) == 97
rng = np.random.default_rng(20261112)
for rows, col, mk, z, key in ((XC, BLUE, "^", 2, "C"), (XL, RED, "o", 4, "LRD")):
    x = np.array([t["x"] for t in rows]) + 45
    y = np.array([t["y"] for t in rows])
    o = np.array([t["id"] in ov for t in rows])
    xj = x + rng.uniform(-0.032, 0.032, len(x))
    a, b, sq = fit2(rows)
    s = sq[2]
    xx = np.linspace(lo + 45 - 0.07, hi + 45 + 0.07, 50)
    yy = a + b * (xx - 45)
    a2.fill_between(
        xx, yy - s, yy + s, color=col, alpha=0.14 if key == "LRD" else 0.08, lw=0, zorder=z - 2
    )
    a2.plot(xx, yy, color=col, lw=1.3, ls="-" if key == "LRD" else (0, (4, 2)), zorder=z)
    sz = 22 if mk == "o" else 20
    a2.scatter(
        xj[~o],
        y[~o],
        s=sz,
        marker=mk,
        facecolor=col,
        edgecolor="white",
        linewidth=0.5,
        zorder=z + 1,
        alpha=1 if key == "LRD" else 0.8,
    )
    a2.scatter(
        xj[o],
        y[o],
        s=sz - 4,
        marker=mk,
        facecolor="white",
        edgecolor=col,
        linewidth=0.9,
        zorder=z + 1,
    )
a2.legend(
    handles=[
        Line2D([], [], marker="o", color=RED, mec="white", ms=5.5, lw=1.3, label="LRDs (48)"),
        Line2D(
            [],
            [],
            marker="^",
            color=BLUE,
            mec="white",
            ms=5.5,
            lw=1.3,
            ls=(0, (4, 2)),
            label="comparison sources (97)",
        ),
        Line2D(
            [],
            [],
            marker="o",
            color="none",
            mec=MUTE,
            mfc="white",
            ms=5,
            label="also in our sample",
        ),
    ],
    loc="lower right",
    handlelength=2.6,
    borderaxespad=0.6,
)
a2.set_xlim(43.55, 45.35)
kms_axis(a2, 540, 7600)
a2.set_xlabel("$\\log\\,(L_{\\rm bol}\\,/\\,{\\rm erg\\,s^{-1}})$, from broad H$\\alpha$")
a2.set_ylabel("broad H$\\alpha$ FWHM  (km s$^{-1}$)")
a2.text(0.04, 0.95, "(b) Baccus & Xu table", transform=a2.transAxes, va="top", fontsize=9.5)
fig.subplots_adjust(left=0.075, right=0.99, bottom=0.135, top=0.975, wspace=0.22)
save(fig, "Figure1")
# Figure 2: scatter under every width definition
r1 = json.load(open(RESULTS / "width_luminosity.json"))
r2 = json.load(open(RESULTS / "robustness.json"))
rows = [
    ("single-Gaussian FWHM", "width_single_gaussian"),
    ("interquartile width $W_{50}$", "width_nonparametric_W50"),
    ("exponential FWHM", "width_exponential"),
    ("Lorentzian FWHM", "width_lorentzian"),
    ("nonparametric FWHM", "width_nonparametric_FWHM"),
    ("two-Gaussian FWHM", "width_double_gaussian"),
]
fig, ax = plt.subplots(figsize=(W1, 3.75))
y = 0
yt = []
yl = []


def head(t, yy):
    ax.text(
        -0.03,
        yy,
        t,
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=8.6,
        fontweight="bold",
    )


def row(lab, qL, qO):
    global y
    for q, col, mk, dy in ((qL, RED, "o", -0.17), (qO, BLUE, "^", 0.17)):
        if q is None:
            continue
        ax.plot([q[0], q[4]], [y + dy] * 2, color=col, lw=0.9, solid_capstyle="butt")
        ax.plot([q[1], q[3]], [y + dy] * 2, color=col, lw=2.6, solid_capstyle="butt")
        ax.plot(
            q[2],
            y + dy,
            marker=mk,
            color=col,
            mec="white",
            mew=0.6,
            ms=5.2 if mk == "o" else 5.6,
            zorder=5,
        )
    yt.append(y)
    yl.append(lab)
    y += 1


head("This work", y)
y += 1
for lab, k in rows:
    row(lab, r1[k]["LRD"]["sigma"], r1[k]["other"]["sigma"])
row(
    "LSF width $\\times0.5$",
    r2["frozen_LSFhalf"]["LRD"]["sigma"],
    r2["frozen_LSFhalf"]["other"]["sigma"],
)
y += 0.35
head("Baccus & Xu", y)
y += 1
row("single-Gaussian FWHM", r2["bx_main"]["LRD"]["sigma"], r2["bx_main"]["other"]["sigma"])
sd = [0.112282, None, 0.12324, None, 0.136016]
ax.axvspan(sd[0], sd[4], color=INK, alpha=0.16, lw=0, zorder=0)
ax.axvline(sd[2], color=INK, lw=0.8, zorder=1)
ax.text(sd[2] + 0.008, 8.35, "SDSS quasars", va="center", ha="left", fontsize=8.3, color=INK)
ax.axvspan(0.05, 0.08, color=RED, alpha=0.07, lw=0, zorder=0)
ax.set_yticks(yt)
ax.set_yticklabels(yl)
ax.set_ylim(y - 0.35, -0.7)
ax.set_xlim(0, 0.47)
ax.set_xlabel("intrinsic scatter $\\sigma_{\\rm int}$ (dex)")
ax.tick_params(axis="y", which="both", left=False, right=False)
ax.legend(
    handles=[
        Line2D([], [], marker="o", color=RED, mec="white", ms=5.2, lw=0.9, label="LRDs"),
        Line2D(
            [],
            [],
            marker="^",
            color=BLUE,
            mec="white",
            ms=5.6,
            lw=0.9,
            label="other broad-line sources",
        ),
    ],
    loc="lower right",
    bbox_to_anchor=(1.0, 1.0),
    ncol=2,
    handlelength=1.8,
    borderaxespad=0.2,
    columnspacing=1.0,
)
fig.subplots_adjust(left=0.43, right=0.965, top=0.935, bottom=0.115)
save(fig, "Figure2")
# Figure 3: what the scatter allows if the widths are orbital
S0 = r1["width_single_gaussian"]["LRD"]["sigma"][2]
S1 = r1["width_single_gaussian"]["LRD"]["sigma"][4]
ori = r2["orient"]
fig, a = plt.subplots(figsize=(W1, 3.1))
th = np.linspace(0, np.pi / 2, 300)


def ell(s, sr=0.0):
    e = math.sqrt(s * s - sr * sr / 4)
    return (2 * e * np.cos(th), e * np.sin(th))


x1, y1 = ell(S1)
x0, y0 = ell(S0)
a.fill_between(x1[::-1], 0, y1[::-1], color=RED, alpha=0.07, lw=0)
a.fill_between(x0[::-1], 0, y0[::-1], color=RED, alpha=0.15, lw=0)
a.plot(x0, y0, color=RED, lw=1.5)
a.plot(x1, y1, color=RED, lw=1.1, ls=(0, (4, 2)))
lab = [
    ("thin disk, $i<45^\\circ$", ori["thin_0_45"], "bottom", "-"),
    ("$H/R=0.5$", ori["h0.5_0_60"], "bottom", (0, (5, 1.5, 1, 1.5))),
    ("$H/R=0.7$", ori["h0.7_0_60"], "bottom", (0, (5, 1.5, 1, 1.5))),
    ("$H/R=1$", ori["h1_0_60"], "top", (0, (5, 1.5, 1, 1.5))),
]
for k, v, va, ls in lab:
    a.axhline(v, color=BLUE, lw=0.8, ls=ls)
    a.text(
        0.236,
        v + (0.003 if va == "bottom" else -0.003),
        k,
        color=BLUE,
        fontsize=8.3,
        ha="right",
        va=va,
    )
a.text(0.01, 0.012, "allowed", color=RED, fontsize=9)
a.legend(
    handles=[
        Line2D([], [], color=RED, lw=1.5, label="$\\sigma_{\\rm int}=%.3f$ dex" % S0),
        Line2D(
            [], [], color=RED, lw=1.1, ls=(0, (4, 2)), label="$%.3f$ dex (95%% upper limit)" % S1
        ),
    ],
    loc="upper left",
    bbox_to_anchor=(0, 0.84),
    handlelength=2.4,
    borderaxespad=0.7,
)
a.set_xlim(0, 0.24)
a.set_ylim(0, 0.235)
a.set_xlabel("scatter of $\\log(M/R)$ at fixed luminosity (dex)")
a.set_ylabel("orientation scatter of $\\log$ FWHM (dex)")
fig.subplots_adjust(left=0.155, right=0.975, bottom=0.14, top=0.975)
save(fig, "Figure3")
# Figures 4 and 5: the photosphere
# v4 restructuring (2026-11-13): Figure 4 keeps only the comparison with v_K (single column). The posterior of the exponents
# moves to the appendix as Figure 5. The added-variable panel of v3 is gone (the regression at fixed line luminosity is
# reported as numbers only).
SXMAX = 0.3


def load(fn):
    rr = list(csv.DictReader(open(PHOTOSPHERE / fn), delimiter="\t"))
    f = lambda k: np.array([float(r[k]) for r in rr])
    return (rr, f("x"), f("t"), f("y"), f("sx"), f("st"), f("sy"))


fig, a = plt.subplots(figsize=(W1, 3.25))
rr, x, t, y, sx, st, sy = load("atmosphere_bx26_widths.tsv")
k = sx < SXMAX
v = np.log10(S.VK0_KMS) + 0.25 * x + t
ev = np.sqrt((0.25 * sx) ** 2 + st**2)
a.errorbar(
    v[k], y[k], xerr=ev[k], yerr=sy[k], fmt="none", ecolor=RED, elinewidth=0.6, alpha=0.6, zorder=3
)
a.scatter(
    v[k],
    y[k],
    s=22,
    marker="o",
    facecolor=RED,
    edgecolor="white",
    linewidth=0.5,
    zorder=4,
    label="atmosphere fits (%d)" % k.sum(),
)
rr, x, t, y, sx, st, sy = load("blackbody_sources.tsv")
own = np.array([r["width_source"] != "BX26_dedup" for r in rr])
v = np.log10(S.VK0_KMS) + 0.25 * x + t
ev = np.sqrt((0.25 * sx) ** 2 + st**2)
a.errorbar(v, y, xerr=ev, yerr=sy, fmt="none", ecolor=GOLD, elinewidth=0.6, alpha=0.6, zorder=3)
a.scatter(
    v[own],
    y[own],
    s=24,
    marker="s",
    facecolor=GOLD,
    edgecolor="white",
    linewidth=0.5,
    zorder=5,
    label="blackbody fits, our widths (%d)" % own.sum(),
)
a.scatter(
    v[~own],
    y[~own],
    s=20,
    marker="s",
    facecolor="white",
    edgecolor=GOLD,
    linewidth=1.0,
    zorder=5,
    label="blackbody fits, B&X widths (%d)" % (~own).sum(),
)
xx = np.array([2.55, 3.25])
for m, ls in ((1, (0, (1, 1.6))), (2, (0, (4, 2))), (3, "-")):
    a.plot(xx, xx + math.log10(m), color=INK, lw=0.9, ls=ls, zorder=1)
for m in (1, 2, 3):
    a.text(
        3.13,
        3.13 + math.log10(m) - 0.014,
        "$%d\\,v_{\\rm K}$" % m if m > 1 else "$v_{\\rm K}$",
        fontsize=8,
        rotation=37,
        ha="center",
        va="top",
    )
kms_axis(a, 430, 1500, (500, 700, 1000, 1400), axis="x")
kms_axis(a, 760, 4600, (1000, 1500, 2000, 3000, 4000))
a.set_xlabel("$v_{\\rm K}(R_{\\rm ph})$ for $\\lambda=1$  (km s$^{-1}$)")
a.set_ylabel("broad H$\\alpha$ FWHM  (km s$^{-1}$)")
a.legend(loc="upper left", handlelength=1.0, borderaxespad=0.5, fontsize=7.6, handletextpad=0.4)
fig.subplots_adjust(left=0.165, right=0.975, bottom=0.135, top=0.975)
save(fig, "Figure4")
fig, b = plt.subplots(figsize=(W1, 3.1))
LAB = dict(A=(0.265, -0.05), B=(0.09, 0.95), C=(0.44, -0.12))
for tag, fn, col in (
    ("A", "atmosphere_bx26_widths.tsv", RED),
    ("B", "atmosphere_our_widths.tsv", BLUE),
    ("C", "blackbody_sources.tsv", GOLD),
):
    rr, x, t, y, sx, st, sy = load(fn)
    k = sx < SXMAX
    g1, g2, s, p = S.post3(x[k], t[k], y[k], sx[k], st[k], sy[k])
    pb = p.sum(2)
    l68, l95 = S.contour_levels(pb)
    b.contour(g1, g2, pb.T, levels=[l95], colors=[col], linewidths=0.9, linestyles=[(0, (4, 2))])
    b.contour(g1, g2, pb.T, levels=[l68], colors=[col], linewidths=1.3)
    b.text(
        *LAB[tag], tag, color=col, fontsize=9, fontweight="bold", ha="center", va="center", zorder=7
    )
b.plot(0.25, 1, marker="*", color=INK, ms=10, mec="white", mew=0.6, zorder=6)
b.annotate(
    "photosphere",
    (0.25, 1),
    xytext=(0.36, 1.6),
    fontsize=8,
    arrowprops=dict(arrowstyle="-", lw=0.5, color=INK),
    ha="left",
)
FS = json.load(open(PHOTOSPHERE / "simulation_results.json"))["results"]
pts = []
for key in ("Curtis_A_own", "Curtis_B_own"):
    h = FS[key]["scenarios"]["measured"]["methods"]["flat"]["Hline_bc"]
    pts.append((h["b_distribution"]["values"][2], h["c_distribution"]["values"][2]))
    b.plot(*pts[-1], marker="P", color=INK, ms=6.5, mec="white", mew=0.6, zorder=6)
    b.annotate(
        key[7],
        pts[-1],
        xytext=(5.5, 1),
        textcoords="offset points",
        fontsize=6.5,
        color=INK,
        va="center",
        zorder=7,
    )
for pt in pts:
    b.annotate("", pt, xytext=(-0.02, -1.25), arrowprops=dict(arrowstyle="-", lw=0.5, color=INK))
b.text(-0.03, -1.3, "$L_{\\rm H\\alpha}$ only", fontsize=8, ha="left", va="top")
b.set_xlim(-0.18, 0.72)
b.set_ylim(-2.9, 2.3)
b.set_xlabel("exponent of $L_{\\rm ph}$")
b.set_ylabel("exponent of $T_{\\rm eff}$")
b.legend(
    handles=[
        Line2D([], [], color=RED, lw=1.3, label="A: atmosphere fits, B&X widths"),
        Line2D([], [], color=BLUE, lw=1.3, label="B: atmosphere fits, our widths"),
        Line2D([], [], color=GOLD, lw=1.3, label="C: blackbody fits (9)"),
    ],
    loc="lower right",
    handlelength=1.4,
    borderaxespad=0.4,
    fontsize=7.6,
)
fig.subplots_adjust(left=0.15, right=0.975, bottom=0.14, top=0.975)
save(fig, "Figure5")
print("figures written to", OUT)
