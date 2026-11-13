"""Write the LaTeX tables of the manuscript: the sample (Table 1, stub and full CSV) and the fits (Table 2).

Usage: python scripts/make_tables.py
"""

import csv
import json

from lrdwidth.paths import DATA, PAPER, RESULTS

T = PAPER / "tables"
J = lambda n: json.load(open(RESULTS / n))
r1 = J("width_luminosity.json")
r2 = J("robustness.json")
src = DATA / "sample" / "table1_sources.csv"
rows = list(csv.DictReader(open(src)))
assert len(rows) == 44
open(T / "Table1_full.csv", "w").write(open(src).read())
key = ["RUBIES-EGS-42046", "CEERS-746", "RUBIES-UDS-119957", "Skyfire-299", "RUBIES-EGS-15825"]
t = [
    "\\begin{deluxetable*}{llcccllcccc}",
    "\\tabletypesize{\\footnotesize}",
    "\\tablewidth{0pt}",
    "\\tablecaption{The sample\\label{tab:sample}}",
    "\\tablehead{\\colhead{Source} & \\colhead{Class} & \\colhead{R.A.} & \\colhead{Decl.} & \\colhead{$z$} & \\colhead{Program} & \\colhead{Grating} & \\colhead{Label} & \\colhead{FWHM} & \\colhead{$\\log L_{\\rm H\\alpha}$} & \\colhead{S/N} \\\\",
    "\\colhead{} & \\colhead{} & \\colhead{(deg)} & \\colhead{(deg)} & \\colhead{} & \\colhead{} & \\colhead{} & \\colhead{} & \\colhead{(km s$^{-1}$)} & \\colhead{(erg s$^{-1}$)} & \\colhead{}}",
    "\\startdata",
]
for k in key:
    r = [x for x in rows if x["name"] == k][0]
    t.append(
        f"{r['name']} & {('LRD' if r['cls'] == 'LRD' else 'other')} & {float(r['ra']):.5f} & ${float(r['dec']):.5f}$ & {float(r['z']):.3f} & {r['program']} & {r['grating']} & {r['label_ref']} & ${r['fwhm']}\\pm{r['fwhm_err']}$ & ${float(r['logL']):.2f}\\pm{float(r['logL_err']):.2f}$ & {r['snr']} \\\\"
    )
t += [
    "\\enddata",
    "\\tablecomments{FWHM is the intrinsic width of the single broad Gaussian and $L_{\\rm H\\alpha}$ the luminosity of that component before absorption, without dust correction. Uncertainties are $1\\sigma$, conditional on the adopted model. Labels: T25 = \\citet{Taylor2025}; K26 = \\citet{Kocevski2026}; H25 = \\citet{Hviding2025}. Five rows are shown. The full table of 44 sources, with the five other width measures and all covariances, is published in machine-readable form.}",
    "\\end{deluxetable*}",
]
open(T / "Table1_stub.tex", "w").write("\n".join(t) + "\n")


def P(q):
    return (
        f"${q[2]:.3f}$ (${q[0]:.3f}$ to ${q[4]:.3f}$)"
        if q[0] < 0
        else f"${q[2]:.3f}$ (${q[0]:.3f}$--${q[4]:.3f}$)"
    )


ND = "\\nodata"


def line(lab, NL, bq, sL, NO, sO, dq):
    bs = (
        (
            f"${bq[2]:.2f}$ (${bq[0]:.2f}$ to ${bq[4]:.2f}$)"
            if bq[0] < 0
            else f"${bq[2]:.2f}$ (${bq[0]:.2f}$--${bq[4]:.2f}$)"
        )
        if bq
        else ND
    )
    return f"{lab} & {NL} & {bs} & {(P(sL) if sL else ND)} & {NO} & {(P(sO) if sO else ND)} & {(P(dq) if dq else ND)} \\\\"


W = lambda k: r1["width_" + k]
lh = r2["frozen_LSFhalf"]
n3 = r2["frozen_noTop3"]
nw = r2["frozen_noWide2"]
bm = r2["bx_main"]
bn = r2["bx_noOverlap"]
w6 = J("sdss_catalog_benchmark.json")
t = [
    "\\centerwidetable",
    "\\begin{deluxetable*}{lclclcl}",
    "\\tabletypesize{\\scriptsize}",
    "\\tabcolsep=1.7pt",
    "\\tablewidth{0pt}",
    "\\tablecaption{Width--luminosity fits\\label{tab:fits}}",
    "\\tablehead{\\colhead{} & \\multicolumn{3}{c}{LRDs} & \\multicolumn{2}{c}{Other broad-line sources} & \\colhead{} \\\\ \\cline{2-4} \\cline{5-6}",
    "\\colhead{Width measure or variant} & \\colhead{$N$} & \\colhead{$b$} & \\colhead{$\\sigma_{\\rm int}$ (dex)} & \\colhead{$N$} & \\colhead{$\\sigma_{\\rm int}$ (dex)} & \\colhead{$\\Delta\\sigma_{\\rm int}$ (dex)}}",
    "\\startdata",
    "\\cutinhead{This work: six measures of width}",
]
for lab, k in (
    ("Single-Gaussian FWHM", "single_gaussian"),
    ("Interquartile width $W_{50}$", "nonparametric_W50"),
    ("Exponential FWHM", "exponential"),
    ("Lorentzian FWHM", "lorentzian"),
    ("Nonparametric FWHM", "nonparametric_FWHM"),
    ("Two-Gaussian FWHM", "double_gaussian"),
):
    t.append(
        line(
            lab,
            29,
            W(k)["LRD"]["slope"],
            W(k)["LRD"]["sigma"],
            15,
            W(k)["other"]["sigma"],
            W(k)["diff"],
        )
    )
t.append("\\cutinhead{This work: single-Gaussian FWHM, variants}")
t.append(
    line(
        "LSF width $\\times0.5$",
        29,
        lh["LRD"]["slope"],
        lh["LRD"]["sigma"],
        15,
        lh["other"]["sigma"],
        lh["diff"],
    )
)
t.append(
    line(
        "Three brightest LRDs removed",
        26,
        n3["LRD"]["slope"],
        n3["LRD"]["sigma"],
        15,
        n3["other"]["sigma"],
        n3["diff"],
    )
)
t.append(
    line(
        "Two widest other sources removed",
        29,
        nw["LRD"]["slope"],
        nw["LRD"]["sigma"],
        13,
        nw["other"]["sigma"],
        nw["diff"],
    )
)
t.append("\\cutinhead{Table of \\citet{BaccusXu2026}}")
t.append(
    line(
        "All records",
        48,
        bm["LRD"]["slope"],
        bm["LRD"]["sigma"],
        97,
        bm["other"]["sigma"],
        bm["diff"],
    )
)
t.append(
    line(
        "Shared records removed",
        bn["LRD"]["N"],
        bn["LRD"]["slope"],
        bn["LRD"]["sigma"],
        bn["other"]["N"],
        bn["other"]["sigma"],
        bn["diff"],
    )
)
t.append(
    "Matched pairs\\tablenotemark{a} & 36 & \\nodata & \\nodata & 36 & \\nodata & $-0.065$ ($-0.105$ to $-0.024$) \\\\"
)
t.append("\\cutinhead{SDSS DR7 quasars \\citep{Shen2011}}")
t.append(
    f"Multi-Gaussian FWHM (catalog) & {ND} & {ND} & {ND} & {w6['N']} & {P(w6['sigma'])} & {ND} \\\\"
)
t.append(
    f"Single-Gaussian FWHM, refitted\\tablenotemark{{b}} & {ND} & {ND} & {ND} & 224 & $0.123$ ($0.112$--$0.136$) & {ND} \\\\"
)
t += [
    "\\enddata",
    "\\tablecomments{Posterior medians, with 95\\% intervals in parentheses. $b$ is the slope of $\\log{\\rm FWHM}$ on $\\log L$, and $\\Delta\\sigma_{\\rm int}$ is the scatter of the LRDs minus that of the other sources. Each width measure is regressed on the luminosity and covariance of its own fit. The luminosity is $L_{\\rm H\\alpha}$ for our sample and for SDSS, and $L_{\\rm bol}$ for the table of \\citet{BaccusXu2026}.}",
    "\\tablenotetext{a}{Pairs matched in grating, redshift ($|\\Delta z|<0.5$), and luminosity ($|\\Delta\\log L|<0.2$). The entry is the difference of the residual standard deviations, with its bootstrap interval.}",
    "\\end{deluxetable*}",
]
t.insert(
    len(t) - 1,
    "\\tablenotetext{b}{A random 300 of the 3755 quasars, refitted with the procedure of Section~\\ref{sec:fit}; 224 pass the rules of Section~\\ref{sec:sample}.}",
)
open(T / "Table2_fits.tex", "w").write("\n".join(t) + "\n")
print("tables written to", T)
