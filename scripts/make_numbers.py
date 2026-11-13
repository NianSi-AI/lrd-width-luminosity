"""Compute every number quoted in the manuscript and write paper/numbers.tex.

Each number is a LaTeX macro. Only the macros that paper/ms.tex uses are written to numbers.tex; all of them are
written to data/results/numbers.json. Inputs are the stored results in data/ (see data/README.md).

Usage: python scripts/make_numbers.py
"""

import csv
import json
import math
import re

from lrdwidth.paths import DATA, PAPER, PHOTOSPHERE, RESULTS


def load(name):
    return json.load(open(RESULTS / name))


r1 = load("width_luminosity.json")
r2 = load("robustness.json")
r6 = load("photosphere_constrained.json")
r11 = load("photosphere_simulations.json")
sdss_cat = load("sdss_catalog_benchmark.json")
X = load("measurement_summary.json")["numbers"]
FS = json.load(open(PHOTOSPHERE / "simulation_results.json"))["results"]
xv = lambda k: X[k]["value"]
M = []  # (macro name, value); ("%", title) entries mark sections


def add(name, val):
    assert name.isalpha(), name
    assert name not in [m[0] for m in M], name
    M.append((name, str(val)))


f3 = lambda v: f"{v:.3f}"
f2 = lambda v: f"{v:.2f}"


def neg(s):
    return "\\ensuremath{-}" + s[1:] if s.startswith("-") else s


def trip(prefix, q, fmt):
    add(prefix, fmt(q[2]))
    add(prefix + "lo", fmt(q[0]))
    add(prefix + "hi", fmt(q[4]))


sec = lambda t: M.append(("%", t))
sec("sample (design and counts)")
for n, v in (
    ("NLRD", 29),
    ("NOther", 15),
    ("NFrozen", 44),
    ("NSources", 115),
    ("NSpectra", 204),
    ("NSpecHa", 177),
    ("NProg", "eight"),
    ("NEGS", 30),
    ("NUDS", 14),
    ("NAbsLRD", "five"),
    ("NAbsOther", "four"),
    ("NLabTLRD", 13),
    ("NLabKLRD", 10),
    ("NLabHLRD", "six"),
):
    add(n, v)
assert r2["N_abs"] == [5, 4]
for n, v in (("cutSN", 5), ("cutBIC", 6), ("cutChi", 2), ("HvidingCut", 1000)):
    add(n, v)
add("matchRadius", "0\\farcs5")
for n, v in (
    ("zLRDmin", "3.36"),
    ("zLRDmax", "7.15"),
    ("zLRDmed", "5.45"),
    ("SNLRDmed", "13.7"),
    ("SNOthmed", "8.2"),
    ("zOthmin", "3.60"),
    ("zOthmax", "6.04"),
):
    add(n, v)
d = r2["frozen_desc"]
add("LLRDmin", f2(d["LRD"]["L"][0]))
add("LLRDmax", f2(d["LRD"]["L"][2]))
add("WLRDmin", round(d["LRD"]["W"][0]))
add("WLRDmax", round(d["LRD"]["W"][2]))
add("WLRDmed", round(d["LRD"]["W"][1]))
sec("line profiles (BIC)")
b = "W1_BIC_counts__LRD_all__"
add("NBICwing", xv(b + "best_counts__lorentzian") + xv(b + "best_counts__exponential"))
add("NBICsix", xv(b + "Gaussian_deltaBIC_gt6"))
add("NBIChiSN", "nine")
add("NBICsixOth", "four")
add("NCoreWing", "seven")
sec("our sample, single-Gaussian FWHM")
m = r1["width_single_gaussian"]
fm = r2["frozen_main"]
trip("Fsig", m["LRD"]["sigma"], f3)
trip("Fb", m["LRD"]["slope"], f2)
add("Fa", f3(fm["LRD"]["a"]))
add("FWnorm", int(round(10 ** fm["LRD"]["a"], -1)))
trip("Fsigo", m["other"]["sigma"], f3)
trip("Fdsig", m["diff"], lambda v: neg(f3(v)))
add("Frms", f3(m["LRD"]["ols_sd"]))
add("Frmso", f3(m["other"]["ols_sd"]))
n3 = r2["frozen_noTop3"]
trip("FbNoTop", n3["LRD"]["slope"], f2)
trip("FsigNoTop", n3["LRD"]["sigma"], f3)
add("NBright", "three")
add("LBright", "43.5")
add("LRestMax", "42.9")
nw = r2["frozen_noWide2"]
trip("FsigoDrop", nw["other"]["sigma"], f3)
trip("FdsigDrop", nw["diff"], lambda v: neg(f3(v)))
add("WOutA", nw["dropped"][0][1])
add("WOutB", nw["dropped"][1][1])
lh = r2["frozen_LSFhalf"]
trip("FsigLsf", lh["LRD"]["sigma"], f3)
trip("FsigoLsf", lh["other"]["sigma"], f3)
trip("FdsigLsf", lh["diff"], lambda v: neg(f3(v)))
sec("six width definitions")
defs = (
    ("Wf", "nonparametric_W50"),
    ("Ex", "exponential"),
    ("Lo", "lorentzian"),
    ("Np", "nonparametric_FWHM"),
    ("Tg", "double_gaussian"),
)
for p, k in defs:
    q = r1["width_" + k]
    trip("Dsig" + p, q["LRD"]["sigma"], f3)
    trip("Dsigo" + p, q["other"]["sigma"], f3)
    trip("Ddsig" + p, q["diff"], lambda v: neg(f3(v)))
    add("Db" + p, f2(q["LRD"]["slope"][2]))
five = [
    r1["width_" + k]["LRD"]["sigma"][2]
    for k in (
        "single_gaussian",
        "nonparametric_W50",
        "exponential",
        "lorentzian",
        "nonparametric_FWHM",
    )
]
add("sigFiveLo", f2(min(five)))
add("sigFiveHi", f2(max(five)))
sec("injection-recovery")
g = lambda sc, inj, meas, what: xv(
    f"W2_final_summary_s{sc}__{sc}_{inj}_M4_Gaussian_QC_{meas}__sigma_median__{what}"
)
add("simInj", "0.05")
add("simNull", "0.15")
add("NReal", 100)
add("recGG", f3(g(3, "single_gaussian", "single_gaussian", "p16_p50_p84__p50")))
add("recEG", f3(g(3, "exponential", "single_gaussian", "p16_p50_p84__p50")))
add("recGW", f3(g(3, "single_gaussian", "nonparametric_W50", "p16_p50_p84__p50")))
add("recEW", f3(g(3, "exponential", "nonparametric_W50", "p16_p50_p84__p50")))
add("recGL", f2(g(3, "single_gaussian", "lorentzian", "p16_p50_p84__p50")))
add("recGE", f2(g(3, "single_gaussian", "exponential", "p16_p50_p84__p50")))
add("nullGG", f3(g(1, "single_gaussian", "single_gaussian", "p16_p50_p84__p50")))
add("nullEG", f3(g(1, "exponential", "single_gaussian", "p16_p50_p84__p50")))
hits = sum(
    (g(1, i, "single_gaussian", "count_le_observed") for i in ("single_gaussian", "exponential"))
)
den = sum((g(1, i, "single_gaussian", "denominator") for i in ("single_gaussian", "exponential")))
add("nullHits", "none" if hits == 0 else int(hits))
add("nullDen", int(den))
add("simA", "0.10")
add("hitEA", 16)
sec("Baccus & Xu table")
for n, v in (
    ("NX", 227),
    ("NXLRD", 48),
    ("NXCtrl", 97),
    ("NXOvRec", 29),
    ("NXOvSrc", 25),
    ("NXLabDiff", "six"),
    ("XWcut", 1000),
    ("XSNcut", 5),
    ("XkLRD", "5.4"),
    ("XkOth", "9.8"),
    ("XzCut", 3),
    ("XLmin", "43.7"),
    ("XLmax", "45.2"),
    ("NXHbCand", 13),
):
    add(n, v)
bm = r2["bx_main"]
assert bm["LRD"]["N"] == 48 and bm["other"]["N"] == 97
trip("Xsig", bm["LRD"]["sigma"], f3)
trip("Xb", bm["LRD"]["slope"], f2)
trip("Xsigo", bm["other"]["sigma"], f3)
trip("Xdsig", bm["diff"], lambda v: neg(f3(v)))
bn = r2["bx_noOverlap"]
add("NXNoOv", bn["LRD"]["N"])
add("NXNoOvCtrl", bn["other"]["N"])
trip("XsigNoOv", bn["LRD"]["sigma"], f3)
trip("XsigoNoOv", bn["other"]["sigma"], f3)
trip("XdsigNoOv", bn["diff"], lambda v: neg(f3(v)))
for n, v in (
    ("NXPairs", 36),
    ("XdMatch", "\\ensuremath{-}0.065"),
    ("XdMatchlo", "\\ensuremath{-}0.105"),
    ("XdMatchhi", "\\ensuremath{-}0.024"),
    ("NPermHits", "eight"),
    ("NPermTrials", "10^{6}"),
    ("NXLeave", 87),
    ("matchDz", "0.5"),
    ("matchDL", "0.2"),
):
    add(n, v)
add("XresM", f2(r2["bx_Mvir_LRD"]["resid_sd"]))
add("XresMo", f2(r2["bx_Mvir_ctrl"]["resid_sd"]))
sec("SDSS DR7 quasars")
w6 = sdss_cat
add("NSDSS", w6["N"])
trip("Ssig", w6["sigma"], f3)
add("Sb", f2(w6["slope"][2]))
add("sigRatio", "four")
add("SsigTwo", f2(w6["sigma"][2]))
add("FsigoTwo", f2(m["other"]["sigma"][2]))
add("sigTwoLo", f2(min(bm["other"]["sigma"][2], m["other"]["sigma"][2])))
add("sigTwoHi", f2(max(bm["other"]["sigma"][2], m["other"]["sigma"][2])))
sec("QSO1 diagnostic")
for n, v in (
    ("QSOsn", "3.4"),
    ("QSOW", 2107),
    ("QSOWerr", 302),
    ("QSOmu", "7.5"),
    ("QSOmuErr", "1.5"),
    ("QSOresCor", "0.17"),
    ("QSOresErr", "0.08"),
):
    add(n, v)
sec("photosphere: sources with sigma(log L_ph) < 0.3 dex")
add("vKnorm", 1493)
add("kappaES", "0.398")
A, B, C = (r6["A"], r6["B"], r6["C"])
E = r11["main_E"]
Ao = r11["Curtis_A_own"]
Bo = r11["Curtis_B_own"]
cx = {
    "PA": "X: ratio 2.9927, HPD 0.8310, b 0.2075 [0.0993, 0.3146], sigma 0.0732 [0.0483, 0.1168]",
    "PB": "X: ratio 2.6126, HPD 0.9042, b 0.1452 [0.0251, 0.2610], sigma 0.0812 [0.0507, 0.1370]",
    "PC": "X: ratio 2.7475, HPD 0.9846, b 0.4007 [0.2128, 0.5662], sigma 0.1146 [0.0669, 0.2310]",
}
for p, r in (("PA", A), ("PB", B), ("PC", C)):
    add("N" + p, "nine" if r["N"] == 9 else r["N"])
    add("rat" + p, f2(r["ratio_median"]))
    trip("b" + p, r["b"], f2)
    trip("c" + p, r["c"], lambda v: neg(f2(v)))
    add("hpd" + p, f"{100 * r['hpd']:.0f}")
    trip("sq" + p, r["sigma_fixed"], f3)
add("NPdrop", "five")
assert len(A["dropped"]) == 5 and len(B["dropped"]) == 3
add("NPAall", 23)
add("NPBall", 18)
add("sxCut", "0.3")
add("sxDropMin", "0.4")
add("sxDropMax", "2.0")
add("sxKeepMax", "0.18")
add("hierDb", "0.01")
assert E["N"] == 5 and Ao["N"] == 13 and (Bo["N"] == 15)
add("NPE", "five")
add("ratPE", f2(E["ratio_p16_50_84"][1]))
trip("sqPE", E["sigma_fixed_photosphere"], f3)
add("NPCbx", "four")
add("NdGall", 146)
add("NdGbb", 47)
add("NdGfit", 14)
add("bPElo", "\\ensuremath{-}0.7")
add("bPEhi", "0.7")
rats = [A["ratio_median"], B["ratio_median"], E["ratio_p16_50_84"][1]]
rl, rh = (min(rats), max(rats))
add("ratLo", f"{rl:.1f}")
add("ratHi", f"{rh:.1f}")
add("sqLo", f2(min(A["sigma_fixed"][2], B["sigma_fixed"][2])))
add("sqHi", f2(max(A["sigma_fixed"][2], B["sigma_fixed"][2])))
up = [100 * (10 ** r_["sigma_fixed"][2] - 1) for r_ in (A, B)]
add("uniPct", f"{5 * round(max(up) / 5):.0f}")
add("sigLamRLo", f2(2 * min(A["sigma_fixed"][2], B["sigma_fixed"][2])))
add("sigLamRHi", f2(2 * max(A["sigma_fixed"][2], B["sigma_fixed"][2])))
add("TminP", int(round(min(A["T_range"][0], C["T_range"][0]), -2)))
add("TmaxP", int(round(max(A["T_range"][1], C["T_range"][1]), -2)))
add("LphMin", f"{A['logL_range'][0]:.1f}")
add("LphMax", f"{A['logL_range'][1]:.1f}")
hb = [
    FS[k]["scenarios"]["measured"]["methods"]["flat"]["Hline_bc"][q + "_distribution"]["values"][2]
    for k in ("Curtis_A_own", "Curtis_B_own")
    for q in ("b", "c")
]
add("lineB", f"{(hb[0] + hb[2]) / 2:.1f}")
add("lineCLo", f"{min(hb[1], hb[3]):.1f}")
add("lineCHi", f"{max(hb[1], hb[3]):.1f}")
cod = lambda k, m: FS[k]["scenarios"]["measured"]["methods"][m]
for p, r, k in (("JE", E, "main_E"), ("JA", Ao, "Curtis_A_own"), ("JB", Bo, "Curtis_B_own")):
    f, h = (cod(k, "flat"), cod(k, "hierarchical"))
    sm = r["sim"]["hierarchical_baseline"]
    so = r["sim"]["own_flat_baseline"]
    note = "X flat %.4f, hierarchical %.4f" % (f["observed_b2"], h["observed_b2"])
    add("N" + p, {5: "five", 13: 13, 15: 15}[r["N"]])
    trip("bp" + p, r["partial"]["b2"], lambda v: neg(f2(v)))
    fr = [sm["p_ph"], so["p_ph"], f["p_ph"]["fraction"], h["p_ph"]["fraction"]]
    add("simPct" + p, f"{100 * sum(fr) / 4:.0f}")
    add("pow" + p, f"{100 * f['power']['fraction']:.0f}")
    add("simMed" + p, f"{sm['Hph_p16_50_84'][1]:.1f}" if p == "JE" else f2(sm["Hph_p16_50_84"][1]))
    trip("sl" + p, r["line_only"]["sigma"], f3)
    trip("so" + p, r["sigma_fixed_photosphere"], f3)
p10 = [
    100 * FS["main_E"]["scenarios"]["0.10"]["methods"][m_]["power"]["fraction"]
    for m_ in ("flat", "hierarchical")
]
add("powTen", f"{5 * round(sum(p10) / 10):.0f}")
add("powLo", f"{5 * math.floor(min((E['sim'][g]['power'] for g in E['sim'])) * 20):.0f}")
add("powHi", f"{5 * math.ceil(max((E['sim'][g]['power'] for g in E['sim'])) * 20):.0f}")
add("kRing", "2")
add("kIso", "1.36")
add("kVir", "0.96")
lam = lambda k, r: (k / r) ** 2
add("lamVirLo", f2(lam(0.964, rh)))
add("lamVirHi", f2(lam(0.964, rl)))
add("lamIsoLo", f2(lam(1.36, rh)))
add("lamIsoHi", f2(lam(1.36, rl)))
add("lamRingLo", f2(lam(2, rh)))
add("lamRingHi", f2(lam(2, rl)))
add("mEddLo", f"{(rl / 2) ** 2:.1f}")
add("mEddHi", f"{(rh / 0.964) ** 2:.0f}")
sec("membership and rules")
inc = load("membership.json")
m3 = json.load(open(DATA / "sample" / "rule_variants.json"))["variants"]["chi2_lt3"]["LRD"]
add("sigChiThree", f3(m3["posterior"]["intrinsic_scatter_dex"][2]))
cx = {
    r_["macro"]: float(r_["value"])
    for r_ in csv.DictReader(open(RESULTS / "membership_full_covariance.tsv"), delimiter="\t")
    if r_["macro"].startswith("sigInc")
}
assert len(cx) == 12
for p, k in (
    ("IncA", "rules_relaxed_same_labels"),
    ("IncB", "any_label_all_rules"),
    ("IncC", "any_label_chi2lt3"),
    ("IncD", "any_label_all"),
):
    r = inc[k]
    add("N" + p, r["N"])
    q = [cx["sig" + p + "lo"], None, cx["sig" + p], None, cx["sig" + p + "hi"]]
    for i in (0, 2, 4):
        assert abs(q[i] - r["sigma"][i]) < 0.003, (p, i)
    trip("sig" + p, q, f3)
assert inc["rules_relaxed_same_labels"]["N"] == 31
xd = [
    e
    for e in inc["extras"]
    if e["id"] in set(inc["any_label_all"]["added"]) - set(inc["any_label_chi2lt3"]["added"])
]
assert len(xd) == 2
add("chiHiA", f"{min((e['chi2'] for e in xd)):.1f}")
add("chiHiB", f"{max((e['chi2'] for e in xd)):.0f}")
add("NIncRules", "two")
sec("comparison with SDSS in like measures")
tg = r1["width_double_gaussian"]["LRD"]["sigma"][2]
npf = r1["width_nonparametric_FWHM"]["LRD"]["sigma"][2]
sd = sdss_cat["sigma"][2]
add("DsigTgTwo", f2(tg))
add("sdssRatLo", f"{sd / tg:.1f}")
add("sdssRatNp", f"{sd / npf:.1f}")
add("sdssRatHi", f"{sd / m['LRD']['sigma'][2]:.1f}")
sec("ratio to v_K by measure of width")
wr = load("width_definition_ratios.json")
for p, k in (("B", "Curtis_B_own"), ("E", "main_E")):
    for q, mm in (
        ("Np", "nonparametric_FWHM"),
        ("Tg", "double_gaussian"),
        ("Lo", "lorentzian"),
        ("Ex", "exponential"),
    ):
        add("rat" + q + p, f2(wr[k][mm]["ratio_median"]))
        if p == "B":
            trip("sq" + q + p, wr[k][mm]["sigma"], f3)
four = [
    wr["Curtis_B_own"][mm]["sigma"][2]
    for mm in ("single_gaussian", "nonparametric_FWHM", "lorentzian", "exponential")
]
add("sqMeasLo", f2(min(four)))
add("sqMeasHi", f2(max(four)))
add(
    "ratExLo",
    f"{min((wr[k]['exponential']['ratio_median'] for k in ('Curtis_B_own', 'main_E'))):.1f}",
)
add(
    "ratWingHi",
    f"{max((wr[k][mm]['ratio_median'] for k in ('Curtis_B_own', 'main_E') for mm in ('lorentzian', 'exponential'))):.1f}",
)
gfr = lambda rat, x: math.erfc(x / (rat / (2 * math.sqrt(2 * math.log(2)))) / math.sqrt(2))
efr = lambda rat, x: math.exp(-x / (rat / (2 * math.log(2))))
rats3 = [A["ratio_median"], B["ratio_median"], E["ratio_p16_50_84"][1]]
rex = wr["Curtis_B_own"]["exponential"]["ratio_median"]
add("fKepLo", f"{100 * gfr(min(rats3), 1):.0f}")
add("fKepHi", f"{100 * gfr(max(rats3), 1):.0f}")
add("fEscLo", f"{100 * gfr(min(rats3), math.sqrt(2)):.0f}")
add("fEscHi", f"{100 * gfr(max(rats3), math.sqrt(2)):.0f}")
add("fKepEx", f"{100 * efr(rex, 1):.0f}")
add("fEscEx", f"{100 * efr(rex, math.sqrt(2)):.0f}")
add("lamExpFac", f"{(B['ratio_median'] / rex) ** 2:.1f}")
sec("selection checks added after the adversarial review of 2026-11-13")
sn = load("snr_threshold.json")
for p, k in (("SNten", "LRD_SNge10"), ("SNtwenty", "LRD_SNge20")):
    add("N" + p, {20: 20, 9: "nine"}[sn[k]["N"]])
    trip("sig" + p, sn[k]["sigma"], f3)
flo = fm["LRD"]["a"] + m["LRD"]["slope"][2] * (d["LRD"]["L"][0] - 43) - 3
add("floorDex", f2(flo))
add("floorSig", f"{flo / m['LRD']['sigma'][2]:.1f}")
add("NTay", 13)
tp = load("photometric_labels_grid.json")
add("FsigTay", f3(tp["medians"]["FsigTay"]["finer_stored"]))
add("FsigTaylo", f3(tp["LRD_interval"]["upstream_finer"][0]))
add("FsigTayhi", f3(tp["LRD_interval"]["upstream_finer"][1]))
add("FsigoTay", f3(tp["medians"]["FsigoTay"]["finer_stored"]))
add("rhoThin", f"{math.sqrt(1 - (m['LRD']['sigma'][2] / r2['orient']['thin_0_45']) ** 2):.2f}")
add("kappaSolar", "0.34")
add("vKdrop", f"{100 * (1 - math.sqrt(0.34 / 0.398)):.0f}")
rlo = wr["Curtis_B_own"]["lorentzian"]["ratio_median"]
lfr = lambda rat, x: 1 - 2 / math.pi * math.atan(2 * x / rat)
add("fKepLor", f"{100 * lfr(rlo, 1):.0f}")
add("fEscLor", f"{100 * lfr(rlo, math.sqrt(2)):.0f}")
sec(
    "SDSS refitted with our single Gaussian; flux fractions by direct integration; hierarchical check"
)
add("NSsgDraw", 300)
add("NSsg", 224)
add("Ssg", "0.123")
add("Ssglo", "0.112")
add("Ssghi", "0.136")
add("SsgTwo", "0.12")
add("SsgCat", "0.185")
add("sdssRatSame", f"{0.12324 / m['LRD']['sigma'][2]:.1f}")
add("fKepNpLo", "38")
add("fKepNpHi", "44")
add("fEscNpLo", "23")
add("fEscNpHi", "28")
add("fEscBicLo", "24")
add("fEscBicHi", "28")
add("hierSig", "0.054")
add("hierB", "0.17")
sec("budget")
s0 = m["LRD"]["sigma"][2]
s1 = m["LRD"]["sigma"][4]
o = r2["orient"]
add("sigLamMax", f2(2 * s0))
add("sigLamMaxUp", f2(2 * s1))
add("oriThin", f2(o["thin_0_45"]))
add("oriHalf", f3(o["h0.5_0_60"]))
add("oriSeven", f3(o["h0.7_0_60"]))
add("oriOne", f3(o["h1_0_60"]))
add("oriThinRatio", "four")
add("sigLamIso", f2(2 * math.sqrt(s0**2 - o["h1_0_60"] ** 2)))
add("RLlo", "0.13")
add("RLhi", "0.19")
add("RLtermLo", f3(0.065))
add("RLtermHi", f3(0.095))
add("virW", "2.06")
add("virL", "0.47")
add("sigMvir", f2(2.06 * s0))
sec("scattering")
add("rusA", 428)
add("rusB", 370)
D_ = lambda t: 428 * t / (428 * t + 370)
add("DHalf", f2(D_(0.5)))
add("DOne", f2(D_(1)))
add("DTwo", f2(D_(2)))
add("sigTauHalf", f2(s0 / D_(0.5)))
add("sigTauOne", f2(s0 / D_(1)))
add("sigTauTwo", f2(s0 / D_(2)))
add("sigTeMax", f2(2 * s0))
add("tauSlopeLo", f"{m['LRD']['slope'][2] / D_(1):.1f}")
add("tauSlopeHi", f"{bm['LRD']['slope'][2] / D_(1):.1f}")
C11 = load("shell_radius_budget.json")["rows"]
for nm_, t_, c_ in zip(("sigRHalf", "sigROne", "sigRTwo", "sigRThree"), (0.5, 1, 2, 3), C11):
    v_ = 2 * s0 / D_(t_)
    assert (
        c_["tau_e"] == t_
        and abs(v_ - c_["sigma_logR_local"]) < 5e-05
        and (f2(v_) == c_["rounded_2dec"])
    ), (nm_, v_, c_)
    add(nm_, f2(v_))
add("tauFixLo", f2(0.5 * D_(0.5)))
add("tauFixHi", f2(0.5 * D_(2)))
add("tauPhLo", f2(0.25 * D_(0.5)))
add("tauPhHi", f2(0.25 * D_(2)))

# write only the macros that the manuscript uses, grouped under their section titles
source = re.sub(r"(?<!\\)%.*", "", (PAPER / "ms.tex").read_text(encoding="utf-8"))
values = {name: value for name, value in M if name != "%"}
used = {name for name in values if re.search(r"\\%s(?![A-Za-z])" % name, source)}
lines = [
    "%% numbers.tex -- generated by scripts/make_numbers.py; do not edit by hand.",
    "%% Write every macro as \\name{} in running text so that TeX does not swallow the following space.",
]
title = None
for name, value in M:
    if name == "%":
        title = value
    elif name in used:
        if title:
            lines += ["", f"%% {title}"]
            title = None
        lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
(PAPER / "numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
json.dump(values, open(RESULTS / "numbers.json", "w"), indent=1, ensure_ascii=False)
print(
    f"{len(values)} numbers computed; {len(used)} used in ms.tex and written to paper/numbers.tex"
)
