"""Dependence of the scatter on the sample rules: sources with conflicting labels, the reduced-chi-square cut, and the
absorption-bound rule.

Sources added to the main sample use their single-Gaussian fit with conditional errors and zero width-luminosity
covariance; data/results/membership_full_covariance.tsv holds the values with the covariance of each fit, which are
the ones quoted in the paper.

Writes data/results/membership.json.
"""

import collections
import csv
import json
import math

import numpy as np

from lrdwidth import stats as S
from lrdwidth.paths import DATA, RESULTS, WIDTHS

f = lambda s: float(s) if s not in ("", "nan", "None") else float("nan")
r = [
    x
    for x in csv.DictReader(open(DATA / "linefits" / "uniform_linefits.tsv"), delimiter="\t")
    if x["line"] == "Halpha"
]
aud = {
    a["file_id"]: a
    for a in json.load(open(DATA / "sample" / "selection_audit.json"))["Ha_file_quality_audit"]
}
W1 = json.load(open(WIDTHS / "regression_inputs.json"))["single_gaussian"]
fro = [t for t in W1 if t["cls"] == "LRD"]
froz_groups = {
    x["physical_source_group"]
    for x in csv.DictReader(open(DATA / "sample" / "frozen_sources.tsv"), delimiter="\t")
}
# file names of the candidate spectra of sources listed by de Graaff et al.
dgfiles = {c["file"] for c in json.load(open(DATA / "sample" / "candidates.json"))}
core = (
    "broad_SNR_ge5",
    "delta_BIC_null_ge6",
    "width_gt_LSF",
    "broad_width_not_at_bound",
    "optimizer_converged",
    "not_fitted",
)
G = collections.defaultdict(list)
for x in r:
    G[x["physical_source_group"]].append(x)
ex = []
for g, xs in G.items():
    if g in froz_groups:
        continue
    cand = []
    for x in xs:
        fr = aud[x["file_id"]]["failure_reasons"]
        if any((any((c in q for c in core)) for q in fr)):
            continue
        cand.append((f(x["main_broad_conditional_SNR"]), x, fr))
    if cand:
        cand.sort(key=lambda t: -t[0])
        sn, x, fr = cand[0]
        indg = any((y["file"] in dgfiles for y in xs))
        ex.append(
            dict(
                id=x["source_id"],
                cls=x["lrd_class_stage2"],
                dis=x["classification_disagreement"] == "True",
                zc=x["redshift_conflict"] == "True",
                dg=indg,
                sn=sn,
                chi2=f(x["main_reduced_chi2"]),
                fr=fr,
                x=math.log10(f(x["main_Halpha_broad_apparent_model_luminosity_Planck18_erg_s"]))
                - 43,
                y=math.log10(f(x["main_FWHM_intrinsic_kms"])),
                sy=f(x["main_FWHM_err_conditional_kms"])
                / f(x["main_FWHM_intrinsic_kms"])
                / math.log(10),
                sx=f(x["main_broad_emitted_flux_err_conditional_erg_s_cm2"])
                / f(x["main_broad_emitted_flux_erg_s_cm2"])
                / math.log(10),
            )
        )


def fit(extra):
    x = np.array([t["x"] for t in fro] + [e["x"] for e in extra])
    y = np.array([t["y"] for t in fro] + [e["y"] for e in extra])
    cxx = np.array([t["cov"][0][0] for t in fro] + [e["sx"] ** 2 for e in extra])
    cyy = np.array([t["cov"][1][1] for t in fro] + [e["sy"] ** 2 for e in extra])
    cxy = np.array([t["cov"][0][1] for t in fro] + [0.0 for e in extra])
    b, s, p, ah = S.post2(x, y, cxx, cyy, cxy)
    qs = S.quant(s, p.sum(0))
    qb = S.quant(b, p.sum(1))
    return dict(N=len(x), sigma=qs, slope=qb, added=[e["id"] for e in extra])


out = {}
out["frozen"] = fit([])
isL = lambda e: e["cls"] == "LRD" and (not e["dis"]) and (not e["zc"])
out["rules_relaxed_same_labels"] = fit([e for e in ex if isL(e)])
anyL = lambda e: (e["cls"] == "LRD" or e["dis"] or e["dg"]) and (not e["zc"])
out["any_label_all_rules"] = fit([e for e in ex if anyL(e) and (not e["fr"])])
out["any_label_chi2lt3"] = fit([e for e in ex if anyL(e) and e["chi2"] < 3])
out["any_label_all"] = fit([e for e in ex if anyL(e)])
out["extras"] = ex
json.dump(out, open(RESULTS / "membership.json", "w"), indent=1)
