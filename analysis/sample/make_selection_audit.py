"""Rebuild data/sample/selection_audit.json and the frozen sample from the released per-spectrum Halpha fits.

Input: the single-Gaussian Halpha main fit (LSF x1) of every spectrum in data/linefits/uniform_linefits.tsv
(broad S/N, Delta BIC against the no-broad null, intrinsic FWHM, default LSF FWHM, parameters at bounds, adopted
absorption, reduced chi^2, convergence). Rule of paper Section 2.3, as in lrdwidth.fit_batch.quality_tests:
bounds reject a fit only for the broad width and, if the absorption screen is adopted, for tau_abs or W_abs.
Then: highest broad S/N per physical source, then exclusion of missing/conflicting class labels and redshift
conflicts.

    python analysis/sample/make_selection_audit.py                  # compare with the stored file (default)
    python analysis/sample/make_selection_audit.py --out FILE.json  # write the rebuilt file
    python analysis/sample/make_selection_audit.py --reduced-chi2-max 3   # counts of the chi^2 < 3 rule variant

The two SHA-256 fields of the stored file identify the original fit-summary table and rule note of the development
run, which are not in the repository; they are carried over as recorded, not recomputed.
"""
import argparse
import csv
import json
import math
import sys
from collections import Counter

from lrdwidth.paths import DATA

AUDIT = DATA / "sample" / "selection_audit.json"
BROAD_WIDTH_PARAMETERS = {"W_broad1", "W_broad2"}
ABSORPTION_BOUND_PARAMETERS = {"tau_abs", "W_abs"}
RULE = "Filter every file,then highest broadSNR per physicalsource,then class/redshift conflict exclusion."
RECORDED_SHA256 = {
    "input_table_sha256": "e6da80bd60bad58b1eaaef48f5bbe251ebcaca51a669d863f07ee1180bbb9f4f",
    "effective_rule_doc_sha256": "9ebc32538a43eac784c7cc7f612e5926839142bddb80817a6ac1eef687c490b2",
}


def num(s):
    return float(s) if s not in ("", "nan", "None") else math.nan


def fit_record(r):
    """Fields of the released main fit used by the quality rule (names as in lrdwidth.linefit/fit_batch)."""
    return {
        "broad_SNR_conditional": num(r["main_broad_conditional_SNR"]),
        "delta_BIC_null": num(r["main_delta_BIC_null"]),
        "FWHM_intrinsic_kms": num(r["main_FWHM_intrinsic_kms"]),
        "bound_parameters": json.loads(r["main_parameter_bounds"]),
        "success": r["main_optimizer_converged"] == "True",
        "chi2": 0.0,  # finite by construction for fitted rows; only convergence is recorded in the table
        "reduced_chi2": num(r["main_reduced_chi2"]),
        "absorption": r["main_absorption_adopted"] == "True",
    }


def quality_tests(f, lsf_fwhm_kms, chi2_max):
    """Same tests and order as lrdwidth.fit_batch.quality_tests (kept here so no calibration files are needed)."""
    bounds = set(f["bound_parameters"])
    return {
        "broad_SNR_ge5": (f["broad_SNR_conditional"] or 0) >= 5,
        "delta_BIC_null_ge6": f["delta_BIC_null"] >= 6,
        "width_gt_LSF": f["FWHM_intrinsic_kms"] > lsf_fwhm_kms,
        "broad_width_not_at_bound": not bounds & BROAD_WIDTH_PARAMETERS,
        "optimizer_converged": f["success"] and math.isfinite(f["chi2"]),
        "reduced_chi2_below_threshold": f["reduced_chi2"] < chi2_max,
        "adopted_absorption_parameters_not_at_bounds": not (f["absorption"] and bounds & ABSORPTION_BOUND_PARAMETERS),
    }


def build(chi2_max=2.0):
    with (DATA / "linefits" / "uniform_linefits.tsv").open() as fp:
        rows = [r for r in csv.DictReader(fp, delimiter="\t") if r["line"] == "Halpha"]
    assert len(rows) == 204
    audit, passing = [], []
    for r in rows:
        if r["status"] != "fitted":
            reasons = ["not_fitted"]
        else:
            tests = quality_tests(fit_record(r), num(r["default_LSF_FWHM_kms"]), chi2_max)
            reasons = [k for k, ok in tests.items() if not ok]
        audit.append({"file_id": r["file_id"], "physical_source_group": r["physical_source_group"],
                      "source_id": r["source_id"], "failure_reasons": reasons})
        if not reasons:
            passing.append(r)
    best = {}
    for r in passing:  # first spectrum wins a tie, as in file order
        g = r["physical_source_group"]
        if g not in best or num(r["main_broad_conditional_SNR"]) > num(best[g]["main_broad_conditional_SNR"]):
            best[g] = r
    clean = [r for r in best.values() if r["lrd_class_stage2"] in ("LRD", "non-LRD")
             and r["classification_disagreement"] != "True" and r["redshift_conflict"] != "True"]
    counts = {"N_pass_files": len(passing), "N_best_groups": len(best), "N_clean": len(clean),
              "by_class": dict(sorted(Counter(r["lrd_class_stage2"] for r in clean).items()))}
    out = {"user_decisions": {"reduced_chi2_threshold": int(chi2_max) if chi2_max == int(chi2_max) else chi2_max,
                              "exclude_adopted_tau_abs_or_W_abs_boundary": True},
           "counts": counts, "Ha_file_quality_audit": audit, **RECORDED_SHA256, "rule": RULE}
    return out, clean


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--out", help="write the rebuilt audit here instead of comparing with the stored file")
    p.add_argument("--reduced-chi2-max", type=float, default=2.0)
    a = p.parse_args()
    out, clean = build(a.reduced_chi2_max)
    text = json.dumps(out, indent=2)
    print("counts:", json.dumps(out["counts"]))
    if a.reduced_chi2_max == 2.0:
        with (DATA / "sample" / "frozen_sources.tsv").open() as fp:
            frozen = [r["file_id"] for r in csv.DictReader(fp, delimiter="\t")]
        same = frozen == [r["file_id"] for r in clean]
        print("selected spectra identical to data/sample/frozen_sources.tsv (same order):", same)
    if a.out:
        with open(a.out, "w") as fp:
            fp.write(text)
        print("written", a.out)
        return
    if a.reduced_chi2_max != 2.0:
        return
    stored = AUDIT.read_text()
    if text == stored:
        print(f"{AUDIT.relative_to(DATA.parent)} reproduced byte for byte")
        return
    old = json.loads(stored)
    diff = [x["file_id"] for x, y in zip(out["Ha_file_quality_audit"], old["Ha_file_quality_audit"]) if x != y]
    print("MISMATCH with stored audit; differing spectra:", diff, "counts stored:", old["counts"])
    sys.exit(1)


if __name__ == "__main__":
    main()
