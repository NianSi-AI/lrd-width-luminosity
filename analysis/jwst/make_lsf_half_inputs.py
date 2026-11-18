"""Regression inputs of the halved-LSF variant (data/widths/regression_inputs_lsf_half.json; lsf_half_inputs.json
is an identical copy) from the per-source single-Gaussian Halpha fits at LSF FWHM x0.5.

Fit records are read either from the stored data/widths/source_fits_lsf_half/ (default) or from a directory written
by `python analysis/jwst/refit_spectra.py --spectra work/jwst_spectra --frozen-only --line Halpha --lsf-scales 0.5
--out work/jwst_lsf_half` (pass --fits work/jwst_lsf_half). Membership is the frozen sample; the quality checks are
those of paper Section 2.3 (analysis/sample/make_selection_audit.py) and are reported, not used to drop sources.

    python analysis/jwst/make_lsf_half_inputs.py                 # compare with the stored file
    python analysis/jwst/make_lsf_half_inputs.py --out FILE.json
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18

from lrdwidth.paths import DATA, ROOT, WIDTHS

sys.path.insert(0, str(ROOT / "analysis" / "sample"))
from make_selection_audit import quality_tests  # noqa: E402

STORED = WIDTHS / "regression_inputs_lsf_half.json"
# key names of the stored rows (as in the injection code) for the Section 2.3 tests
CHECK_NAMES = [("broadSNR", "broad_SNR_ge5"), ("deltaBIC", "delta_BIC_null_ge6"), ("widthLSF", "width_gt_LSF"),
               ("widthBoundary", "broad_width_not_at_bound"), ("chi2", "reduced_chi2_below_threshold"),
               ("converged", "optimizer_converged"), ("absBoundary", "adopted_absorption_parameters_not_at_bounds")]


def lsf_half_main(path):
    a = json.loads(path.read_text())
    if "main" in a:  # stored record
        return a["main"]
    ha = next(x for x in a["lines"] if x["line"] == "Halpha")  # lrdwidth.fit_batch output
    return next(c for c in ha["cases"] if abs(c["lsf_FWHM_scale"] - 0.5) < 1e-9)["main"]


def build(fits_dir):
    frozen = list(csv.DictReader((DATA / "sample" / "frozen_sources.tsv").open(), delimiter="\t"))
    old_y = {r["file_id"]: r["y"] for r in json.loads((WIDTHS / "regression_inputs.json").read_text())["single_gaussian"]}
    rows = []
    for r in frozen:
        m = lsf_half_main(fits_dir / f"{r['file_id']}.json")
        names, p = m["names"], np.array(m["params"])
        ii = [names.index("F_broad1"), names.index("W_broad1")]
        J = np.diag(1 / (np.log(10) * p[ii]))
        cov = J @ np.array(m["covariance_conditional"])[np.ix_(ii, ii)] @ J.T
        L = m["broad_flux_erg_s_cm2"] * 4 * np.pi * Planck18.luminosity_distance(float(r["z_fit"])).to_value("cm") ** 2
        tests = quality_tests(m, m["LSF_FWHM_kms"], 2.0)
        old = 10 ** old_y[r["file_id"]]
        rows.append({
            "source_id": r["source_id"], "file_id": r["file_id"], "cls": r["lrd_class_stage2"],
            "x": float(np.log10(L) - 43), "y": float(np.log10(m["FWHM_intrinsic_kms"])), "cov": cov.tolist(),
            "width_kms": m["FWHM_intrinsic_kms"], "luminosity_erg_s": float(L), "flux_erg_s_cm2": m["broad_flux_erg_s_cm2"],
            "quality_checks": {k: bool(tests[t]) for k, t in CHECK_NAMES}, "quality_pass": bool(all(tests.values())),
            "delta_logwidth_dex": float(np.log10(m["FWHM_intrinsic_kms"]) - np.log10(old)), "old_FWHM_kms": float(old),
            "absorption": m["absorption"], "bound_parameters": m["bound_parameters"],
        })
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--fits", type=Path, default=WIDTHS / "source_fits_lsf_half")
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    text = json.dumps(build(a.fits), indent=2)
    if a.out:
        a.out.write_text(text)
        print("written", a.out)
    elif text == STORED.read_text():
        print(f"{STORED.relative_to(ROOT)} reproduced byte for byte from {a.fits}")
    else:
        sys.exit(f"MISMATCH with {STORED.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
