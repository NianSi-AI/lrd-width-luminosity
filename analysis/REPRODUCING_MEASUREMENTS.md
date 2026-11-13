# Recomputing measurements and simulations

Run commands from the repository root after `pip install -e .` (or set `PYTHONPATH=.`).
The existing `make paper` path uses stored results. The commands below expose the
original spectral and simulation generators. They retain the original likelihoods,
selection, conditional covariance, grids and seed sequences. Code provenance is in
`reproduction_provenance.json`; file names there identify the development versions.

## Short validation

```sh
python scripts/check_reproduction.py
```

This checks a complete frozen recovery replicate against the corrected grid
posterior, the hierarchical log density against the independent NumPy target,
measurement-noise covariance, and all 459 macro values/332 manuscript definitions.
It does not run every production fit or all 24,000 hierarchical datasets.

## 1. JWST spectra and batch fits

DJA spectra remain external, with exact URLs and SHA256 in `data/jwst/targets.json`.
The file stores only source identity, redshift, path, URL and checksum, not the
third-party classification catalogue. Instrument tables are fetched separately.

```sh
bash scripts/fetch_calibration.sh
python scripts/fetch_jwst.py --limit 1              # list only
python scripts/fetch_jwst.py --download            # exact 204-spectrum inventory
python analysis/jwst/refit_spectra.py --spectra work/jwst_spectra --limit 1 --line Halpha
python analysis/jwst/refit_spectra.py --spectra work/jwst_spectra
```

Fits are written to `work/jwst/`; existing fit records are reused. The single-Gaussian
main width, alternative profiles, phenomenological absorption, no-broad null,
three LSF scales and nonparametric diagnostics follow the original batch routine.
The optional development pilot cache was removed, so every requested fit is
computed directly. Flux is in erg/s/cm^2; intrinsic widths are km/s.

## 2. Seeded injection recovery

`data/jwst/injection_parameters.json` contains the 44 author-fitted template
parameters and spectrum/error hashes. Pixel values are read from external spectra.

```sh
python scripts/fetch_jwst.py --download --frozen-only
export LRDWIDTH_SPECTRA="$PWD/work/jwst_spectra"
python analysis/injection/gaussian_recovery.py --limit 1  # small execution check
python analysis/injection/gaussian_recovery.py           # 2 profiles x 4 scenarios x 100 x 44
python analysis/injection/width_recovery.py --scenario 1
python analysis/injection/width_recovery.py --scenario 3
python analysis/injection/finalize.py --scenario 1
python analysis/injection/finalize.py --scenario 3
```

Products live under `work/injection`. Gaussian recoveries are prerequisites for
other definitions. Scenarios 0/1/2 inject log-width scatter 0.10/0.15/0.20 dex;
scenario 3 injects the luminosity relation with 0.05 dex scatter. The same
SeedSequence `[20261112, profile, scenario, replicate, template]` is retained.
The multi-definition experiments use the 29 LRD templates and scenarios 1/3.
The outputs `W2_final_summary_s1.json` and `W2_final_summary_s3.json` correspond
to `data/widths/injection_recovery_s1.json` and `injection_recovery_s3.json`.
These stored companions have been aligned to the final cell-boundary convention
already used by the manuscript; earlier summaries remain in Git history.
The finalizer retains failed definitions and fixed selection membership, with
at least 160/200 usable nonparametric bootstrap draws. It uses the final
cell-boundary CDF correction; empirical replicate quantiles remain NumPy linear.
`data/widths/injection_observed.json` holds the matching corrected observed
thresholds. Full runs are expensive; production denominators are not reduced by
using a small execution check.

## 3. Hierarchical fixed-design simulations

```sh
python analysis/hierarchical/run_simulations.py --selfcheck-only
python analysis/hierarchical/run_simulations.py --out work/hierarchical --limit 1 --workers 1
# Use a different output directory for a production run after a small check:
python analysis/hierarchical/run_simulations.py --out work/hierarchical_full --n 2000 --workers 6
python analysis/hierarchical/summarize.py work/hierarchical_full
```

The three datasets, generating baselines and observed thresholds are in
`data/photosphere/hierarchical/`. There are 3 samples x 4 scatter scenarios x
2 hypotheses x 2000 = 48,000 generated datasets, with multiple fitted statistics
per dataset. (The two-hypothesis comparison has 24,000 datasets per hypothesis.)
The output summary is `work/hierarchical_full/F_simulation_results.json`, the
source of `data/photosphere/simulation_results.json`. All failed fits and planned
denominators are retained; the final summary requires 2000 draws per cell.

Numba is optional (`pip install -e '.[simulations]'`) for acceleration. Without it
the same arithmetic runs in Python and can be much slower. Same seed does not
guarantee identical Monte Carlo posterior summaries across numerical-library and
JIT versions; compare diagnostics and sampling uncertainty. Stored production
results are retained, not replaced by short-check outputs. The design conditions
on observed predictor medians and does not simulate targeting or selection.

## 4. Measurement-summary construction

```sh
python scripts/make_measurement_summary.py
python scripts/make_numbers.py
```

The summary is now built from BIC counts and the **corrected** recovery aggregates
in `data/measurements/`, using exact JSON pointers and units in `summary_map.json`.
Only the 887 fields used by the manuscript generator are retained; unrelated
historic registry entries (including publication fees) have been removed.
For new injection runs, copy the regenerated `W2_final_summary_s1.json` and
`W2_final_summary_s3.json` into a separate inputs directory alongside
`W1_BIC_counts.json`, and pass `--inputs DIRECTORY --out OUTPUT.json`.
This step extracts aggregate measurements; it does not itself rerun spectral fits.
The original counts can also be rebuilt from the per-source profile BICs (see
`make_profile_counts.py`).
