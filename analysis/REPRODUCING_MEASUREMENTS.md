# Recomputing measurements and simulations

Run commands from the repository root after `pip install -e .` (or set `PYTHONPATH=.`).
All statistical analyses (`make analysis`, then `make paper`) reproduce from the stored
per-source measurements in `data/`. The commands below expose the original spectral and
simulation generators; they need the external DJA spectra (`scripts/fetch_jwst.py`) and
the msaexp dispersion tables (`scripts/fetch_calibration.sh`, pinned and SHA-256 checked),
and they retain the original likelihoods, conditional covariance, grids and seed sequences.
Code provenance is in `reproduction_provenance.json`; file names there identify the
development versions. Section 6 lists the steps that are not scripted here.

## Short validation

```sh
python scripts/check_reproduction.py
```

This checks a complete frozen recovery replicate against the corrected grid
posterior, the hierarchical log density against the independent NumPy target,
measurement-noise covariance, and all 767 computed macro values and 497 manuscript/table definitions.
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
Quality flags follow Section 2.3: a parameter at a bound rejects a fit only if it is
the broad width or, with absorption adopted, `tau_abs` or `W_abs` (the development
batch script flagged any bound; the stored selection already used this rule).
QC-passing LSF x1 fits get 200 nonparametric bootstrap draws (development script: 80),
the number of draws behind the stored nonparametric covariances (Section 6).

Sensitivity variants (one output directory each; records store their `fit_settings`):

```sh
S="--spectra work/jwst_spectra --frozen-only --line Halpha"
python analysis/jwst/refit_spectra.py $S --lsf-scales 0.5 --out work/jwst_lsf_half
python analysis/jwst/refit_spectra.py $S --window-kms 6800 --out work/jwst_win6800
python analysis/jwst/refit_spectra.py $S --window-kms 8000 --mask-sii --out work/jwst_win8000
python analysis/jwst/make_lsf_half_inputs.py --fits work/jwst_lsf_half --out lsf_half.json
```

`make_lsf_half_inputs.py` without arguments rebuilds `data/widths/regression_inputs_lsf_half.json`
(`lsf_half_inputs.json` is an identical copy) byte for byte from `data/widths/source_fits_lsf_half/`.
The [S II] mask half-width of the stored +-8000 km/s run is not recorded (`--mask-halfwidth-kms`,
default 500 km/s), and no script aggregates the window runs into `window_and_wing_checks.json`.

The quality audit and the frozen sample (29 LRDs, 15 other sources) follow from the
released fits in `data/linefits/uniform_linefits.tsv`:

```sh
python analysis/sample/make_selection_audit.py   # reproduces data/sample/selection_audit.json byte for byte
```

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

## 5. Parent-sample narrow-core refits and 49-template injection suite

```sh
python scripts/fetch_jwst.py --download --qso1     # 204 spectra plus A2744-QSO1
export LRDWIDTH_SPECTRA="$PWD/work/jwst_spectra"
python analysis/jwst/narrow_core_refit.py
python analysis/jwst/frozen_narrow_core_audit.py
python analysis/injection/parent49_recovery.py --smoke-test
python analysis/injection/parent49_recovery.py --work work/parent49_injection
```

They overwrite `narrow_degeneracy_refit.json`, `narrow_refit_11sources.tsv`,
`frozen_fixed_narrow_audit.json` and `parent49_injection.json` in `data/results/`, and
`data/jwst/failed_templates_20valid.json`; run them in this order. Apart from the spectra,
their inputs are in the repository: released Halpha fits (`data/linefits/uniform_linefits.tsv`),
main-fit parameters (`data/jwst/injection_parameters.json`), the audit of the 21 QC-failing
LRDs (`data/widths/parent_and_repeat_covariances.json`) and the A2744-QSO1 fit
(`data/qso1/qso1_fit.json`). Each script stops before fitting if a spectrum is missing or
differs from its recorded SHA-256; `narrow_core_refit.py` then checks that it reproduces the
released Skyfire-23289 and A2744-QSO1 fits, and `parent49_recovery.py` checks its 29-template
subset against `data/widths/sloped_injection_recovery.json`, whose generating script is not included.
For failed templates without a published broad Halpha luminosity, the injected luminosity
is the narrow-only (no-broad) Halpha flux of the fixed-narrow-width refit; no empirical
LRD draw is made (cf. `data/jwst/parent49_criteria.json`).

## 6. Not scripted in this repository

- Tabulation of the per-spectrum fit records into `data/linefits/uniform_linefits.tsv`,
  `data/sample/` (except `selection_audit.json`) and `data/widths/` (`source_fits/`,
  `regression_inputs.json`, `source_widths_six_definitions.tsv`, `profile_bic_sources.tsv`,
  `parent_and_repeat_covariances.json`, `window_and_wing_checks.json`), and the A2744-QSO1
  fit `data/qso1/qso1_fit.json`. The 200-draw nonparametric covariances in `source_fits/`
  use the algorithm of `analysis/injection/measure.py:np_bootstrap`; the run on the observed
  spectra and its seed are not included.
- `parent_and_repeat_covariances.json` cannot be rebuilt from stored per-fit data: the
  flux-width covariance is not stored for 3 of its 43 QC spectra and 11 of the 21 failed
  spectra, nor the narrow-line widths and in-band error ratios of the failed spectra.
- The 29-template sloped injection suites (`data/widths/sloped_injection_recovery.json`), and
  `membership_full_covariance.tsv`, `photometric_labels_grid.json` and `shell_radius_budget.json`
  in `data/results/`.
- SDSS refits (`analysis/sdss/refit_spectra.py`) need the 300 DR17 spectra in
  `data/sdss/spectra/` and the list `data/sdss/download_successes.json`, which are not
  included (URLs and SHA-256 are in `third_party/external_data.tsv`; each list row is the
  `input` record of the matching `data/sdss/fits/*.json`). `analysis/sdss/catalog_benchmark.py`
  needs the Shen et al. (2011) catalogue, which is not redistributed.
