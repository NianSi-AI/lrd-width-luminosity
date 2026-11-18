# Data

Line luminosities are in erg/s, widths in km/s, covariances in dex² (base-10 logarithms of luminosity and width).
In the regression inputs `x` is log L − 43 for Hα luminosities and log L − 45 for bolometric or photospheric
luminosities, `y` is log FWHM, and `t` is log(T_eff / 5000 K). "Frozen" marks the main sample of 44 sources (29 LRDs,
15 other broad-line sources), which was fixed before the analysis.

| Directory | Contents |
|---|---|
| `linefits/` | `uniform_linefits.tsv`: single-Gaussian fits of every line in every spectrum of the parent sample |
| `sample/` | The frozen sample (`frozen_sources.tsv`, `table1_sources.csv`), the audit of the quality rules for every spectrum (`selection_audit.json`, rebuilt from `linefits/uniform_linefits.tsv` by `analysis/sample/make_selection_audit.py`), fits under alternative rules (`rule_variants.json`), and the candidates added in the membership tests (`candidates.json`, `candidate_source_table.json`, `candidate_fits/`) |
| `widths/` | Six measures of width for the frozen sample: per-source table (`source_widths_six_definitions.tsv`), regression inputs (`regression_inputs.json`; `regression_inputs_lsf_half.json` with the instrumental width halved, rebuilt from `source_fits_lsf_half/` by `analysis/jwst/make_lsf_half_inputs.py`; `lsf_half_inputs.json` is an identical copy), per-source fit records (`source_fits/`, `source_fits_lsf_half/`), profile comparison by BIC, fit-window (±6800 and ±8000 km/s) and wing-extrapolation summaries (`window_and_wing_checks.json`; `B15_nonparam_flux_below_vesc_pct_median` is the median percentage of nonparametric broad flux below the escape speed in sample B), covariances of the parent and repeat spectra, and the injection–recovery summaries |
| `bx26/` | Regression inputs built from the machine-readable table of Baccus & Xu (2026), and its overlap with our sample |
| `jwst/` | Spectrum inventory with URLs and SHA-256 (`targets.json`), released main-fit parameters and hashes of the 44 injection templates (`injection_parameters.json`), the 20 in-band QC-failing parent templates (`failed_templates_20valid.json`), and the pre-registered criteria of the 49-template suite (`parent49_criteria.json`) |
| `qso1/` | Fit of A2744-QSO1 with the procedure of the paper |
| `photosphere/` | Sources with published photospheric parameters: atmosphere fits of Curtis et al. (2026) with BX26 widths (`atmosphere_bx26_widths.tsv`, sample A) and with our widths (`atmosphere_our_widths.tsv`, sample B), blackbody fits of de Graaff et al. (2026) (`blackbody_sources.tsv`, sample C); inputs and results of the simulations |
| `sdss/` | SDSS DR7 quasars refitted with the single-Gaussian procedure: the random sample of 300, the spectra used, one fit record per source (`fits/`), the measurement table and the regression |
| `tables/` | Machine-readable versions of Table 1 (44 sources × six widths, with covariances) and of the photospheric samples, in CSV and ECSV |
| `results/` | Output of the scripts in `analysis/`, plus `numbers.json` (every number in the manuscript) |

Other files in `results/`:

- `sdss_catalog_benchmark.json`: scatter of the SDSS DR7 quasars from the catalogue widths of Shen et al. (2011)
  (`analysis/sdss/catalog_benchmark.py`; the catalogue is not redistributed).
- `membership_full_covariance.tsv`: scatter of the inclusive samples with the covariance of each fit.
- `measurement_summary.json`: summary statistics of the line-profile measurements (profile counts, flux fractions,
  injection–recovery medians), keyed by the name of the quantity.
- `photometric_labels_grid.json`: scatter for the sample defined by photometric labels alone.
- `shell_radius_budget.json`: the shell-radius scatter budget of Section 5.3 and its checks.
- `scattering_kernel.json`: broadening and optical-depth sensitivity for Gaussian and top-hat cores convolved with an exponential scattering kernel (`analysis/10_scattering_kernel.py`).
- `repeat_calibration.json`: empirical per-spectrum error floors from the 9 multi-spectrum LRDs, posterior refits with floors, one-sided 95% upper limits, BX26 sky deduplication, and refitted SDSS width offset (`analysis/11_repeat_calibration.py`).
- `ancillary_checks.json`: BIC-preferred width regression, MAD and leave-one-out checks, and wide-window sensitivity summary (`analysis/12_ancillary_checks.py`).
- `parent_covariance_audit.json`, `repeat_random_effects.json`, `exterior_potential_check.json`: in-band parent-sample regressions (`N = 29..39`), repeat-pair reweightings, 43-spectrum bivariate random-effects models (additive, multiplicative, and free G395H offset), and exterior-envelope potential bounds (`analysis/13_parent_and_repeat_checks.py`).
- `review_checks.json`: 19-variant summary with one-sided 95% upper limits, [O III]-anchored refit with the absorption rule re-applied, constant-width model, photospheric partial slope without RUBIES-UDS-29813, error-inflation tests, grating and primary-spectrum offsets, SDSS without the chi^2 rule, and shell-slope predictions (`analysis/14_review_checks.py`).
- `narrow_degeneracy_refit.json`, `narrow_refit_11sources.tsv`, `frozen_fixed_narrow_audit.json`, `parent49_injection.json`: narrow-core refits of the parent and frozen samples and the 49-template parent injection–recovery suite (`analysis/jwst/`, `analysis/injection/parent49_recovery.py`).
