# A tight broad-line width–luminosity relation in little red dots

Data, analysis code and manuscript sources for

> Haoyu Wang & Xinmiao Zhao, *A Tight Broad-line Width–Luminosity Relation in Little Red Dots and Its Consequence for
> Virial Black Hole Masses* (in preparation).

We refit public JWST/NIRSpec grating spectra of 29 little red dots (LRDs) and 15 other broad-line sources with one
procedure, measure the intrinsic scatter of the broad Hα width at fixed line luminosity, compare it with other
broad-line samples and with SDSS quasars fitted in the same way, and compare the widths with the Keplerian speed at
the photosphere of an Eddington-limited source.

## Layout

| Path | Contents |
|---|---|
| `lrdwidth/` | Python package: grid posteriors (`stats.py`), regressions (`regression.py`, `photosphere.py`), the emission-line fitting engine (`linefit.py`) |
| `analysis/` | Scripts that turn the per-source measurements in `data/` into the results in `data/results/` |
| `scripts/` | Scripts that turn the results into the numbers, tables and figures of the manuscript |
| `data/` | Per-source measurements, samples, and stored results; see [`data/README.md`](data/README.md) |
| `paper/` | AASTeX source of the manuscript |
| `third_party/` | Licenses of bundled third-party files and links to data that are not redistributed |

## Installation

Python 3.10 or later.

```sh
pip install -e .
```

This installs `numpy`, `scipy`, `matplotlib` and `astropy`. The versions used for the paper are pinned in
`requirements.txt`.

## Reproducing the paper

Everything quoted in the manuscript is regenerated from the stored results:

```sh
python scripts/make_numbers.py    # paper/numbers.tex: every number in the text, as a LaTeX macro
python scripts/make_tables.py     # paper/tables/
python scripts/make_figures.py    # paper/figures/
cd paper && pdflatex ms && bibtex ms && pdflatex ms && pdflatex ms
```

or `make paper`. The manuscript contains no hand-typed results: the text refers to macros in `numbers.tex`.

The stored results are regenerated from the per-source measurements with the scripts in `analysis/`
(`make analysis` runs all but `08_photosphere_simulations.py`, which runs 2000 simulations per sample and is slow):

| Script | Result file in `data/results/` | Contents |
|---|---|---|
| `01_width_luminosity.py` | `width_luminosity.json` | the relation for six measures of width; widths against the photospheric Keplerian speed |
| `02_robustness.py` | `robustness.json` | sample variants, halved instrumental width, the BX26 table, virial-mass residuals, orientation scatter |
| `03_joint_line_photosphere.py` | `joint_line_photosphere.json` | joint regression on line and photospheric luminosity (not quoted in the paper) |
| `04_joint_robustness.py` | `joint_robustness.json` | robustness of the joint regression (not quoted in the paper) |
| `05_photosphere_constrained.py` | `photosphere_constrained.json` | photospheric comparison for sources with constrained photospheric luminosity |
| `06_membership.py` | `membership.json` | dependence on the sample rules |
| `07_snr_threshold.py` | `snr_threshold.json` | dependence on the signal-to-noise threshold |
| `08_photosphere_simulations.py` | `photosphere_simulations.json` | partial regressions and simulations under the two hypotheses |
| `09_width_definition_ratios.py` | `width_definition_ratios.json` | ratio to the Keplerian speed for each measure of width |

Regenerated files agree with the committed ones to floating-point precision, so `git diff data/results` is the check.

## What is and is not included

Included: the line measurements of every source (six measures of width, with covariances), the samples and the record
of the selection rules, the regression inputs, and all derived results.

Not included, and needed only to repeat the measurements from the spectra:

- **Spectra.** The JWST spectra are the public reductions of the DAWN JWST Archive and the SDSS spectra are from DR17;
  `third_party/external_data.tsv` lists the address of every file. `scripts/fetch_calibration.sh` downloads the
  NIRSpec dispersion tables that `lrdwidth.linefit` needs.
- **Third-party catalogues without a redistribution license**, in particular the SDSS DR7 quasar catalogue of
  Shen et al. (2011) and the table of Kocevski et al. (2026). Results that depend on them
  (`data/results/sdss_catalog_benchmark.json`, the labels in `data/sample/`) are stored;
  `analysis/sdss/catalog_benchmark.py` recomputes the first from the catalogue.
The spectral batch drivers, seeded injection recovery, hierarchical simulations,
and measurement-summary construction are available under `analysis/jwst/`,
`analysis/injection/`, `analysis/hierarchical/`, and `scripts/`.
See [Recomputing measurements and simulations](analysis/REPRODUCING_MEASUREMENTS.md)
for inputs, commands, validation and limits. `python scripts/check_reproduction.py`
runs the short validation without downloading spectra or running production simulations.

## Use of AI tools

The code in this repository was written with the assistance of AI tools (Claude, Codex and ChatGPT) under the
authors' direction, as stated in the manuscript. The authors are responsible for its content.

## License and citation

Code: MIT (`LICENSE`). Measurements, tables, figures and documentation: CC BY 4.0 (`LICENSE-DATA.md`).
Third-party material keeps its own terms; see `ATTRIBUTION.md` and `third_party/`.

Please cite the paper and the data release, doi:10.5281/zenodo.24182569 (reserved; it resolves once the release is
published). See `CITATION.cff`.
