# A tight broad-line width–luminosity relation in little red dots

Data, analysis code, and manuscript sources for

> Haoyu Wang & Xinmiao Zhao, *A Tight Broad-line Width–Luminosity Relation in Little Red Dots and Its Consequence for
> Virial Black Hole Masses* (submitted).

We uniformly refit public JWST/NIRSpec grating spectra of 29 little red dots (LRDs) and 15 other broad-line sources,
calibrate unmodeled per-spectrum variance from repeat exposures, measure the intrinsic scatter of broad Hα widths at
fixed line luminosity across six width definitions, compare with published JWST measurements and identically refitted
SDSS DR7 quasars, and test the widths against both bound orbital motion and electron-scattering envelopes around an
Eddington-limited photosphere.

## Layout

| Path | Contents |
|---|---|
| `lrdwidth/` | Python package: grid posteriors (`stats.py`), regressions (`regression.py`, `photosphere.py`), emission-line fitting engine (`linefit.py`, `fit_batch.py`) |
| `analysis/` | Scripts turning the per-source measurements in `data/` into the derived results in `data/results/` |
| `scripts/` | Scripts turning `data/results/` into every number (`paper/numbers.tex`), table (`paper/tables/`), and figure (`paper/figures/`) of the manuscript |
| `data/` | Per-source measurements, repeat-spectrum pairs, sample records, and stored results; see [`data/README.md`](data/README.md) |
| `paper/` | AASTeX v7.0.1 source (`ms.tex`), generated macros/tables/figures, and compiled PDF (`ms.pdf`, `ms_arxiv.pdf`) |
| `third_party/` | Licenses of bundled third-party files and URLs/hashes of external spectra not redistributed |

## Installation

Python 3.10 or later:

```sh
pip install -e .
```

This installs `numpy`, `scipy`, `matplotlib`, and `astropy`. Exact package versions used for the paper are pinned in
`requirements.txt`.

## Reproducing the paper

Every numerical result, table, and figure in the manuscript is generated automatically from the stored results:

```sh
python scripts/make_tables.py     # paper/tables/ (Tables 1-3 and Appendix Table B1)
python scripts/make_numbers.py    # paper/numbers.tex (all 366 macros used in ms.tex and tables)
python scripts/make_figures.py    # paper/figures/ (Figures 1-4)
cd paper && pdflatex ms && bibtex ms && pdflatex ms && pdflatex ms
```

or simply `make paper`. No numbers are hand-typed in `paper/ms.tex`.

The stored results in `data/results/` are regenerated from the per-source measurements via `make analysis`
(which runs all fast scripts; `08_photosphere_simulations.py` runs 2000 simulations per sample separately):

| Script | Result file in `data/results/` | Contents |
|---|---|---|
| `01_width_luminosity.py` | `width_luminosity.json` | regressions across all six measures of width and photospheric comparisons |
| `02_robustness.py` | `robustness.json` | sample variants, halved LSF width, BX26 comparisons, virial-mass residuals, disk orientation scatters |
| `03_joint_line_photosphere.py` | `joint_line_photosphere.json` | joint regression on line and photospheric luminosity |
| `04_joint_robustness.py` | `joint_robustness.json` | robustness of the joint regression |
| `05_photosphere_constrained.py` | `photosphere_constrained.json` | photospheric regressions on sources with well-constrained photospheric luminosity |
| `06_membership.py` | `membership.json` | sensitivity to sample quality and label rules |
| `07_snr_threshold.py` | `snr_threshold.json` | sensitivity to the broad-line signal-to-noise threshold |
| `08_photosphere_simulations.py` | `photosphere_simulations.json` | partial regressions and hierarchical/flat simulations under line vs. photosphere hypotheses |
| `09_width_definition_ratios.py` | `width_definition_ratios.json` | ratio to the photospheric Keplerian speed across width definitions |
| `10_scattering_kernel.py` | `scattering_kernel.json` | broadening and optical-depth sensitivity for Gaussian and top-hat cores convolved with an exponential scattering kernel |
| `11_repeat_calibration.py` | `repeat_calibration.json` | repeat-spectrum error floors from 9 multi-spectrum LRDs (`data/widths/repeat_pairs_22.tsv`), one-sided 95% upper limits, high-S/N scatter, BX26 sky deduplication, and SDSS width offset |
| `12_ancillary_checks.py` | `ancillary_checks.json` | BIC-preferred width regression, MAD/leave-one-out diagnostics, residual-vs-S/N check, and wide-window summary |

Run `python scripts/check_reproduction.py` for a fast end-to-end verification of posteriors, hierarchical likelihoods,
and SHA-256 reference hashes of `data/results/numbers.json` and `paper/numbers.tex`.

## What is and is not included

Included: all per-source broad Hα measurements (six width definitions with full covariances), repeat-spectrum pair
measurements, sample selection audits, injection–recovery and hierarchical simulation code and summaries, and the full
manuscript build pipeline.

Not redistributed (needed only to re-fit raw detector spectra from scratch; see
[`analysis/REPRODUCING_MEASUREMENTS.md`](analysis/REPRODUCING_MEASUREMENTS.md)):

- **Raw 1D FITS spectra** from the public DAWN JWST Archive and SDSS DR17 (`third_party/external_data.tsv` lists every URL
  and SHA-256 checksum; `scripts/fetch_calibration.sh` downloads the NIRSpec dispersion tables).
- **Third-party catalogues without a redistribution license**, specifically Shen et al. (2011) and Kocevski et al. (2026).

## Use of AI tools

Portions of the analysis pipeline and manuscript text were developed and edited with the assistance of large language
models (Gemini, Claude, and ChatGPT/Codex), as disclosed in the manuscript. All scientific choices, calculations, and
conclusions were verified by the authors, who take full responsibility for the work.

## License and citation

Code: MIT (`LICENSE`). Measurements, tables, figures, and documentation: CC BY 4.0 (`LICENSE-DATA.md`).
Third-party material retains its original terms; see `ATTRIBUTION.md` and `third_party/`.

Archive DOI: `10.5281/zenodo.24182569` (see `CITATION.cff` and `.zenodo.json`).
