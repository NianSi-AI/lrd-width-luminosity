# Third-party data

The CC BY 4.0 license of the author's measurements does not replace the terms or the attribution of upstream
products. Please cite the original sources.

- **de Graaff et al. (2026)**, *DJA/NIRSpec LRD sample and SED properties*, version 2,
  [doi:10.5281/zenodo.21977747](https://doi.org/10.5281/zenodo.21977747), CC BY 4.0. Used: temperatures, luminosities,
  redshifts and identifiers (`data/photosphere/blackbody_sources.tsv`, `main_sample_sources.json`,
  `data/sample/candidates.json`).
- **Curtis et al. (2026)**, [arXiv:2609.09265](https://arxiv.org/abs/2609.09265), data at
  [o-curtis/lrd-quasistar-fits](https://github.com/o-curtis/lrd-quasistar-fits), CC BY 4.0 (data). Used: photospheric
  parameters (`data/photosphere/atmosphere_*.tsv`, `samples.json`).
- **Baccus & Xu (2026)**, ApJ, 1006, 165, [doi:10.3847/1538-4357/ae7de7](https://doi.org/10.3847/1538-4357/ae7de7),
  CC BY 4.0. Used: the machine-readable Table 1 (`data/bx26/`, `data/photosphere/atmosphere_bx26_widths.tsv`).
- **Taylor et al. (2025)**, [arXiv:2409.06772](https://arxiv.org/abs/2409.06772), and **Hviding et al. (2025)**,
  [doi:10.5281/zenodo.16758549](https://doi.org/10.5281/zenodo.16758549), both CC BY 4.0. Used: source labels.
- **Kocevski et al. (2026)**: labels only; the catalogue has no explicit license and is not redistributed.
- **Shen et al. (2011)**, ApJS, 194, 45: the catalogue is not redistributed; only aggregate results and our own fits
  of the SDSS spectra are included.
- **DAWN JWST Archive (DJA)** and **msaexp**: the widths, fluxes and covariances here are new measurements made on
  the public DJA reductions. The spectra themselves are not redistributed.
- **SDSS**: the fits in `data/sdss/` were made on DR17 Legacy spectra, which are in the public domain.

`third_party/external_data.tsv` gives the address and checksum of every external file used. The AASTeX class,
bibliography style and ORCID icon in `paper/` keep their own terms; see `third_party/AASTeX_NOTICE.md`.
