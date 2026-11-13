"""Width-luminosity scatter of SDSS DR7 quasars from the catalogue widths of Shen et al. (2011).

Usage: python analysis/sdss/catalog_benchmark.py <table.tsv>
The table has one row per quasar with columns x (log L_Halpha - 43), y (log FWHM in km/s) and cov (JSON 2x2
covariance). It is built from the Shen et al. catalogue, which is not redistributed here (see third_party/).

Writes data/results/sdss_catalog_benchmark.json.
"""

import csv
import json
import sys

import numpy as np

from lrdwidth.paths import RESULTS
from lrdwidth.regression import ols_sd, post2, quant

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
x = np.array([float(r["x"]) for r in rows])
y = np.array([float(r["y"]) for r in rows])
C = np.array([json.loads(r["cov"]) for r in rows])
b, s, p = post2(
    x,
    y,
    C[:, 0, 0],
    C[:, 1, 1],
    C[:, 0, 1],
    b_lo=-0.05,
    b_hi=0.2,
    nb=250,
    s_lo=0.18,
    s_hi=0.24,
    ns=300,
)
edge = float(max(p.sum(0)[0], p.sum(0)[-1], p.sum(1)[0], p.sum(1)[-1]))
sd, sl = ols_sd(x, y)
out = dict(
    N=len(rows),
    sigma=quant(s, p.sum(0)),
    slope=quant(b, p.sum(1)),
    ols_sd=sd,
    ols_slope=sl,
    edge_mass=edge,
    x_range=[float(x.min()), float(x.max())],
    fwhm_min=float(10 ** y.min()),
)
json.dump(out, open(RESULTS / "sdss_catalog_benchmark.json", "w"), indent=1)
