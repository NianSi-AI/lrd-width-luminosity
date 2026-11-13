"""Ratio of the line width to the photospheric Keplerian speed, and its scatter, for each measure of width.

Writes data/results/width_definition_ratios.json.
"""

import json
import math

import numpy as np

from lrdwidth.paths import PHOTOSPHERE, RESULTS
from lrdwidth.photosphere import one_dim

J = json.load(open(PHOTOSPHERE / "samples.json"))
out = {}
for key in ("main_E", "Curtis_B_own", "Curtis_A_own"):
    rows = J[key]
    out[key] = {}
    for m in (
        "single_gaussian",
        "nonparametric_FWHM",
        "double_gaussian",
        "lorentzian",
        "exponential",
        "nonparametric_W50",
    ):
        q = []
        v = []
        for r in rows:
            w = r.get("widths", {}).get(m)
            if not w or not w.get("width_kms") or (not np.isfinite(w["width_kms"])):
                continue
            lvk = r["y"] - r["q_dex"]
            q.append(math.log10(w["width_kms"]) - lvk)
            v.append(w["cov"][1][1] + (0.25 * r["sx"]) ** 2 + r["st"] ** 2)
        q = np.array(q)
        v = np.array(v)
        if len(q) < 3:
            continue
        sig, _ = one_dim(q, v)
        out[key][m] = dict(
            N=len(q),
            ratio_median=float(10 ** np.median(q)),
            ratio_p16_84=[float(10 ** np.quantile(q, 0.16)), float(10 ** np.quantile(q, 0.84))],
            sigma=sig,
        )
    r0 = rows[0]
json.dump(out, open(RESULTS / "width_definition_ratios.json", "w"), indent=1)
