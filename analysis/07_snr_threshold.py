"""Intrinsic scatter as the signal-to-noise threshold on the broad component is raised.

Writes data/results/snr_threshold.json.
"""

import json

import numpy as np

from lrdwidth import stats as S
from lrdwidth.paths import RESULTS, WIDTHS

W1 = json.load(open(WIDTHS / "regression_inputs.json"))["single_gaussian"]
out = {}
for cls in ("LRD", "non-LRD"):
    rows = [t for t in W1 if (t["cls"] == "LRD") == (cls == "LRD")]
    for thr in (5, 8, 10, 15, 20):
        r = [t for t in rows if t["broad_SNR"] >= thr]
        if len(r) < 6:
            continue
        x = np.array([t["x"] for t in r])
        y = np.array([t["y"] for t in r])
        C = np.array([t["cov"] for t in r])
        b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
        qs = S.quant(s, p.sum(0))
        qb = S.quant(b, p.sum(1))
        out[f"{cls}_SNge{thr}"] = dict(N=len(r), sigma=qs, slope=qb)
json.dump(out, open(RESULTS / "snr_threshold.json", "w"), indent=1)
