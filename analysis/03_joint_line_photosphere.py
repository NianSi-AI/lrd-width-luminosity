"""Width regressed jointly on line luminosity and photospheric luminosity for the same sources:
y = a + b_line x_line + b_ph x_ph.

Writes data/results/joint_line_photosphere.json.
"""

import csv
import json
import math

import numpy as np

from lrdwidth.paths import DATA, PHOTOSPHERE, RESULTS
from lrdwidth.regression import hpd_level, load_line_inputs, post3g, quant

W1, byf, Q, byx = load_line_inputs(DATA)
out = {}
for tag, fn in (("A", "atmosphere_bx26_widths.tsv"), ("B", "atmosphere_our_widths.tsv")):
    rows = list(csv.DictReader(open(PHOTOSPHERE / fn), delimiter="\t"))
    xph = np.array([float(r["x"]) for r in rows])
    t = np.array([float(r["t"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    sph = np.array([float(r["sx"]) for r in rows])
    st = np.array([float(r["st"]) for r in rows])
    sy = np.array([float(r["sy"]) for r in rows])
    if tag == "A":
        xl = np.array([byx[r["width_id"]]["logLbol"] - 45 for r in rows])
        sl = np.array([byx[r["width_id"]]["sx"] for r in rows])
    else:
        xl = np.array([byf[r["width_id"]]["x"] for r in rows])
        sl = np.array([math.sqrt(byf[r["width_id"]]["cov"][0][0]) for r in rows])
    b1, b2, s, p = post3g(xl, xph, y, sl, sph, sy)
    pbc = p.sum(2)
    r = dict(
        N=len(y),
        b_line=quant(b1, p.sum((1, 2))),
        b_ph=quant(b2, p.sum((0, 2))),
        sigma=quant(s, p.sum((0, 1))),
        P_bph_gt0=float(p.sum((0, 2))[b2 > 0].sum()),
        P_bline_gt0=float(p.sum((1, 2))[b1 > 0].sum()),
        hpd_line0_ph025=hpd_level(b1, b2, pbc, 0.0, 0.25),
        hpd_ph0=None,
    )
    # contour through the maximum of the slice b_ph = 0, and through that of b_ph = 0.25
    j0 = int(np.argmin(abs(b2)))
    i0 = int(np.argmax(pbc[:, j0]))
    r["hpd_ph0"] = float(pbc[pbc > pbc[i0, j0]].sum())
    j25 = int(np.argmin(abs(b2 - 0.25)))
    i25 = int(np.argmax(pbc[:, j25]))
    r["hpd_ph025_best_line"] = float(pbc[pbc > pbc[i25, j25]].sum())
    r["bline_at_ph025"] = float(b1[i25])
    out[tag] = r
    c1 = np.polyfit(xl, y, 1)
    ry = y - np.polyval(c1, xl)
    c2 = np.polyfit(xl, xph, 1)
    rph = xph - np.polyval(c2, xl)
json.dump(out, open(RESULTS / "joint_line_photosphere.json", "w"), indent=1)
