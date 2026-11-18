"""Partial regressions on line and photospheric luminosity, with simulations of widths that follow the photospheric
scaling or the line luminosity alone.

Usage: python analysis/08_photosphere_simulations.py [number of simulations, default 2000]
Writes data/results/photosphere_simulations.json.
"""

import json
import sys

from lrdwidth.paths import PHOTOSPHERE, RESULTS
from lrdwidth.photosphere import analyse, arr, load_main, load_sample

nsim = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
base = json.load(open(PHOTOSPHERE / "simulation_baselines.json"))
main = load_main(PHOTOSPHERE / "main_sample_sources.json")
stored = load_sample(PHOTOSPHERE / "samples.json", "main_E")
keys = ("xl", "cll", "cly", "cyy", "y", "xp", "sx", "t", "st", "logvk")
res = {"input_check_main_E": {k: float(max(abs(arr(main, k) - arr(stored, k)))) for k in keys}}
res["main_E"] = analyse(main, nsim, 1, base["main_E"])
for i, key in enumerate(("Curtis_B_own", "Curtis_A_own")):
    res[key] = analyse(load_sample(PHOTOSPHERE / "samples.json", key), nsim, 2 + i, base[key])
json.dump(res, open(RESULTS / "photosphere_simulations.json", "w"), indent=1)
