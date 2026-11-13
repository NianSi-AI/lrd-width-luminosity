#!/usr/bin/env bash
# Download the NIRSpec dispersion tables distributed with msaexp (https://github.com/gbrammer/msaexp) into data/calibration/.
# They are needed only by lrdwidth.linefit, i.e. to refit spectra; the numbers, tables and figures do not use them.
set -euo pipefail
dir="$(cd "$(dirname "$0")/.." && pwd)/data/calibration"
mkdir -p "$dir"
for g in g395m g395h g235m; do
    curl -fsSL -o "$dir/jwst_nirspec_${g}_disp.fits" \
        "https://raw.githubusercontent.com/gbrammer/msaexp/main/msaexp/data/jwst_nirspec_${g}_disp.fits"
done
echo "dispersion tables written to $dir"
