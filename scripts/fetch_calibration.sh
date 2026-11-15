#!/usr/bin/env bash
# Download the NIRSpec dispersion tables distributed with msaexp (https://github.com/gbrammer/msaexp) into data/calibration/.
# They are needed only by lrdwidth.linefit, i.e. to refit spectra; the numbers, tables and figures do not use them.
#
# The tables are pinned to msaexp commit 9d342cd8 ("fix units of dispersion tables", 2023-02-17), the last commit that
# changed them. Its three files reproduce the SHA-256 checksums recorded for the measurements in
# third_party/external_data.tsv; every download is verified against them and the script stops on any mismatch.
set -euo pipefail
MSAEXP_COMMIT=9d342cd85f2dae2435f14c5201733a5e98d8d9a5
dir="$(cd "$(dirname "$0")/.." && pwd)/data/calibration"
mkdir -p "$dir"

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

while read -r g sha; do
    out="$dir/jwst_nirspec_${g}_disp.fits"
    curl -fsSL -o "$out.part" \
        "https://raw.githubusercontent.com/gbrammer/msaexp/${MSAEXP_COMMIT}/msaexp/data/jwst_nirspec_${g}_disp.fits"
    got="$(sha256_of "$out.part")"
    if [ "$got" != "$sha" ]; then
        rm -f "$out.part"
        echo "checksum mismatch for jwst_nirspec_${g}_disp.fits: expected $sha, got $got" >&2
        exit 1
    fi
    mv "$out.part" "$out"
done <<'EOF'
g395m 2160a0c7b6630028d4ca1766c570b1aaa5bfd26c3b42060d81bfbcb874be0db4
g395h 46df8738580a1386b705801add7e6433ad678e1b23476762a2f37e398fe5434c
g235m 376505b675280ee852fbd045ca035eda6d9124bbb6c526ef30396dbab00f5fc9
EOF
echo "dispersion tables (msaexp ${MSAEXP_COMMIT:0:8}, SHA-256 verified) written to $dir"
