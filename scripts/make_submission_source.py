"""Build the comment-free source tree for arXiv and the journal.
arXiv publishes the TeX source, so every comment in ms.tex, numbers.tex and the table files is removed.
Usage: python scripts/make_submission_source.py [out dir]   (default: paper/submission_source). Run after the final
pdflatex/bibtex pass, which produces paper/ms.bbl.
Full-line comments are deleted; a trailing comment is deleted together with the space before it. \\% is kept.
"""

import os
import re
import shutil
import sys

from lrdwidth.paths import PAPER

H = str(PAPER)
O = os.path.join(H, sys.argv[1] if len(sys.argv) > 1 else "submission_source")


def strip(text):
    out = []
    for ln in text.split("\n"):
        if ln.lstrip().startswith("%"):
            continue
        m = re.search(r"(?<!\\)%", ln)
        if m:
            assert ln[: m.start()].strip(), ln  # trailing comment after real content
            ln = ln[: m.start()].rstrip()
        out.append(ln)
    s = re.sub(r"\n{3,}", "\n\n", "\n".join(out))
    assert not re.search(r"(?<!\\)%", s) and not re.search(
        r"[^\x00-\x7f]", s
    ), "comment or non-ASCII character left"
    return s


for d in ("", "figures", "tables"):
    os.makedirs(os.path.join(O, d), exist_ok=True)
for f in ("ms.tex", "numbers.tex", "tables/Table1_stub.tex", "tables/Table2_fits.tex"):
    open(os.path.join(O, f), "w", encoding="ascii").write(
        strip(open(os.path.join(H, f), encoding="utf-8").read())
    )
for f in [
    "ms.bbl",
    "refs.bib",
    "refs_ads.bib",
    "aastex701.cls",
    "aasjournalv7.bst",
    "orcid-ID.png",
] + ["figures/Figure%d.pdf" % i for i in range(1, 6)]:
    shutil.copyfile(os.path.join(H, f), os.path.join(O, f))
left = [
    l
    for l in open(os.path.join(O, "ms.tex"))
    if "TODO" in l or "example.com" in l or "Author Name" in l
]
print("written to", O)
print("placeholders still in ms.tex:", len(left))
[print("  ", l.strip()[:100]) for l in left]
