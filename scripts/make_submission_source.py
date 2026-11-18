"""Build the comment-free source trees for ApJ (named and dual-anonymous) and arXiv.
Usage: python scripts/make_submission_source.py [out dir]
Run after the final pdflatex/bibtex pass in paper/.
"""
import os
import re
import shutil
import sys

from lrdwidth.paths import PAPER

H = str(PAPER)
O = os.path.join(H, sys.argv[1] if len(sys.argv) > 1 else "submission_source")
O_ARX = os.path.join(H, "submission_source_arxiv")
O_ANON = os.path.join(H, "submission_source_apj_anonymous")


def strip(text):
    out = []
    for ln in text.split("\n"):
        if ln.lstrip().startswith("%"):
            continue
        m = re.search(r"(?<!\\)%", ln)
        if m:
            assert ln[: m.start()].strip(), ln
            ln = ln[: m.start()].rstrip()
        out.append(ln)
    s = re.sub(r"\n{3,}", "\n\n", "\n".join(out))
    assert not re.search(r"(?<!\\)%", s) and not re.search(r"[^\x00-\x7f]", s), "comment or non-ASCII character left"
    return s


TEX_FILES = (
    "ms.tex",
    "numbers.tex",
    "tables/Table1_stub.tex",
    "tables/Table2_fits.tex",
    "tables/Table3_readings.tex",
    "tables/TableB1_injections.tex",
)


def build_tree(target_dir, mode="apj_named"):
    for d in ("", "figures", "tables"):
        os.makedirs(os.path.join(target_dir, d), exist_ok=True)
    for f in TEX_FILES:
        txt = strip(open(os.path.join(H, f), encoding="utf-8").read())
        if f == "ms.tex":
            if mode == "arxiv":
                txt = txt.replace(r"\documentclass[twocolumn,linenumbers]{aastex701}", r"\documentclass[twocolumn]{aastex701}")
            elif mode == "apj_anonymous":
                txt = txt.replace(r"\documentclass[twocolumn,linenumbers]{aastex701}", r"\documentclass[twocolumn,linenumbers,anonymous]{aastex701}")
                txt = re.sub(r",pdfauthor=\{[^}]+\}", "", txt)
                txt, n_doi = re.subn(r"archived at\s*(?:\\dataset\[doi:[^\]]+\]\{[^}]+\}|\\url\{[^}]+\})", "archived in a public repository (URL withheld for review)", txt)
                assert n_doi == 1 and "zenodo" not in txt.lower() and "niansi" not in txt.lower(), "Repository URL/DOI left in the anonymous build"
        open(os.path.join(target_dir, f), "w", encoding="ascii").write(txt)
    for f in ("refs.bib", "refs_ads_20261113.bib"):
        t = "\n".join(l for l in open(os.path.join(H, f), encoding="utf-8").read().split("\n") if not l.lstrip().startswith("%"))
        open(os.path.join(target_dir, f), "w", encoding="utf-8").write(re.sub(r"\n{3,}", "\n\n", t).lstrip("\n"))
    copy_list = ["aastex701.cls", "aasjournalv7.bst", "orcid-ID.png", "tables/Table1_full.csv"] + [f"figures/Figure{i}.pdf" for i in range(1, 5)]
    for f in copy_list:
        shutil.copyfile(os.path.join(H, f), os.path.join(target_dir, f))
    bbl = open(os.path.join(H, "ms.bbl"), encoding="utf-8").read()
    bbl = re.sub(r"\\bibitem\[\{\s+\{", r"\\bibitem[{{", bbl)
    open(os.path.join(target_dir, "ms.bbl"), "w", encoding="utf-8").write(bbl)
    bad = [f for f in os.listdir(target_dir) if f.endswith(".bib") and re.search(r"%", open(os.path.join(target_dir, f), encoding="utf-8").read())]
    assert not bad, bad


build_tree(O, mode="apj_named")
build_tree(O_ARX, mode="arxiv")
build_tree(O_ANON, mode="apj_anonymous")
print("written to", O, ",", O_ARX, ", and", O_ANON)
