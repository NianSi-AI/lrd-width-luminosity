"""Build the comment-free source trees for ApJ (with linenumbers) and arXiv (without linenumbers) for v11 (2026-11-14).
Usage: python3 make_submission_source.py [out dir]   (default: submission_source; also writes submission_source_arxiv).
Run after the final pdflatex/bibtex pass.
"""
import os, re, shutil, sys

from lrdwidth.paths import PAPER
H = str(PAPER)
O = os.path.join(H, sys.argv[1] if len(sys.argv) > 1 else 'submission_source')
O_ARX = os.path.join(H, 'submission_source_arxiv')

def strip(text):
    out = []
    for ln in text.split('\n'):
        if ln.lstrip().startswith('%'):
            continue
        m = re.search(r'(?<!\\)%', ln)
        if m:
            assert ln[:m.start()].strip(), ln
            ln = ln[:m.start()].rstrip()
        out.append(ln)
    s = re.sub(r'\n{3,}', '\n\n', '\n'.join(out))
    assert not re.search(r'(?<!\\)%', s) and not re.search(r'[^\x00-\x7f]', s), 'comment or non-ASCII character left'
    return s

TEX_FILES = (
    'ms.tex',
    'numbers.tex',
    'tables/Table1_stub.tex',
    'tables/Table2_fits.tex',
    'tables/Table3_readings.tex',
    'tables/TableB1_injections.tex',
)

def build_tree(target_dir, drop_linenumbers=False):
    for d in ('', 'figures', 'tables'):
        os.makedirs(os.path.join(target_dir, d), exist_ok=True)
    for f in TEX_FILES:
        txt = strip(open(os.path.join(H, f), encoding='utf-8').read())
        if drop_linenumbers and f == 'ms.tex':
            txt = txt.replace(r'\documentclass[twocolumn,linenumbers]{aastex701}', r'\documentclass[twocolumn]{aastex701}')
        open(os.path.join(target_dir, f), 'w', encoding='ascii').write(txt)
    for f in ('refs.bib', 'refs_ads_20261113.bib'):
        t = '\n'.join(l for l in open(os.path.join(H, f), encoding='utf-8').read().split('\n') if not l.lstrip().startswith('%'))
        open(os.path.join(target_dir, f), 'w', encoding='utf-8').write(re.sub(r'\n{3,}', '\n\n', t).lstrip('\n'))
    copy_list = ['ms.bbl', 'aastex701.cls', 'aasjournalv7.bst', 'orcid-ID.png', 'tables/Table1_full.csv'] + [f'figures/Figure{i}.pdf' for i in range(1, 5)]
    for f in copy_list:
        shutil.copyfile(os.path.join(H, f), os.path.join(target_dir, f))
    bad = [f for f in os.listdir(target_dir) if f.endswith('.bib') and re.search(r'%', open(os.path.join(target_dir, f), encoding='utf-8').read())]
    assert not bad, bad

build_tree(O, drop_linenumbers=False)
build_tree(O_ARX, drop_linenumbers=True)
print('written to', O, 'and', O_ARX)
