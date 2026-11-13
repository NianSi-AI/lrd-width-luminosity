PY ?= python

.PHONY: analysis numbers tables figures paper

analysis:
	for s in analysis/0[1-7]_*.py analysis/09_*.py; do $(PY) $$s || exit 1; done

numbers:
	$(PY) scripts/make_numbers.py

tables:
	$(PY) scripts/make_tables.py

figures:
	$(PY) scripts/make_figures.py

paper: numbers tables figures
	cd paper && pdflatex ms && bibtex ms && pdflatex ms && pdflatex ms
