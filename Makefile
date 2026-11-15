PY ?= python

.PHONY: analysis numbers tables figures paper check-reproduction measurement-summary

analysis:
	for s in analysis/0[1-7]_*.py analysis/09_*.py analysis/11_*.py analysis/10_*.py analysis/12_*.py analysis/13_*.py; do $(PY) $$s || exit 1; done

numbers:
	$(PY) scripts/make_numbers.py

tables:
	$(PY) scripts/make_tables.py

figures:
	$(PY) scripts/make_figures.py

paper: tables numbers figures
	cd paper && pdflatex -interaction=nonstopmode ms && bibtex ms && pdflatex -interaction=nonstopmode ms && pdflatex -interaction=nonstopmode ms

check-reproduction:
	$(PY) scripts/check_reproduction.py

measurement-summary:
	$(PY) scripts/make_measurement_summary.py
