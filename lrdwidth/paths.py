"""Locations of the data, the stored results, and the manuscript."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = DATA / "results"
PAPER = ROOT / "paper"
WIDTHS = DATA / "widths"
PHOTOSPHERE = DATA / "photosphere"
