"""Locations of inputs and outputs, relative to the repository root.

Every location can be overridden by an environment variable, so the code also runs
against data kept elsewhere.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def _p(var, default):
    return Path(os.environ.get(var, str(ROOT / default)))


SAMPLES_DIR = _p('SYNTH_SAMPLES_DIR', 'samples')               # samples/<formula>/buttons/<NNNN>/: scans
SEARCH_DIR = _p('SYNTH_SEARCH_DIR', 'databases/searchmatch')   # inputs of the search-match library
LIBRARY_DIR = _p('SYNTH_LIBRARY_DIR', 'results/searchmatch/library')   # generated library cache
INSTPRM_DIR = _p('SYNTH_INSTPRM_DIR', 'instruments')            # instrument profiles
# Optional ICDD PDF-card CIFs supplied by the user (not distributed with this repository):
# any '*.cif' below this folder whose elements lie in the chemical space is added to the library.
PDF_CARD_DIR = Path(os.environ['SYNTH_PDF_CARD_DIR']) if os.environ.get('SYNTH_PDF_CARD_DIR') else None


def button_dir(sample_number):
    """Folder of one button, samples/<formula>/buttons/<NNNN>/, or None."""
    hits = sorted(SAMPLES_DIR.glob(f'*/buttons/{int(sample_number):04d}'))
    return hits[0] if hits else None
