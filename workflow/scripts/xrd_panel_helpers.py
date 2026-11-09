"""Helpers shared by the XRD figure scripts (from the lab repository's plot_xrd_formed.py)."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PHASEID = Path('.')

COLOR_FIT = '#D32F2F'          # red fit line over black data crosses (Rietveld convention)
COLOR_TICK_OTHER = '0.45'
REPORTED_MARK = {'yes': '*', 'composition': r'$\dagger$', 'series': r'$\ddagger$', 'no': ''}


def parse_phases(text):
    """'key=display;key=display' -> list of (key, display)."""
    out = []
    for item in text.split(';'):
        key, display = item.split('=', 1)
        out.append((key.strip(), display.strip()))
    return out



def load_fit(path):
    """GSAS-II export: two_theta obs calc background diff, plus the tick file."""
    x, yo, yc, bg, diff = np.loadtxt(path, unpack=True)
    ticks = json.loads(path.with_name(path.stem + '_ticks.json').read_text())
    return x, yo, yc, bg, diff, {k: np.asarray(v, float) for k, v in ticks.items()}



def lookup_fit_stats(row):
    """(wR, {key: weight fraction}) from the table that produced the fit."""
    s = int(row['fit_sample'])
    table, key = row['source_table'], str(row['source_key'])
    if table == 'refinements':
        df = pd.read_csv(PHASEID / f'{s:04d}' / 'refinements.csv')
        hit = df[(df['rank'] == int(key)) & (df['mode'] == 'free') & (df['failed'] != True)]  # noqa: E712
        if hit.empty:
            raise KeyError(f'{s:04d}: no converged free row with rank {key} in refinements.csv')
        rec = hit.iloc[0]
        keys = _phase_keys(rec['phases'])
        wts = [float(w) for w in str(rec['wt']).split(' / ')]
        return float(rec['wR']), dict(zip(keys, wts))
    if table == 'verdict_extra':
        df = pd.read_csv(PHASEID / 'verdict_extra.csv')
        hit = df[(df['sample'] == s) & (df['label'] == key) & (df['failed'] != True)]  # noqa: E712
        if hit.empty:
            raise KeyError(f'{s:04d}: no converged row labelled {key!r} in verdict_extra.csv')
        rec = hit.iloc[0]
        keys = _phase_keys(rec['phases'])
        wts = [float(w) for w in str(rec['wt']).split(' / ')]
        return float(rec['wR']), dict(zip(keys, wts))
    if table == 'summary_json':
        rows = json.loads((PHASEID / f'{s:04d}' / 'verdict' / 'summary.json').read_text())
        hit = [r for r in rows if r['label'] == key and not r.get('failed')]
        if not hit:
            raise KeyError(f'{s:04d}: no converged entry labelled {key!r} in verdict/summary.json')
        rec = hit[0]
        # The json carries no phase names; they are the keys of the tick file,
        # written in the same phase order.
        ticks = json.loads(Path(str(row['_fit_path'])).with_name(
            Path(str(row['_fit_path'])).stem + '_ticks.json').read_text())
        return float(rec['wR']), dict(zip(ticks.keys(), [float(w) for w in rec['wt']]))
    raise ValueError(f'unknown source_table {table!r}')


