"""The 34 targets with everything the tables and the figure print about them.

Joins, on the sample number:
  tables/paper_samples.csv  list (A/B/C), category, alt_category, confidence, reported,
                                  table_row (see tables/README.md)
  tables/outcomes.csv       one_or_three ('yes': category 3 where the ordered structure,
                                  category 1, cannot be excluded; drawn as a 3/1 split marker),
                                  outcome (Table 1; ' | ' forces the line break) and
                                  phases (Supplementary phases table)
  databases/snapshot/targets.csv       formula and the OQMD entry the target was selected from
  computed/p_disorder.csv      P_disorder

Rows come out in table order (category, then decreasing P_disorder) with ``row``
(the printed number, checked against table_row), ``alt`` (alternative category
or <NA>) and ``split``.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .formula import ordered

SAMPLES = Path('tables/paper_samples.csv')
OUTCOMES = Path('tables/outcomes.csv')
TARGETS = Path('databases/snapshot/targets.csv')
P_DISORDER = Path('computed/p_disorder.csv')

CATEGORY_NAME = {
    1: 'On-hull structure',
    2: 'Above-hull polymorph',
    3: 'Disordered parent of on-hull structure',
    4: 'Other solid solution',
    5: 'Multiphase, on-hull structure not formed',
    6: 'On-hull structure did not form but structure(s) unresolved',
}


def load_targets() -> pd.DataFrame:
    s = pd.read_csv(SAMPLES, dtype={'alt_category': 'Int64'})
    s = s.rename(columns={'sample': 'sample_number', 'formula': 'printed_formula'})
    o = pd.read_csv(OUTCOMES).rename(columns={'sample': 'sample_number'})
    t = pd.read_csv(TARGETS)[['sample_number', 'formula', 'oqmd_entry_id']]
    p = pd.read_csv(P_DISORDER)[['sample_number', 'disorder_probability']]
    df = t.merge(s, on='sample_number', how='outer', validate='1:1')
    df = df.merge(o, on='sample_number', how='outer', validate='1:1').merge(p, on='sample_number', how='left')
    if len(df) != 34 or df[['formula', 'category', 'outcome', 'disorder_probability']].isna().any().any():
        raise ValueError('the sample tables, the snapshot targets and P_disorder do not cover the same 34 targets')
    bad = df[df['formula'].map(ordered) != df['printed_formula']]
    if not bad.empty:
        raise ValueError(f'formula mismatch: {bad[["sample_number", "formula", "printed_formula"]].values.tolist()}')
    one_or_three = df['one_or_three'].eq('yes')
    if (one_or_three & df['category'].ne(3)).any():
        raise ValueError('one_or_three rows must be category 3')
    df['alt'] = df['alt_category'].mask(one_or_three, 1)
    df['split'] = df['alt'].notna()
    df = df.sort_values(['category', 'disorder_probability', 'sample_number'],
                        ascending=[True, False, True]).reset_index(drop=True)
    df['row'] = range(1, len(df) + 1)
    if not df['row'].equals(df['table_row']):
        raise ValueError('row order differs from table_row in paper_samples.csv')
    return df
