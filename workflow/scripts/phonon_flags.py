#!/usr/bin/env python3
"""Dynamical-stability cross of Fig. 2, from the frozen MatterSim phonon verdicts.

  in:  computed/phonons/candidates.csv         every OQMD, MP and Alexandria PBE entry at a
                                                  target composition within 50 meV/atom of the
                                                  hull of the other compositions, grouped into
                                                  distinct structures
       computed/phonons/phonon_results_5M.csv  verdict per structure (MatterSim-v1.0.0-5M +
                                                  phonopy, unstable below -0.1 THz)
       databases/snapshot/                             to re-check each entry's hull energy
  out: results/phonon_flags.csv

A target gets the cross when it lies below the surrounding hull in at least one
database and none of its below-hull structures is dynamically stable.  The
hull energies of the candidates were computed when the phonon run was set up;
they are recomputed here from the snapshot and must agree.
"""
import sys
from pathlib import Path

sys.path.insert(0, 'workflow/src')
import pandas as pd
from pymatgen.core import Composition

from synthesizability.hull import DATABASES, hull_energy_of_rest, system_entries

FROZEN = Path('computed/phonons')
OUT = Path('results/phonon_flags.csv')


def check_hull_energies(cands):
    hull = {}
    for (db, formula), g in cands.groupby(['db', 'formula']):
        els = [e.symbol for e in Composition(formula).elements]
        hull[(db, formula)] = 1000 * hull_energy_of_rest(system_entries(DATABASES[db], els), formula)
    recomputed = cands['formation_meV'] - [hull[k] for k in zip(cands['db'], cands['formula'])]
    worst = (recomputed - cands['hull_rest_meV']).abs().max()
    if worst > 0.11:
        raise ValueError(f'candidate hull energies differ from the snapshot by up to {worst:.2f} meV/atom')
    print(f'candidate hull energies reproduced from the snapshot (max difference {worst:.2f} meV/atom)')


def main():
    cands = pd.read_csv(FROZEN / 'candidates.csv')
    check_hull_energies(cands)
    # results are keyed by the representative entry; group ids come from candidates.csv
    res = pd.read_csv(FROZEN / 'phonon_results_5M.csv').drop_duplicates('key', keep='last')
    res = res.drop(columns='group').merge(cands[['key', 'group']], on='key')
    ent = cands.merge(res[['group', 'band_min_THz', 'status_paper']], on='group', how='left')
    if ent['status_paper'].isna().any():
        raise ValueError(f'no verdict for {ent.loc[ent.status_paper.isna(), "key"].tolist()}')
    rows = []
    for sample, g in ent.groupby('sample'):
        below = g[g['hull_rest_meV'] <= 0].drop_duplicates('group').sort_values('hull_rest_meV')
        n_stable = int(below['status_paper'].eq('Stable').sum())
        rows.append(dict(sample_number=sample, formula=g['formula'].iloc[0],
                         n_below_hull=len(below), n_stable_below_hull=n_stable,
                         below_hull_structures='; '.join(
                             f'{r.db} sg{r.spg} {r.hull_rest_meV:+.0f} {r.status_paper}' for r in below.itertuples()),
                         all_below_hull_unstable=bool(len(below) and n_stable == 0)))
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f'{len(out)} targets -> {OUT}; crossed: ' + ', '.join(out.loc[out.all_below_hull_unstable, 'formula']))


if __name__ == '__main__':
    main()
