#!/usr/bin/env python3
"""Energy of each target relative to the convex hull in every cached database.

  data/processed/hull_distances.csv

OQMD reports, for a phase on the hull, how far *below* the hull of all other
compositions it sits (a negative "stability"), which is why its hull distances
can be negative.  MP and Alexandria report only the energy above the hull, zero
for stable phases.  To put the four databases on the same footing this script
recomputes the OQMD-style quantity for each of them from the cached entries:

  1. collect every cached entry of the ternary system and its binary and
     elemental subsystems, with energy = formation energy per atom x atoms;
  2. take the lowest-energy entry at the target composition;
  3. remove every entry at that composition, build the convex hull of the rest
     with pymatgen, and report E_f(target) - E_hull(composition).

Negative means the database predicts the target to be stable, positive by how
much it misses.  Checked against OQMD's own stability field: 32 of the 33
targets agree to 0.1 meV/atom, Hf2MoIr differs by 6.6 meV/atom (OQMD's internal
hull evidently differs from the cached subsystem entries there).

Columns per database (oqmd, mp, alex_pbe, alex_pbesol):
  <db>_n_entries        entries cached at the target composition (0 = none)
  <db>_formation_meV    lowest formation energy at the composition, meV/atom
  <db>_hull_meV         distance to the hull of the other compositions, meV/atom
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, 'src')

import pandas as pd
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition

FORMULAS_CSV = Path('data/processed/synthesis_data_no_disorder.csv')
OUT_CSV = Path('data/processed/hull_distances.csv')

DATABASES = {
    'oqmd': Path('data/external/oqmd_ternary_phases'),
    'mp': Path('data/external/mp_ternary_phases'),
    'alex_pbe': Path('data/external/alexandria_pbe_ternary_phases'),
    'alex_pbesol': Path('data/external/alexandria_pbesol_ternary_phases'),
}


def system_entries(db_dir, elements):
    """PDEntry objects for every cached entry of the system and its subsystems,
    plus zero-energy elemental references where the cache has none."""
    entries = []
    for n in (1, 2, 3):
        for combo in combinations(sorted(elements), n):
            f = db_dir / ('-'.join(combo) + '.json')
            if not f.exists():
                continue
            for e in json.load(open(f)).get('entries', []):
                comp = Composition(e['composition_id'])
                entries.append(PDEntry(comp, e['delta_e'] * comp.num_atoms,
                                       name=str(e.get('entry_id', e.get('mp_id')))))
    have = {e.composition.reduced_formula for e in entries if len(e.composition) == 1}
    for el in elements:
        if el not in have:
            entries.append(PDEntry(Composition(el), 0.0, name=f'{el}-reference'))
    return entries


def hull_distance(entries, formula):
    """(n_entries, formation energy, distance to the hull of other compositions),
    energies in eV/atom; (0, None, None) if the database has no entry."""
    target = Composition(formula).reduced_composition
    same = [e for e in entries if e.composition.reduced_composition == target]
    if not same:
        return 0, None, None
    best = min(same, key=lambda e: e.energy_per_atom)
    others = [e for e in entries if e.composition.reduced_composition != target]
    hull = PhaseDiagram(others).get_hull_energy_per_atom(target)
    return len(same), best.energy_per_atom, best.energy_per_atom - hull


def main():
    formulas = pd.read_csv(FORMULAS_CSV)['formula'].dropna().unique()
    rows = []
    for formula in formulas:
        elements = [el.symbol for el in Composition(formula).elements]
        if len(elements) != 3:
            continue
        row = {'formula': formula}
        for db, db_dir in DATABASES.items():
            n, ef, dist = hull_distance(system_entries(db_dir, elements), formula)
            row[f'{db}_n_entries'] = n
            row[f'{db}_formation_meV'] = None if ef is None else round(ef * 1000, 1)
            row[f'{db}_hull_meV'] = None if dist is None else round(dist * 1000, 1)
        rows.append(row)
        print(f"  {formula:<12} " + '  '.join(
            f"{db}={row[f'{db}_hull_meV']:+.1f}" if row[f'{db}_hull_meV'] is not None else f"{db}=--"
            for db in DATABASES))
    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f'Saved {OUT_CSV} ({len(out)} formulas)')


if __name__ == '__main__':
    main()
