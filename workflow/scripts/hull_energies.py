#!/usr/bin/env python3
"""Delta H_f - C(x) of the 34 targets in OQMD, MP, Alexandria PBE and PBEsol.

  in:  databases/snapshot/  (targets.csv and the four database extracts)
  out: results/hull_energies.csv

Construction in src/synthesizability/hull.py.  Columns per database
(oqmd, mp, alex_pbe, alex_pbesol): <db>_n_entries (entries at the target
composition), <db>_formation_meV (lowest formation energy there) and
<db>_hull_meV (Delta H_f - C(x), meV/atom; blank = no entry).
oqmd_stored_meV is the stability OQMD itself stores for the selected entry,
kept only for comparison: the paper uses the recomputed value everywhere.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, 'workflow/src')
import pandas as pd
from pymatgen.core import Composition

from synthesizability.hull import DATABASES, hull_distance, system_entries

OUT = Path('results/hull_energies.csv')


def stored_oqmd_stability(system, entry_id):
    for e in json.loads((DATABASES['oqmd'] / f'{system}.json').read_text())['entries']:
        if e['entry_id'] == entry_id:
            return round(1000 * e['stability'], 1)
    return None


def main():
    targets = pd.read_csv('databases/snapshot/targets.csv')
    rows = []
    for _, t in targets.iterrows():
        elements = [el.symbol for el in Composition(t.formula).elements]
        row = {'sample_number': t.sample_number, 'formula': t.formula}
        for db, folder in DATABASES.items():
            n, ef, dist = hull_distance(system_entries(folder, elements), t.formula)
            row[f'{db}_n_entries'] = n
            row[f'{db}_formation_meV'] = None if ef is None else round(1000 * ef, 1)
            row[f'{db}_hull_meV'] = None if dist is None else round(1000 * dist, 1)
        row['oqmd_stored_meV'] = stored_oqmd_stability(t.system, int(t.oqmd_entry_id))
        rows.append(row)
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    diff = (out['oqmd_hull_meV'] - out['oqmd_stored_meV']).abs()
    print(f'{len(out)} targets -> {OUT}')
    print(f'OQMD recomputed vs stored: {int((diff <= 0.1).sum())} of {diff.notna().sum()} agree to 0.1 meV/atom; '
          + ', '.join(f'{f} {a:+.1f} vs {b:+.1f}' for f, a, b, d in
                      zip(out.formula, out.oqmd_hull_meV, out.oqmd_stored_meV, diff) if d > 0.1))


if __name__ == '__main__':
    main()
