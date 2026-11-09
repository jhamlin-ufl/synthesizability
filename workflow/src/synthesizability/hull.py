"""Hull energy Delta H_f - C(x) of a target in one database, from databases/snapshot/.

For a target and a database: collect every entry of the ternary system and its
binary and elemental subsystems (energy = formation energy per atom x atoms),
take the lowest entry at the target composition, remove every entry at that
composition, build the convex hull of the rest with pymatgen and return
E_f(target) - E_hull(target composition).  Negative = the database predicts the
target stable, by that margin; positive = how far above the hull it lies.
This is OQMD's "stability" construction, applied uniformly to all databases
(MP and Alexandria report only the energy above the hull, zero when stable).
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition

SNAPSHOT = Path('databases/snapshot')
DATABASES = {   # column prefix: snapshot folder
    'oqmd': SNAPSHOT / 'oqmd_v1.7',
    'mp': SNAPSHOT / 'mp',
    'alex_pbe': SNAPSHOT / 'alexandria_pbe_2025.07.02',
    'alex_pbesol': SNAPSHOT / 'alexandria_pbesol_optimade',
}


def system_entries(db_dir, elements):
    """PDEntry objects for every entry of the system and its subsystems, plus a
    zero-energy reference for any element the database file lacks."""
    entries = []
    for n in (1, 2, 3):
        for combo in combinations(sorted(elements), n):
            f = Path(db_dir) / ('-'.join(combo) + '.json')
            if not f.exists():
                continue
            for e in json.loads(f.read_text()).get('entries', []):
                comp = Composition(e['composition_id'])
                entries.append(PDEntry(comp, e['delta_e'] * comp.num_atoms,
                                       name=str(e.get('entry_id', e.get('mp_id')))))
    have = {e.composition.reduced_formula for e in entries if len(e.composition) == 1}
    entries += [PDEntry(Composition(el), 0.0, name=f'{el}-reference') for el in elements if el not in have]
    return entries


def hull_energy_of_rest(entries, formula):
    """Hull energy (eV/atom) at `formula` from all entries at other compositions."""
    target = Composition(formula).reduced_composition
    others = [e for e in entries if e.composition.reduced_composition != target]
    return PhaseDiagram(others).get_hull_energy_per_atom(target)


def hull_distance(entries, formula):
    """(n entries at the composition, lowest formation energy, Delta H_f - C(x)),
    energies in eV/atom; (0, None, None) if the database has no entry there."""
    target = Composition(formula).reduced_composition
    same = [e for e in entries if e.composition.reduced_composition == target]
    if not same:
        return 0, None, None
    best = min(e.energy_per_atom for e in same)
    return len(same), best, best - hull_energy_of_rest(entries, formula)
