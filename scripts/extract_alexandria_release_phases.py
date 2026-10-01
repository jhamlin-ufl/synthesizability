#!/usr/bin/env python3
"""Build the Alexandria PBE phase caches from the bulk 2025.07.02 release files.

Replaces the OPTIMADE pull (query_alexandria_phases.py) for PBE.  The OPTIMADE
server serves only part of the database (it misses on-hull entries such as
TiPt5Si and TiPt5Ge that the 2025.07.02 release contains), whereas the release
files hold every entry with its structure, formation energy, energy above hull
and Alexandria's own "phase separation" energy (the energy relative to the hull
of all other compositions, the same construction as compute_hull_distances.py).

Inputs
  $ALEXANDRIA_DIR/pbe_2025.07.02/alexandria_*.json.bz2   (58 files, ~3.4 GB)
  data/processed/synthesis_data_no_disorder.csv           (formulas -> chemical spaces)

Outputs
  data/external/alexandria_pbe_ternary_phases/<space>.json   one file per unary,
      binary and ternary space of every target system, in the schema shared with
      the OQMD and MP caches (entry_id, composition_id, delta_e, stability, icsd)
      plus e_above_hull, e_phase_separation, spg, nsites, mat_id, source.
  data/external/alexandria_hull_structures/pbe_2025.07.02/<formula>_<mat_id>.cif
      every entry at a target composition (on or off the hull), and
  data/external/alexandria_hull_structures/index.csv  (one row per structure above)
  data/external/alexandria_pbe_ternary_phases/<space>/cifs/*.cif  the same
      structures in the layout the dashboard links to.

Each chunk is streamed and entries are pre-filtered on their "elements" list
with a regex before any JSON decoding, so the whole release takes minutes with
a few processes.  Run from the repository root.
"""
from __future__ import annotations

import bz2
import json
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

sys.path.insert(0, 'src')

import pandas as pd
from pymatgen.core import Composition, Structure
from pymatgen.io.cif import CifWriter

RELEASE = 'pbe_2025.07.02'
ALEXANDRIA_DIR = Path(os.environ.get('ALEXANDRIA_DIR', os.path.expanduser('~/data/alexandria'))) / RELEASE
FORMULAS_CSV = Path('data/processed/synthesis_data_no_disorder.csv')
OUT_DIR = Path('data/external/alexandria_pbe_ternary_phases')
CIF_DIR = Path('data/external/alexandria_hull_structures') / RELEASE
INDEX_CSV = Path('data/external/alexandria_hull_structures/index.csv')
HEADER = '{"@module": "pymatgen.entries.computed_entries", "@class": "ComputedStructureEntry"'
ELEMENTS_RE = re.compile(r'"elements": \[([^\]]*)\]')
N_PROC = int(os.environ.get('ALEXANDRIA_PROCS', '6'))


def target_spaces(formulas):
    """All unary/binary/ternary element sets spanned by the target systems."""
    spaces = set()
    targets = {}
    for f in formulas:
        comp = Composition(f)
        els = tuple(sorted(e.symbol for e in comp.elements))
        if len(els) != 3:
            continue
        targets[comp.reduced_formula] = f
        for n in (1, 2, 3):
            for sub in combinations(els, n):
                spaces.add(sub)
    return spaces, targets


def scan_file(path, spaces, targets):
    """Return (kept entries, structures at target compositions) from one chunk."""
    decoder = json.JSONDecoder()
    kept, structs = [], []
    buf = ''
    n = 0
    with bz2.open(path, 'rt') as fh:
        while True:
            chunk = fh.read(1 << 24)
            buf += chunk
            parts = buf.split(HEADER)
            buf = parts[-1] if chunk else ''
            for part in (parts[1:-1] if chunk else parts[1:]):
                n += 1
                m = ELEMENTS_RE.search(part)
                if not m:
                    continue
                els = tuple(sorted(s.strip().strip('"') for s in m.group(1).split(',') if s.strip()))
                if els not in spaces:
                    continue
                try:
                    obj, _ = decoder.raw_decode(HEADER + part)
                except json.JSONDecodeError:
                    continue
                d = obj['data']
                comp = Composition(obj['composition'])
                natoms = comp.num_atoms
                rec = {
                    'entry_id': d['mat_id'],
                    'composition_id': ' '.join(f'{el}{int(round(c))}' for el, c in sorted(comp.get_el_amt_dict().items())),
                    'delta_e': d.get('e_form'),
                    'stability': d.get('e_above_hull'),
                    'icsd': False,
                    'e_above_hull': d.get('e_above_hull'),
                    'e_phase_separation': d.get('e_phase_separation'),
                    'energy_per_atom': obj['energy'] / natoms,
                    'spg': d.get('spg'),
                    'nsites': d.get('nsites'),
                    'prototype_id': d.get('prototype_id'),
                    'source': RELEASE,
                    'space': '-'.join(els),
                }
                kept.append(rec)
                red = comp.reduced_formula
                if red in targets:
                    structs.append((targets[red], d['mat_id'], obj['structure'], rec))
            if not chunk:
                break
    return path.name, n, kept, structs


def main():
    formulas = pd.read_csv(FORMULAS_CSV)['formula'].dropna().unique()
    spaces, targets = target_spaces(formulas)
    files = sorted(ALEXANDRIA_DIR.glob('alexandria_*.json.bz2'))
    if not files:
        sys.exit(f'no release files in {ALEXANDRIA_DIR}')
    print(f'{len(files)} release files, {len(spaces)} chemical spaces, {len(targets)} target compositions, {N_PROC} processes')
    t0 = time.time()
    by_space = defaultdict(list)
    structs = []
    total = 0
    with ProcessPoolExecutor(N_PROC) as ex:
        futs = {ex.submit(scan_file, f, spaces, targets): f for f in files}
        for fut in as_completed(futs):
            name, n, kept, st = fut.result()
            total += n
            for rec in kept:
                by_space[rec['space']].append(rec)
            structs.extend(st)
            print(f'  {name}: {n} entries scanned, {len(kept)} kept, {len(st)} at target compositions '
                  f'({time.time() - t0:.0f} s)', flush=True)
    print(f'{total} entries scanned in {time.time() - t0:.0f} s')

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for old in OUT_DIR.glob('*.json'):
        old.unlink()
    for old in OUT_DIR.glob('*/cifs/*.cif'):
        old.unlink()
    for space in sorted('-'.join(s) for s in spaces):
        entries = sorted(by_space.get(space, []), key=lambda r: (r['composition_id'], r['stability'] if r['stability'] is not None else 9))
        payload = {'space': space, 'elements': space.split('-'), 'source': RELEASE, 'entries': entries}
        (OUT_DIR / f'{space}.json').write_text(json.dumps(payload, indent=1))
    n_empty = sum(1 for s in spaces if not by_space.get('-'.join(s)))
    print(f'wrote {len(spaces)} space files to {OUT_DIR} ({n_empty} empty)')

    CIF_DIR.mkdir(parents=True, exist_ok=True)
    for old in CIF_DIR.glob('*.cif'):
        old.unlink()
    rows = []
    for formula, mat_id, sdict, rec in structs:
        struct = Structure.from_dict(sdict)
        cif = CIF_DIR / f'{formula}_{mat_id}.cif'
        CifWriter(struct).write_file(cif)
        rows.append({'release': RELEASE, 'formula': formula, 'mat_id': mat_id,
                     'e_form_eV_atom': rec['delta_e'], 'e_above_hull_eV_atom': rec['e_above_hull'],
                     'e_phase_separation_eV_atom': rec['e_phase_separation'],
                     'energy_eV_atom': rec['energy_per_atom'], 'n_sites': len(struct), 'spg': rec['spg'], 'cif': str(cif)})
        # Same file a second time in the per-space layout the dashboard links to
        # (<space>/cifs/<compact composition>_<id>_stab+<meV>meV.cif).
        mev = round(rec['e_above_hull'] * 1000)
        stab = f'stab+{mev}meV' if mev >= 0 else f'stab{mev}meV'
        space_cif = OUT_DIR / rec['space'] / 'cifs' / f"{rec['composition_id'].replace(' ', '')}_{mat_id}_{stab}.cif"
        space_cif.parent.mkdir(parents=True, exist_ok=True)
        CifWriter(struct).write_file(space_cif)
    index = pd.DataFrame(rows).sort_values(['formula', 'e_above_hull_eV_atom']).reset_index(drop=True)
    index.to_csv(INDEX_CSV, index=False)
    print(f'wrote {len(index)} structures at target compositions to {CIF_DIR} and per-space cifs/; index {INDEX_CSV}')


if __name__ == '__main__':
    main()
