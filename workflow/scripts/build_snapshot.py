#!/usr/bin/env python3
"""Build databases/snapshot/ from the database caches of the lab repository.

This is how the frozen database extract was made; it is not part of `snakemake all`
and needs the lab repository (SOURCE, default ~/repos/synthesizability) and, for the
List C PBEsol extract, the Alexandria file convex_hull_ps_2023.12.29.json.bz2.

For the 34 targets it copies, per database, the cached entries of the 22 ternary
systems and all their binary and elemental subsystems (104 files each), plus the
structures the database table compares: the OQMD prediction of each target, the
MP and Alexandria PBE entries at the target compositions.  It writes
databases/snapshot/targets.csv and a checksum manifest.  Run from the repository root:

    python workflow/scripts/build_snapshot.py            # databases + structures
    python workflow/scripts/build_snapshot.py --list-c   # also the PBEsol 2023.12.29 extract (needs ijson)
"""
from __future__ import annotations

import argparse
import bz2
import hashlib
import os
import shutil
import subprocess
from itertools import combinations
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

SOURCE = Path(os.environ.get('SOURCE', Path.home() / 'repos/synthesizability'))
OUT = Path('databases/snapshot')
DATABASES = {   # snapshot folder: source cache
    'oqmd_v1.7': 'data/external/oqmd_ternary_phases',
    'mp': 'data/external/mp_ternary_phases',
    'alexandria_pbe_2025.07.02': 'data/external/alexandria_pbe_ternary_phases',
    'alexandria_pbesol_optimade': 'data/external/alexandria_pbesol_ternary_phases',
}
PBESOL_HULL = Path.home() / 'data/alexandria/convex_hull_ps_2023.12.29.json.bz2'
LIST_C_CSV = 'alexandria_pbesol_2023.12.29_list_c.csv'
LIST_C = ['B2NbNi2', 'C2CrSc', 'C2LaRu', 'CoSiY', 'CoSi3Y', 'LaNiSi3', 'Al3Cu2Y', 'CuSnY']


def system(formula):
    return '-'.join(sorted(e.symbol for e in Composition(formula).elements))


def subsystems(space):
    els = space.split('-')
    return ['-'.join(c) for n in (1, 2, 3) for c in combinations(els, n)]


def last_commit_date(path):
    out = subprocess.run(['git', '-C', str(SOURCE), 'log', '-1', '--format=%ad', '--date=short', '--', path],
                         capture_output=True, text=True).stdout.strip()
    return out


def targets():
    cats = pd.read_csv(SOURCE / 'results/phaseid/review/categories_option_T.csv')
    synth = pd.read_csv(SOURCE / 'data/processed/synthesis_data.csv')
    t = cats[['sample_number']].merge(synth[['sample_number', 'formula', 'oqmd_entry_id']], on='sample_number')
    t['oqmd_entry_id'] = t['oqmd_entry_id'].astype(int)
    t['system'] = t['formula'].map(system)
    return t.sort_values('sample_number')


def copy_databases(t, manifest):
    spaces = sorted({s for sp in t['system'] for s in subsystems(sp)})
    for folder, src in DATABASES.items():
        (OUT / folder).mkdir(parents=True, exist_ok=True)
        for sp in spaces:
            rel = f'{src}/{sp}.json'
            shutil.copy2(SOURCE / rel, OUT / folder / f'{sp}.json')
            manifest.append((f'{folder}/{sp}.json', rel, last_commit_date(rel)))
    print(f'{len(spaces)} subsystems x {len(DATABASES)} databases')


def copy_structures(t, manifest):
    d = OUT / 'structures'
    for sub in ('oqmd', 'mp', 'alexandria_pbe'):
        (d / sub).mkdir(parents=True, exist_ok=True)
    for _, r in t.iterrows():
        rel = f'data/external/oqmd_structures/{r.formula}/{r.oqmd_entry_id}.cif'
        name = f'oqmd/{r.formula}_{r.oqmd_entry_id}.cif'
        shutil.copy2(SOURCE / rel, d / name)
        manifest.append((f'structures/{name}', rel, last_commit_date(rel)))
        red = Composition(r.formula).reduced_formula
        mpdir = SOURCE / 'data/external/mp_ternary_phases' / r.system / 'cifs'
        for cif in sorted(mpdir.glob('*.cif')) if mpdir.exists() else []:
            if Composition(cif.name.split('_')[0]).reduced_formula == red:
                rel = str(cif.relative_to(SOURCE))
                shutil.copy2(cif, d / 'mp' / cif.name)
                manifest.append((f'structures/mp/{cif.name}', rel, last_commit_date(rel)))
    idx = pd.read_csv(SOURCE / 'data/external/alexandria_hull_structures/index.csv')
    idx = idx[(idx['release'] == 'pbe_2025.07.02') & idx['formula'].isin(t['formula'])].copy()
    for _, r in idx.iterrows():
        name = Path(r['cif']).name
        shutil.copy2(SOURCE / r['cif'], d / 'alexandria_pbe' / name)
        manifest.append((f'structures/alexandria_pbe/{name}', r['cif'], last_commit_date(r['cif'])))
    idx['cif'] = idx['cif'].map(lambda p: f'databases/snapshot/structures/alexandria_pbe/{Path(p).name}')
    idx.drop(columns='release').to_csv(d / 'alexandria_pbe_index.csv', index=False)


def extract_list_c_pbesol(manifest):
    """Entries of the 2023.12.29 PBEsol convex-hull file in the List C systems
    (the hull the List C selection used), without structures."""
    import ijson
    keep = {frozenset(s.split('-')) for f in LIST_C for s in subsystems(system(f))}
    rows = []
    with bz2.open(PBESOL_HULL, 'rb') as fh:
        for e in ijson.items(fh, 'entries.item', use_float=True):
            if frozenset(e['composition']) in keep:
                d = e['data']
                rows.append(dict(mat_id=d['mat_id'], formula=d['formula'],
                                 composition=' '.join(f'{k}{v:g}' for k, v in e['composition'].items()),
                                 spg=d['spg'], nsites=d['nsites'], energy_eV=e['energy'],
                                 e_form_eV_atom=d['e_form'], e_above_hull_eV_atom=d['e_above_hull'],
                                 e_phase_separation_eV_atom=d['e_phase_separation']))
    out = OUT / LIST_C_CSV
    pd.DataFrame(rows).sort_values(['formula', 'mat_id']).to_csv(out, index=False)
    manifest.append((out.name, str(PBESOL_HULL), 'file dated 2026-09-28'))
    print(f'{len(rows)} PBEsol 2023.12.29 entries in the List C systems')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--list-c', action='store_true')
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t = targets()
    t.to_csv(OUT / 'targets.csv', index=False)
    manifest = []
    copy_databases(t, manifest)
    copy_structures(t, manifest)
    if args.list_c:
        extract_list_c_pbesol(manifest)
    elif (OUT / LIST_C_CSV).exists():   # keep the extract of an earlier --list-c run
        manifest.append((LIST_C_CSV, str(PBESOL_HULL), 'file dated 2026-09-28'))
    seen, rows = set(), []
    for path, src, date in manifest:
        if path in seen:
            continue
        seen.add(path)
        sha = hashlib.sha256((OUT / path).read_bytes()).hexdigest()
        rows.append(dict(path=path, source=src, source_last_commit=date, sha256=sha))
    pd.DataFrame(rows).sort_values('path').to_csv(OUT / 'MANIFEST.csv', index=False)
    print(f'{len(rows)} files in {OUT}/MANIFEST.csv')


if __name__ == '__main__':
    main()
