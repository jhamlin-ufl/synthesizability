#!/usr/bin/env python3
"""Record the provenance of every cached database the analysis depends on.

  data/external/snapshot_provenance.json

The hull distances, structure comparisons and the OQMD-derived targets all rest
on cached copies of OQMD, the Materials Project and Alexandria that were pulled
at different times from different kinds of source.  This file pins down, for
each cache, where it came from, when, how many files it holds and a checksum of
their contents, so that a later re-pull can be diffed against this state and the
paper can state exactly which database versions it used.

Sources as of the first snapshot (2026-10-01):
  oqmd         local MySQL dump `qmdb`, OQMD v1.7, queried by src/synthesizability/oqmd.py
  mp           Materials Project REST API (mp_api), April 2026
  alex_pbe     until 2026-10-01: Alexandria OPTIMADE API, PBE (partial);
               from 2026-10-01: bulk PBE 2025.07.02 release files (complete)
  alex_pbesol  Alexandria OPTIMADE API, PBEsol, re-pulled 2026-09-28
  alexandria_local  convex-hull release files downloaded to ~/data/alexandria
                    (PBE 2025.07.02; PBE and PBEsol 2023.12.29), on-hull entries
                    with full structures

Run from the repository root.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import subprocess
from pathlib import Path

OUT = Path('data/external/snapshot_provenance.json')

CACHES = {
    'oqmd': ('data/external/oqmd_ternary_phases', 'local MySQL dump qmdb (OQMD v1.7) via src/synthesizability/oqmd.py'),
    'mp': ('data/external/mp_ternary_phases', 'Materials Project REST API (mp_api)'),
    'alex_pbe': ('data/external/alexandria_pbe_ternary_phases', 'Alexandria bulk PBE 2025.07.02 release (extract_alexandria_release_phases.py)'),
    'alex_pbesol': ('data/external/alexandria_pbesol_ternary_phases', 'Alexandria OPTIMADE API, PBEsol'),
    'oqmd_structures': ('data/external/oqmd_structures', 'CIFs of the predicted (lowest-stability) OQMD entry per target'),
}
ALEXANDRIA_DIR = Path(os.environ.get('ALEXANDRIA_DIR', os.path.expanduser('~/data/alexandria')))
ALEXANDRIA_FILES = {
    'convex_hull.json.bz2': 'Alexandria PBE 2025.07.02 convex-hull release',
    'convex_hull_pbe_2023.12.29.json.bz2': 'Alexandria PBE 2023.12.29 convex-hull release',
    'convex_hull_ps_2023.12.29.json.bz2': 'Alexandria PBEsol 2023.12.29 convex-hull release',
}


def sha256_of_files(files):
    h = hashlib.sha256()
    for f in files:
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def sha256_file(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def git(*args):
    try:
        return subprocess.run(['git', *args], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def mysql(query):
    try:
        out = subprocess.run(['mysql', '-N', '-e', f'USE qmdb; {query}'], capture_output=True, text=True,
                             check=True, timeout=120).stdout.strip()
        return out or None
    except Exception:
        return None


def describe_cache(path, source):
    root = Path(path)
    files = sorted(p for p in root.rglob('*') if p.is_file() and p.name != '.gitkeep')
    jsons = [p for p in files if p.suffix == '.json']
    mtimes = [p.stat().st_mtime for p in files]
    return {
        'path': str(root),
        'source': source,
        'n_files': len(files),
        'n_json': len(jsons),
        'bytes': sum(p.stat().st_size for p in files),
        'oldest_file': dt.datetime.fromtimestamp(min(mtimes)).date().isoformat() if mtimes else None,
        'newest_file': dt.datetime.fromtimestamp(max(mtimes)).date().isoformat() if mtimes else None,
        'last_commit': git('log', '-1', '--format=%h %cs %s', '--', str(root)),
        'sha256_json': sha256_of_files(jsons),
        'git_clean': git('status', '--porcelain', '--', str(root)) == '',
    }


def main():
    record = {
        'recorded': dt.datetime.now().isoformat(timespec='seconds'),
        'repository_commit': git('rev-parse', '--short', 'HEAD'),
        'caches': {k: describe_cache(p, s) for k, (p, s) in CACHES.items()},
        'oqmd_dump': {
            'database': 'qmdb (MySQL)',
            'version': '1.7',
            'entries': mysql('SELECT COUNT(*) FROM entries;'),
            'formation_energy_rows': mysql('SELECT COUNT(*) FROM formation_energies;'),
            'note': ('The stored `stability` column is a snapshot and can disagree with the dump\'s own '
                     'energies (Hf2MoIr: stored -18.6, recomputed -12.0 meV/atom); hull distances are '
                     'recomputed by scripts/compute_hull_distances.py.'),
        },
        'alexandria_local': {},
    }
    chunks = sorted((ALEXANDRIA_DIR / 'pbe_2025.07.02').glob('alexandria_*.json.bz2'))
    record['alexandria_local']['pbe_2025.07.02_chunks'] = {
        'description': 'Alexandria PBE 2025.07.02 full release, alexandria_*.json.bz2',
        'path': str(ALEXANDRIA_DIR / 'pbe_2025.07.02'), 'n_files': len(chunks),
        'bytes': sum(f.stat().st_size for f in chunks),
        'sha256_of_listing': hashlib.sha256('\n'.join(f'{f.name} {f.stat().st_size}' for f in chunks).encode()).hexdigest(),
    }
    for name, desc in ALEXANDRIA_FILES.items():
        f = ALEXANDRIA_DIR / name
        record['alexandria_local'][name] = (
            {'description': desc, 'path': str(f), 'bytes': f.stat().st_size,
             'file_date': dt.datetime.fromtimestamp(f.stat().st_mtime).date().isoformat(),
             'sha256': sha256_file(f)}
            if f.exists() else {'description': desc, 'path': str(f), 'missing': True})
    OUT.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ('n_json', 'newest_file', 'last_commit', 'git_clean')}
                      for k, v in record['caches'].items()}, indent=2))
    print(f'Saved {OUT}')


if __name__ == '__main__':
    main()
