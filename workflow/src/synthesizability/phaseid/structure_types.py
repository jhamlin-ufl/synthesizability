"""Structure types, framework groups and provenance of the library phases.

Every phase of a space's library (``build_library``) is given

* ``method``: the route that produced the cell (database entry,
  disordered twin, Vegard prototype, pair-merged cell, MatterSim-relaxed
  cell, PDF card, COD entry, predicted cell);
* ``origin``: where the database entry itself came from (for OQMD the
  sub-library of the entry: ICSD import, prototype decoration, the ternary
  decoration library, SQS cells, ...; read from
  ``data/external/oqmd_entry_origins.tsv``, built from the local dump);
* ``group``: a framework group.  Two cells share a group when their
  frameworks match with species ignored (pymatgen ``StructureMatcher`` with
  ``FrameworkComparator``), so the three copies of one Heusler from OQMD, MP
  and Alexandria, its disordered twin, and the P-4m2 derivative of an A15
  cell all land in the group of their framework;
* an AFLOW prototype label for the group (pymatgen's
  ``AflowPrototypeMatcher``, Encyclopedia of Crystallographic Prototypes),
  found on the highest-symmetry ordered member; lower-symmetry members are
  'ordered derivatives' of it.

Results are cached per space in ``results/searchmatch/library/<space>/structure_groups.csv``.
"""
from __future__ import annotations

import csv
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from .library import LIB_DIR, ROOT

ORIGINS_FILE = ROOT / 'databases' / 'searchmatch' / 'oqmd_entry_origins.tsv'

# OQMD sub-library (second path component under /home/oqmd/libraries) -> description
OQMD_LIBRARIES = {
    'icsd': 'ICSD structure (OQMD icsd import)',
    'icsd2': 'ICSD structure (OQMD icsd2 import)',
    'prototypes': 'OQMD prototype decoration (binary/ternary prototype library)',
    'junjie': "OQMD ternary decoration sub-library ('junjie', Trinary_Compound; hypothetical cells)",
    'griesemer': 'OQMD ternary prototype decorations (Griesemer)',
    'nate': 'OQMD SQS solid-solution cells (nate)',
    'dongshu': 'OQMD sub-library (dongshu)',
    'jiangang': 'OQMD sub-library (jiangang)',
    'dohun': 'OQMD sub-library (dohun)',
    'data_mining': 'OQMD data-mined substitutions',
    'marques': 'OQMD sub-library (marques)',
    'libraries_v1_overflow_ext': 'OQMD v1 overflow library',
    'tecca_competing': 'OQMD sub-library (tecca competing phases)',
}


def _origins():
    """entry_id -> (library, duplicate_of_id) from the dump extract."""
    if not ORIGINS_FILE.exists():
        return {}
    out = {}
    with open(ORIGINS_FILE) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            path = row.get('path') or ''
            m = re.search(r'libraries/([^/]+)', path)
            lib = m.group(1) if m else path.split('/')[2] if path.count('/') > 2 else path
            dup = row.get('duplicate_of_id') or ''
            out[str(row['entry_id'])] = (lib, dup if dup not in ('NULL', '') else None)
    return out


def method_of(phase):
    """Generation route of a library phase, as text."""
    k, src = phase.key, phase.source
    if k.endswith('_dis') and phase.parent:
        return f'disordered twin of {phase.parent} (every site at the overall composition)'
    if src == 'proto':
        return 'Vegard random-occupancy prototype cell (bcc/fcc/hcp at the nominal composition)'
    if src == 'pm':
        return 'pair-merged Vegard cell between two known isostructural binaries'
    if src == 'p2':
        if '_twin_' in k:
            return 'sublattice-merged parent of an OQMD decoration (MatterSim-relaxed)'
        if '_random' in k:
            return 'random prototype solution (MatterSim-relaxed)'
        return 'MatterSim-relaxed disordered cell'
    if src == 'pdf':
        return 'ICDD PDF card (known phase)'
    if src == 'cod':
        return 'COD entry (known phase)'
    if src == 'pred':
        return 'predicted cell of the selection list'
    return 'database entry'


def origin_of(phase, origins):
    src = phase.source
    if src == 'oqmd':
        lib, dup = origins.get(str(phase.entry_id), (None, None))
        text = OQMD_LIBRARIES.get(lib, f'OQMD sub-library ({lib})' if lib else 'OQMD entry')
        if phase.icsd:
            text += '; ICSD flag'
        if dup:
            text += f'; duplicate of OQMD {dup}'
        return f'OQMD {phase.entry_id}: {text}'
    if src == 'alex':
        return f'Alexandria {phase.entry_id}'
    if src == 'mp':
        return f'Materials Project {phase.entry_id}'
    if src == 'pdf':
        return f'PDF card {phase.entry_id}'
    if src == 'cod':
        return f'COD {phase.entry_id}'
    if src == 'p2':
        return f'relaxed disordered cell {phase.entry_id}'
    if src == 'pred':
        return f'selection list {phase.entry_id}'
    return f'{src} {phase.entry_id}'


def _dummy(structure):
    """Same cell, one dummy species on every site (framework only)."""
    s = structure.copy()
    for i in range(len(s)):
        s.replace(i, 'Si')
    return s


SYMPREC = 0.3   # A; relaxed bcc-derived cells (Immm, 4-site primitive at 0.1 A) hash to bcc at 0.3, C37 and tP30 cells stay intact


def _framework_key(structure):
    """(space group of the species-free cell, primitive site count) at the
    loose tolerance, so that a relaxed, slightly distorted decoration hashes
    to its parent framework."""
    d = _dummy(structure)
    for prec in (SYMPREC, 0.5):
        try:
            sga = SpacegroupAnalyzer(d, symprec=prec, angle_tolerance=8)
            return sga.get_space_group_number(), len(sga.find_primitive())
        except Exception:
            continue
    return 0, len(d)


def _primitive(structure):
    """Primitive cell at the loose tolerance (real species kept), the
    structure itself when spglib fails."""
    for prec in (SYMPREC, 0.5):
        try:
            return SpacegroupAnalyzer(structure, symprec=prec, angle_tolerance=8).find_primitive()
        except Exception:
            continue
    return structure


_AFLOW = None


def _aflow_library():
    """[(reduced prototype structure, anonymised formula, n sites, tags)], loaded once."""
    global _AFLOW
    if _AFLOW is None:
        from pymatgen.analysis.prototypes import AFLOW_PROTOTYPE_LIBRARY, AflowPrototypeMatcher
        lib = []
        for dct in AFLOW_PROTOTYPE_LIBRARY:
            s = AflowPrototypeMatcher._preprocess_structure(dct['snl'].structure)
            lib.append((s, s.composition.anonymized_formula, len(s), dct['tags']))
        _AFLOW = lib
    return _AFLOW


def _aflow_label(structure):
    """(aflow, pearson, mineral, strukturbericht) of the best AFLOW prototype
    match, or Nones.  Same algorithm as pymatgen's AflowPrototypeMatcher
    (anonymous fit, tolerances shrunk until one prototype remains) but the
    library is loaded once and only prototypes with the same anonymised
    formula and primitive site count are tried."""
    from pymatgen.analysis.prototypes import AflowPrototypeMatcher
    from pymatgen.analysis.structure_matcher import StructureMatcher
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            red = AflowPrototypeMatcher._preprocess_structure(structure)
            anon, n = red.composition.anonymized_formula, len(red)
            cands = [(s, tags) for s, a, m, tags in _aflow_library() if a == anon and m == n]
            if not cands:
                return None, None, None, None
            sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True)
            hits = [tags for s, tags in cands if sm.fit_anonymous(s, red, skip_structure_reduction=True)]
            while len(hits) > 1 and sm.ltol >= 0.01:
                sm.ltol *= 0.8; sm.stol *= 0.8; sm.angle_tol *= 0.8
                hits = [tags for s, tags in cands if sm.fit_anonymous(s, red, skip_structure_reduction=True)]
    except Exception:
        return None, None, None, None
    if not hits:
        return None, None, None, None
    t = hits[0]
    return t.get('aflow'), t.get('pearson'), t.get('mineral') or None, t.get('strukturbericht') or None


def _match_into(buckets, fk, s, matcher, n_groups):
    """Group id of structure ``s`` with framework key ``fk``: the first bucket
    representative it matches, else a new group."""
    for g, rep in buckets.get(fk, []):
        try:
            if matcher.fit(rep, s):
                return g, n_groups
        except Exception:
            continue
    buckets.setdefault(fk, []).append((n_groups, s))
    return n_groups, n_groups + 1


def group_label(ms, tries=4):
    """AFLOW label of a framework group from its members [(key, disordered,
    spg, structure)]: ordered members tried in order of decreasing symmetry,
    more species first (AFLOW's library rarely holds the species-blind
    elemental frameworks), then smaller cells."""
    ordered = sorted([m for m in ms if not m[1] and m[3] is not None],
                     key=lambda m: (-m[2], -len(m[3].composition.elements), len(m[3])))
    lab = (None, None, None, None)
    for m in ordered[:tries]:
        lab = _aflow_label(m[3])
        if lab[0]:
            break
    return lab


def framework_groups(phases, space, cache=True, verbose=True):
    """DataFrame, one row per phase: group, framework_spg, nprim, aflow,
    pearson, mineral, strukturbericht, role, method, origin, icsd.

    Every cell, disordered ones included, goes through the matcher on its
    own: a fully disordered twin reduces to the one-site primitive of its
    framework, so the twins of two different ordered parents, the Vegard
    prototype cell and a random cell at the same composition land
    in one group (inheriting the parent's group would split them).

    The cache ``results/searchmatch/library/<space>/structure_groups.csv`` is the
    union over every sample of the space: cells of other samples (their
    prototype, predicted and pair-merged cells) are re-matched into the new
    groups and keep their metadata, so the table never loses keys.
    """
    out_dir = LIB_DIR / space
    f = out_dir / 'structure_groups.csv'
    keys = [p.key for p in phases]
    old = None
    if cache and f.exists():
        old = pd.read_csv(f)
        if set(old['key']) >= set(keys):
            return old.set_index('key').loc[keys].reset_index()
    from pymatgen.analysis.structure_matcher import FrameworkComparator, StructureMatcher
    matcher = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, comparator=FrameworkComparator(),
                               primitive_cell=False, scale=True, attempt_supercell=False)
    origins = _origins()
    buckets, group_of, fkeys, n_groups = {}, {}, {}, 0
    structures = {}
    for p in phases:
        try:
            s = Structure.from_file(str(p.cif))
        except Exception:
            s = None
        structures[p.key] = s
        if s is None:
            group_of[p.key], fkeys[p.key] = n_groups, (0, 0)
            n_groups += 1
            continue
        fk = _framework_key(s)
        fkeys[p.key] = fk
        group_of[p.key], n_groups = _match_into(buckets, fk, _primitive(s), matcher, n_groups)
    # cells of the space that belong to other samples: re-match, keep metadata
    extra = {}
    if old is not None:
        for r in old[~old['key'].isin(set(keys))].to_dict('records'):
            cif = out_dir / f"{r['key']}.cif"
            if not cif.exists():
                continue
            try:
                s = Structure.from_file(str(cif))
            except Exception:
                continue
            fk = _framework_key(s)
            g, n_groups = _match_into(buckets, fk, _primitive(s), matcher, n_groups)
            r.update(group=g, framework_spg=fk[0], nprim=fk[1])
            extra[r['key']] = (r, s, bool(r.get('role') == 'disordered member'), int(r.get('spg', 0)))
    # AFLOW label of each group from its highest-symmetry ordered member
    members = {}
    for p in phases:
        s = structures[p.key]
        members.setdefault(group_of[p.key], []).append((p.key, bool(p.disordered or (s is not None and not s.is_ordered)), p.spg, s))
    for k, (r, s, dis, spg) in extra.items():
        members.setdefault(r['group'], []).append((k, dis, spg, s))
    labels = {g: group_label(ms) for g, ms in members.items()}
    tops = {g: max((m[2] for m in ms if not m[1]), default=0) for g, ms in members.items()}

    def role(dis, spg, g):
        if dis:
            return 'disordered member'
        if spg < tops[g]:
            return f'ordered derivative (spg {spg} of {tops[g]})'
        return 'representative symmetry'

    rows = []
    for p in phases:
        g = group_of[p.key]
        fspg, nprim = fkeys[p.key]
        aflow, pearson, mineral, sb = labels[g]
        s = structures[p.key]
        dis = bool(p.disordered or (s is not None and not s.is_ordered))
        rows.append(dict(key=p.key, formula=p.formula, source=p.source, spg=p.spg, group=g, framework_spg=fspg, nprim=nprim,
                         aflow=aflow, pearson=pearson, mineral=mineral, strukturbericht=sb, role=role(dis, p.spg, g),
                         method=method_of(p), origin=origin_of(p, origins), icsd=bool(p.icsd), n_in_group=len(members[g])))
    for k, (r, s, dis, spg) in extra.items():
        g = r['group']
        aflow, pearson, mineral, sb = labels[g]
        r.update(aflow=aflow, pearson=pearson, mineral=mineral, strukturbericht=sb, role=role(dis, spg, g), n_in_group=len(members[g]))
        rows.append(r)
    df = pd.DataFrame(rows)
    if cache:
        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(f, index=False)
    if verbose:
        print(f'{space}: {len(phases)} phases (+{len(extra)} of other samples) in {n_groups} framework groups, '
              f'{sum(1 for g in labels if labels[g][0])} with an AFLOW prototype label')
    return df.set_index('key').loc[keys].reset_index()


def type_label(row):
    """Short structure-type text for a groups row: AFLOW label or framework symmetry."""
    if isinstance(row.get('aflow'), str) and row['aflow']:
        name = row['mineral'] or row['strukturbericht'] or ''
        lab = f"{name + ', ' if name else ''}{row['pearson']} ({row['aflow']})"
    else:
        lab = f"framework spg {int(row['framework_spg'])}, {int(row['nprim'])} sites per primitive cell"
    if str(row.get('role', '')).startswith('ordered derivative'):
        lab += f"; ordered derivative in spg {int(row['spg'])}"
    elif row.get('role') == 'disordered member':
        lab += '; disordered'
    return lab


def describe(key, groups, with_members=True, max_members=6):
    """One-line provenance and structure-type description of a library phase."""
    g = groups.set_index('key') if 'key' in groups.columns else groups
    if key not in g.index:
        return key
    r = g.loc[key]
    parts = [type_label(r), r['method'], r['origin']]
    if with_members:
        same = g[(g['group'] == r['group']) & (g.index != key)]
        if len(same):
            items = []
            for k2, r2 in same.head(max_members).iterrows():
                tag = r2['origin'].split(':')[0]
                items.append(f"{r2['formula']} {tag}{' (ICSD)' if r2['icsd'] else ''}{' dis' if r2['role'] == 'disordered member' else ''}")
            more = f' and {len(same) - max_members} more' if len(same) > max_members else ''
            parts.append('same framework: ' + ', '.join(items) + more)
    return '; '.join(parts)


def assemblage_groups(phases_str, groups):
    """Sorted tuple of framework-group ids of a ' + '-joined phase list, with
    'd' appended to disordered members, so that Mn2VSi from three databases
    collapse to one assemblage while ordered and disordered stay apart."""
    g = groups.set_index('key') if 'key' in groups.columns else groups
    out = []
    for k in str(phases_str).split(' + '):
        if k in g.index:
            r = g.loc[k]
            out.append(f"{int(r['group'])}{'d' if r['role'] == 'disordered member' else ''}")
        else:
            out.append(k)
    return tuple(sorted(out))


def collapse(df, groups, phases_col='phases', score_col='wR'):
    """Best row per structural assemblage; ``variants`` counts the rows it
    replaced and ``variant_keys`` lists their phase strings."""
    d = df.copy()
    d['assemblage_groups'] = [assemblage_groups(p, groups) for p in d[phases_col]]
    d = d.sort_values(score_col)
    best = d.groupby('assemblage_groups', sort=False).head(1).copy()
    counts = d.groupby('assemblage_groups').size()
    others = d.groupby('assemblage_groups')[phases_col].apply(lambda s: ' | '.join(s.iloc[1:4]))
    best['variants'] = best['assemblage_groups'].map(counts)
    best['variant_keys'] = best['assemblage_groups'].map(others)
    return best.sort_values(score_col)
