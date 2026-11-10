"""Candidate phase library for one chemical space.

Sources (``data/external/<db>_ternary_phases/<space>/cifs``), all three with
element, binary and ternary coverage once ``extract_alexandria_space_cifs.py``
and ``extract_mp_space_cifs.py`` have run:

* OQMD cache (``oqmd``)
* Alexandria PBE 2025.07.02 cache (``alex``)
* Materials Project cache (``mp``)
* ``data/temp/predicted_structures_cifs/<sample>_<formula>``: the structures
  the paper calls "predicted" (``pred``; no energy of their own)

Formation energies are only comparable within one database, so each phase
keeps its source and the hull geometry is evaluated per source
(``mixtures.enumerate_all``).

Every CIF in the caches is written in P1 by pymatgen.  Jade and GSAS-II then
refine six free cell parameters per phase, which wrecked the student's first
multiphase batch, so each structure is symmetrised here (spglib through
pymatgen) and rewritten with its real space group before anything else sees
it.  For every ordered phase with at least two species a fully disordered
parent (all sites carry the average composition) is added as well, because the
solid-solution outcomes of the study are disordered variants of ordered
database structures.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

import numpy as np
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Composition, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ..paths import ROOT, SEARCH_DIR, LIBRARY_DIR, PDF_CARD_DIR
DB_DIRS = {
    'oqmd': SEARCH_DIR / 'oqmd_ternary_phases',
    'alex': SEARCH_DIR / 'alexandria_pbe_ternary_phases',
    'mp': SEARCH_DIR / 'mp_ternary_phases',
}
PRED_DIR = SEARCH_DIR / 'predicted_structures'
LIB_DIR = LIBRARY_DIR

CU_KA1 = 1.540598
SYMPRECS = (0.05, 0.1, 0.01)


@dataclass
class Phase:
    key: str                    # unique within the library, used in file names
    source: str                 # oqmd | alex | mp | pred
    entry_id: str
    composition: Composition
    x: np.ndarray               # atomic fractions in the space's element order
    e_form: float | None        # eV/atom, on the source database's scale
    stab: float | None          # hull distance in the source database, eV/atom
    icsd: bool
    cif: Path                   # symmetrised CIF
    spg: int
    nsites: int
    at_target: bool
    disordered: bool = False
    parent: str | None = None   # key of the ordered ordered decoration of a disordered parent
    d_hkl: np.ndarray = field(default_factory=lambda: np.zeros(0), repr=False)
    intensity: np.ndarray = field(default_factory=lambda: np.zeros(0), repr=False)

    @property
    def formula(self):
        return self.composition.reduced_formula

    @property
    def order(self):
        return len(self.composition.elements)


def space_elements(formula):
    return sorted(el.symbol for el in Composition(formula).elements)


def fractions(comp, elements):
    return np.array([comp.get_atomic_fraction(el) for el in elements])


def symmetrise(structure):
    """Conventional standard cell with the detected space group, or None."""
    for symprec in SYMPRECS:
        try:
            sga = SpacegroupAnalyzer(structure, symprec=symprec)
            refined = sga.get_refined_structure()
            return refined, sga.get_space_group_number(), symprec
        except Exception:
            continue
    return None, None, None


def write_symmetrised(structure, path, symprec):
    path.parent.mkdir(parents=True, exist_ok=True)
    CifWriter(structure, symprec=symprec).write_file(str(path))


def disordered_twin(structure):
    """Every site carries the overall composition; None if already disordered."""
    if not structure.is_ordered or len(structure.composition.elements) < 2:
        return None
    occ = {sp: f for sp, f in structure.composition.fractional_composition.items()}
    s = structure.copy()
    for i in range(len(s)):
        s.replace(i, occ)
    return s


def sticks(structure, two_theta_range=(5, 135)):
    """d-spacings and relative intensities (max 100) for Cu K-alpha1."""
    pat = XRDCalculator(wavelength=CU_KA1).get_pattern(structure, two_theta_range=two_theta_range)
    y = np.asarray(pat.y, float)
    if y.size and y.max() > 0:
        y = y / y.max() * 100
    return np.asarray(pat.d_hkls, float), y


def _save_cached(out_dir, key, spg, nsites, d, y):
    np.savez(out_dir / f'{key}.npz', spg=spg, nsites=nsites, d_hkl=d, intensity=y)


def _load_cached(out_dir, key):
    """(spg, nsites, d_hkl, intensity) if both the symmetrised CIF and the
    stick pattern of ``key`` are cached, else None."""
    npz, cif = out_dir / f'{key}.npz', out_dir / f'{key}.cif'
    if not (npz.exists() and cif.exists()):
        return None
    try:
        z = np.load(npz)
        return int(z['spg']), int(z['nsites']), z['d_hkl'], z['intensity']
    except Exception:
        # Truncated by an interrupted build; rebuild this phase
        npz.unlink(missing_ok=True)
        return None


def db_entries(source, elements):
    """(source, id, Composition, e_form, stab, icsd, cif) for every cached
    entry of the space and its subsystems that has a CIF."""
    d = DB_DIRS[source]
    out = []
    for n in (1, 2, 3):
        for combo in combinations(elements, n):
            space = '-'.join(combo)
            f = d / f'{space}.json'
            if not f.exists():
                continue
            cif_dir = d / space / 'cifs'
            if not cif_dir.exists():
                continue
            cifs = {}
            for p in cif_dir.glob('*.cif'):
                m = re.match(r'.+?_(.+)_stab', p.name)
                if m:
                    cifs[m.group(1)] = p
            for e in json.load(open(f))['entries']:
                eid = str(e.get('entry_id') or e.get('mp_id'))
                cif = cifs.get(eid)
                if cif is None:
                    continue
                out.append((source, eid, Composition(e['composition_id']), e.get('delta_e'),
                            e.get('stability'), bool(e.get('icsd')), cif))
    return out


# MatterSim-relaxed disordered structures: sublattice-merged (partially disordered) parents of
# OQMD decorations and prototype solutions, cells relaxed with MatterSim.
# Index CSV columns: space, name, composition, entry_id (source OQMD entry or
# blank), pair, e_above_hull (eV/atom, OQMD scale) or dE_ord_meV, cif (path).
P2_DIR = Path(os.environ.get('PHASEID_P2_DIR', str(SEARCH_DIR / 'partially_disordered')))
P2_CELL_SCALE = 0.98   # MatterSim cells run about 2 % large against the measured ones
# Formation-energy scale a source shares for hull geometry: PDF-card phases
# are given the OQMD hull energy at their composition, so they join OQMD's hull;
# The relaxed pair-merged parents carry OQMD formation energy plus their ordering energy.
SOURCE_GROUPS = {'oqmd': ('oqmd', 'pdf', 'cod', 'p2', 'proto', 'pm'), 'alex': ('alex',), 'mp': ('mp',)}


def _pdf_cards(elements):
    """ICDD PDF-card CIFs supplied by the user (SYNTH_PDF_CARD_DIR) whose elements lie in the space.

    The cards are licensed by the ICDD and are not part of this repository; without
    SYNTH_PDF_CARD_DIR the library is built from the open sources only.

    These are known phases the student tried in Jade; several (Y3Cu4Sn4,
    Nb2Ni21B6) exist in no DFT database, so without them the search cannot
    find them.  Energies are assigned later (OQMD hull at the composition).
    """
    out, seen = [], set()
    els = set(elements)
    if PDF_CARD_DIR is None or not PDF_CARD_DIR.is_dir():
        return out
    for cif in sorted(PDF_CARD_DIR.rglob('*.cif')):
        card = cif.stem.replace('PDF Card - ', '')
        if card in seen:
            continue
        try:
            s = Structure.from_file(str(cif))
        except Exception:
            continue
        if len(s) == 0 or not {e.symbol for e in s.composition.elements} <= els:
            continue
        seen.add(card)
        out.append(('pdf', card, s.composition, None, 0.0, True, cif))
    return out


# ---------------------------------------------------------------------------
# Crystallography Open Database entries: experimental cells with
# partial occupancies, downloaded per space by workflow/scripts/cod_download.py
# ---------------------------------------------------------------------------
COD_DIR = SEARCH_DIR / 'cod'
COD_MAX_SITES = 48
COD_DISORDERED = {}      # key -> has partial occupancy (from the index), for the cached branch


def _cod_entries(elements):
    """COD structures of the space as (source 'cod', id, Composition, None, 0.0, True, cif).

    Near-duplicate determinations (same reduced formula, space-group number
    and cell lengths within 0.5 %) collapse onto the first one in the index;
    giant cells (> COD_MAX_SITES sites) are skipped.  Energies are assigned
    later like the PDF cards (OQMD hull at the composition).
    """
    import pandas as pd
    space = '-'.join(elements)
    index = COD_DIR / space / 'index.csv'
    if not index.exists():
        return []
    df = pd.read_csv(index)
    out, seen = [], []
    for r in df.itertuples(index=False):
        cif = COD_DIR / space / str(r.file)
        if not cif.exists():
            continue
        try:
            s = Structure.from_file(str(cif))
        except Exception:
            continue
        if len(s) == 0 or len(s) > COD_MAX_SITES or not {e.symbol for e in s.composition.elements} <= set(elements):
            continue
        sig = (s.composition.reduced_formula, int(r.sg_number) if r.sg_number == r.sg_number else 0)
        abc = np.array(sorted([s.lattice.a, s.lattice.b, s.lattice.c]))
        # elements: one determination per space group (the many Pt, Fe, Ni entries
        # are pressure and temperature series); compounds: cells within 1 % collapse
        tol = 1.0 if len(s.composition.elements) == 1 else 0.01
        if any(sg == sig and np.all(np.abs(abc / cell - 1) < tol) for sg, cell in seen):
            continue
        seen.append((sig, abc))
        key = f'{s.composition.reduced_formula}_cod_{r.cod_id}'
        COD_DISORDERED[key] = bool(not s.is_ordered)
        out.append(('cod', str(r.cod_id), s.composition, None, 0.0, True, cif))
    return out


def _project2_disordered(elements, oqmd_entries):
    """The relaxed disordered structures for this space, as
    (source 'p2', name, Composition, e_form on the OQMD scale, stab, icsd,
    cif).  The cell is scaled by P2_CELL_SCALE when loaded (see
    ``build_library``); e_form = OQMD formation energy of the source entry +
    ordering energy when both are known, else the OQMD hull energy at the
    composition (prior 0, as for PDF cards)."""
    index = P2_DIR / 'index.csv'
    if not index.exists():
        return []
    import pandas as pd
    from .mixtures import hull_at
    df = pd.read_csv(index)
    space = '-'.join(elements)
    df = df[df['space'].astype(str) == space]
    if df.empty:
        return []
    by_id = {str(e[1]): e for e in oqmd_entries}
    Xo = np.array([fractions(e[2], elements) for e in oqmd_entries])
    Eo = np.array([e[3] for e in oqmd_entries])
    out = []
    for r in df.to_dict('records'):
        cif = Path(str(r.get('file') or r.get('cif')))
        if not cif.is_absolute():
            cif = P2_DIR / cif
        if not cif.exists():
            cif = P2_DIR / space / cif.name
            if not cif.exists():
                continue
        comp = Composition(str(r['composition']))
        try:
            hull = hull_at(Xo, Eo, fractions(comp, elements))[0]
        except Exception:
            hull = None
        eah = r.get('e_above_oqmd_hull_eV', r.get('e_above_hull'))
        de = r.get('dE_ord_meV')
        src = by_id.get(str(r.get('source_entry_id', r.get('entry_id', ''))))
        if eah is not None and not pd.isna(eah) and hull is not None:
            e_form = hull + float(eah)          # their value is already relative to the OQMD hull
        elif src is not None and de is not None and not pd.isna(de):
            e_form = src[3] + float(de) / 1000
        else:
            e_form = hull                       # prototype solutions: prior 0, like PDF cards
        stab = float(eah) if eah is not None and not pd.isna(eah) else 0.0
        # Only the prototype solutions have MatterSim-relaxed random cells that
        # run 2 % large; the pair-merged parents carry the relaxed volume of
        # the ordered entry and must not be rescaled
        scale = str(r.get('kind', '')).lower() != 'twin'
        P2_SCALE_BY_NAME[str(r['name'])] = scale
        out.append(('p2', str(r['name']), comp, e_form, stab, False, cif))
    return out


P2_SCALE_BY_NAME = {}

PROTOTYPES = {   # name: (space group number, fractional coordinates of the one site, atoms per cell, c/a)
    'bcc_A2': (229, [[0.0, 0.0, 0.0]], 2, None),
    'fcc_A1': (225, [[0.0, 0.0, 0.0]], 4, None),
    'hcp_A3': (194, [[1 / 3, 2 / 3, 1 / 4]], 2, 1.633),
}


def _element_volumes(raw, elements):
    """Atomic volume per element (A^3) from the most stable elemental entry in
    the raw database list, falling back to a sphere of the atomic radius."""
    from pymatgen.core import Element
    vol = {}
    for el in elements:
        cands = sorted([r for r in raw if r[0] in ('oqmd', 'alex', 'mp') and r[2].reduced_formula == el],
                       key=lambda r: (r[4] if r[4] is not None else 1.0))
        for r in cands:
            try:
                s = Structure.from_file(str(r[6]))
                vol[el] = s.volume / len(s)
                break
            except Exception:
                continue
        if el not in vol:
            vol[el] = 4 / 3 * np.pi * float(Element(el).atomic_radius or 1.5) ** 3 / 0.74
    return vol


# ---------------------------------------------------------------------------
# Pair-merged parents: interpolated solid solutions between
# isostructural known binaries that share one element
# ---------------------------------------------------------------------------
PM_FRACTIONS = (0.25, 0.5, 0.75)


def _pair_merged(elements, raw, out_dir, verbose=True):
    """Vegard-interpolated cells along the tie line between two experimentally
    known binaries of the same structure type that share one element with the
    same stoichiometry (PtGe/NiGe -> (Pt,Ni)Ge, ScPd/ScCu -> Sc(Pd,Cu), ...).
    The library's twins keep each entry's own composition; these are the
    cells in between, which no database contains.  Known = OQMD icsd flag or a
    PDF card.  Pairs are cached in ``out_dir/pair_merged.json``; the cells join
    the OQMD hull at the hull energy of their composition (prior 0, like PDF
    cards).  Returns (source, id, Composition, e_form, stab, icsd, cif) rows."""
    import json
    from itertools import combinations
    from pymatgen.analysis.structure_matcher import FrameworkComparator, StructureMatcher
    from pymatgen.core import Lattice
    from .mixtures import hull_at
    pm_dir = out_dir / 'pm'
    pm_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / 'pair_merged.json'
    oq = [r for r in raw if r[0] == 'oqmd' and r[3] is not None]
    Xo = np.array([fractions(r[2], elements) for r in oq]) if oq else None
    Eo = np.array([r[3] for r in oq]) if oq else None
    if cache.exists():
        pairs = json.loads(cache.read_text())
    else:
        known = [r for r in raw if r[0] in ('oqmd', 'pdf') and (r[5] or r[0] == 'pdf') and len(r[2].elements) == 2]
        structs = {}
        for r in known:
            try:
                s = Structure.from_file(str(r[6]))
            except Exception:
                continue
            if len(s) > 24:      # dilute 64-atom cells and other big decorations are not solution end members
                continue
            structs[str(r[6])] = (r, s)
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, comparator=FrameworkComparator(), primitive_cell=True, scale=True, attempt_supercell=False)
        pairs = []
        seen = set()
        for (c1, (r1, s1)), (c2, (r2, s2)) in combinations(structs.items(), 2):
            e1, e2 = set(r1[2].elements), set(r2[2].elements)
            shared = e1 & e2
            if len(shared) != 1 or r1[2].reduced_formula == r2[2].reduced_formula:
                continue
            sh = next(iter(shared))
            if abs(r1[2].get_atomic_fraction(sh) - r2[2].get_atomic_fraction(sh)) > 1e-3:
                continue
            tag = tuple(sorted([r1[2].reduced_formula, r2[2].reduced_formula]))
            if tag in seen:
                continue
            try:
                if sm.fit(s1, s2):
                    pairs.append(dict(f1=r1[2].reduced_formula, f2=r2[2].reduced_formula, cif1=c1, cif2=c2, shared=str(sh)))
                    seen.add(tag)
            except Exception:
                continue
        cache.write_text(json.dumps(pairs, indent=1))
        if verbose and pairs:
            print(f'  pair-merged parents: {len(pairs)} isostructural known-binary pairs ({", ".join(p["f1"] + "/" + p["f2"] for p in pairs)})')
    out = []
    for p in pairs:
        try:
            s1, s2 = Structure.from_file(p['cif1']), Structure.from_file(p['cif2'])
        except Exception:
            continue
        sh = p['shared']
        o1 = [str(e) for e in s1.composition.elements if str(e) != sh][0]
        o2 = [str(e) for e in s2.composition.elements if str(e) != sh][0]
        # primitive cells of both for the Vegard volume; the cell of s1 carries the sites
        v1, v2 = s1.volume / len(s1), s2.volume / len(s2)
        for x in PM_FRACTIONS:
            v = (1 - x) * v1 + x * v2
            s = s1.copy()
            s.scale_lattice(v * len(s))
            for i, site in enumerate(s):
                if str(site.specie) == o1:
                    s.replace(i, {o1: 1 - x, o2: x})
            comp = s.composition
            name = f"{p['f1']}-{p['f2']}_x{int(round(100 * x)):02d}"
            cif = pm_dir / f'{name}.cif'
            if not cif.exists():
                CifWriter(s).write_file(str(cif))
            e_form = None
            if Xo is not None:
                try:
                    e_form, _ = hull_at(Xo, Eo, fractions(comp, elements))
                except Exception:
                    e_form = None
            out.append(('pm', name, comp, e_form, 0.0, False, cif))
    return out


def _prototype_cells(elements, target, raw, oqmd_entries):
    """Random-occupancy bcc, fcc and hcp cells at the nominal composition,
    cell size from Vegard's rule on the elemental volumes (source 'proto',
    prior 0 like the PDF cards: the hull energy at the composition).  They
    are the solid-solution hypotheses the databases never contain."""
    from pymatgen.core.lattice import Lattice
    from .mixtures import hull_at
    space = '-'.join(elements)
    out_dir = LIB_DIR / space / 'proto'
    out_dir.mkdir(parents=True, exist_ok=True)
    x = {el: target.get_atomic_fraction(el) for el in target.elements}
    vols = _element_volumes(raw, elements)
    v = sum(x[el] * vols[el.symbol] for el in x)
    species = {el: float(f) for el, f in x.items() if f > 0}
    Xo = np.array([fractions(e[2], elements) for e in oqmd_entries])
    Eo = np.array([e[3] for e in oqmd_entries])
    try:
        e_hull = hull_at(Xo, Eo, fractions(target, elements))[0]
    except Exception:
        e_hull = None
    out = []
    for name, (spg, coords, n, c_over_a) in PROTOTYPES.items():
        if c_over_a is None:
            a = (n * v) ** (1 / 3)
            lat = Lattice.cubic(a)
        else:
            a = (n * v / (c_over_a * np.sqrt(3) / 2)) ** (1 / 3)
            lat = Lattice.hexagonal(a, c_over_a * a)
        from pymatgen.symmetry.groups import SpaceGroup
        sg = SpaceGroup.from_int_number(spg)
        s = Structure.from_spacegroup(sg.symbol, lat, [species], coords)
        label = f'{target.reduced_formula}_{name}'
        cif = out_dir / f'{label}.cif'
        if not cif.exists():
            s.to(filename=str(cif))
        out.append(('proto', name, target, e_hull, 0.0, False, cif))
    return out


def _predicted(sample_number, target):
    out = []
    for d in PRED_DIR.glob(f'{sample_number:04d}_*'):
        for cif in sorted(d.glob('*.cif')):
            m = re.search(r'hull([+-]?\d+)meV', cif.name)
            stab = int(m.group(1)) / 1000 if m else None
            eid = cif.stem.split('_')[0]
            out.append(('pred', eid, Composition(target), None, stab, False, cif))
    return out


def build_library(formula, sample_number=None, sources=('oqmd', 'alex', 'mp'),
                  stab_cutoff=0.25, twins=True, pdf_cards=True, p2=True, prototypes=True, pair_merged=True, cod=True, verbose=True):
    """Symmetrised phases of the chemical space of ``formula``.

    ``stab_cutoff`` drops individual phases more than that far above their own
    database's hull (eV/atom); the mixture ranking later applies its own prior.
    Returns (elements, phases).
    """
    elements = space_elements(formula)
    target = Composition(formula)
    space = '-'.join(elements)
    out_dir = LIB_DIR / space
    raw = []
    for source in sources:
        raw += db_entries(source, elements)
    if sample_number is not None:
        raw += _predicted(sample_number, target)
    if pdf_cards:
        cards = _pdf_cards(elements) + (_cod_entries(elements) if cod else [])
        if cards and 'oqmd' in sources:
            # Known phases join OQMD's hull at the hull energy of their composition
            from .mixtures import hull_at
            oq = [r for r in raw if r[0] == 'oqmd']
            Xo = np.array([fractions(r[2], elements) for r in oq])
            Eo = np.array([r[3] for r in oq])
            fixed = []
            for src, card, comp, _, stab, icsd, cif in cards:
                try:
                    e_hull, _ = hull_at(Xo, Eo, fractions(comp, elements))
                except Exception:
                    e_hull = None
                fixed.append((src, card, comp, e_hull, stab, icsd, cif))
            cards = fixed
        raw += cards
    if p2 and 'oqmd' in sources:
        raw += _project2_disordered(elements, [r for r in raw if r[0] == 'oqmd'])
    if prototypes and 'oqmd' in sources:
        raw += _prototype_cells(elements, target, raw, [r for r in raw if r[0] == 'oqmd'])
    if pair_merged and 'oqmd' in sources:
        raw += _pair_merged(elements, raw, out_dir, verbose=verbose)

    phases, seen = [], set()
    at_target = lambda c: c.reduced_formula == target.reduced_formula  # noqa: E731
    for source, eid, comp, e_form, stab, icsd, cif in raw:
        if stab is not None and stab > stab_cutoff:
            continue
        key = f'{comp.reduced_formula}_{source}_{eid}'
        if key in seen:
            continue
        seen.add(key)
        x = fractions(comp, elements)
        cached = _load_cached(out_dir, key)
        if cached is not None:
            # Both the phase and (if it exists) its twin were cached together
            spg, nsites, d, y = cached
            phases.append(Phase(key, source, eid, comp, x, e_form, stab, icsd, out_dir / f'{key}.cif',
                                spg, nsites, at_target(comp), disordered=(source in ('proto', 'pm') or COD_DISORDERED.get(key, False)), d_hkl=d, intensity=y))
            tw = _load_cached(out_dir, key + '_dis')
            if twins and tw is not None:
                phases.append(Phase(key + '_dis', source, eid, comp, x, e_form, stab, icsd,
                                    out_dir / f'{key}_dis.cif', tw[0], tw[1], at_target(comp),
                                    disordered=True, parent=key, d_hkl=tw[2], intensity=tw[3]))
            continue
        try:
            s = Structure.from_file(str(cif))
            if source == 'p2' and P2_SCALE_BY_NAME.get(eid, True):
                s.scale_lattice(s.volume * P2_CELL_SCALE ** 3)
        except Exception as exc:
            if verbose:
                print(f'  skip {cif.name}: {exc}')
            continue
        sym, spg, symprec = symmetrise(s)
        if sym is None:
            if verbose:
                print(f'  skip {cif.name}: symmetrisation failed')
            continue
        path = out_dir / f'{key}.cif'
        write_symmetrised(sym, path, symprec)
        d, y = sticks(sym)
        _save_cached(out_dir, key, spg, len(sym), d, y)
        phases.append(Phase(key, source, eid, comp, x, e_form, stab, icsd, path, spg, len(sym),
                            at_target(comp), disordered=(source in ('proto', 'pm') or not sym.is_ordered), d_hkl=d, intensity=y))
        if twins and source not in ('proto', 'pm') and sym.is_ordered:
            tw = disordered_twin(sym)
            if tw is None:
                continue
            tsym, tspg, tprec = symmetrise(tw)
            if tsym is None:
                continue
            tkey = key + '_dis'
            tpath = out_dir / f'{tkey}.cif'
            write_symmetrised(tsym, tpath, tprec)
            td, ty = sticks(tsym)
            _save_cached(out_dir, tkey, tspg, len(tsym), td, ty)
            phases.append(Phase(tkey, source, eid, comp, x, e_form, stab, icsd, tpath,
                                tspg, len(tsym), at_target(comp), disordered=True, parent=key,
                                d_hkl=td, intensity=ty))
    if verbose:
        n_dis = sum(p.disordered for p in phases)
        by_src = {s: sum(1 for p in phases if p.source == s and not p.disordered)
                  for s in list(sources) + ['pdf', 'pred', 'p2']}
        print(f'{space}: {len(phases)} phases ({len(phases) - n_dis} ordered {by_src}, '
              f'{n_dis} disordered parents)')
    return elements, phases
