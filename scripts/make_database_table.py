#!/usr/bin/env python3
"""Supplementary table: how the three databases see each target.

  results/publication_ready/database_table.tex
  results/publication_ready/database_agreement.csv

Per target and per database: the recomputed hull distance (negative = the
database predicts the compound stable; see compute_hull_distances.py) and
whether the database's own ground-state structure at that composition is the
same structure as the OQMD prediction the target was selected from.

Structures: MP from the cached CIFs at the target compositions
(data/external/mp_ternary_phases/<space>/cifs); Alexandria from every entry at
the target compositions in the PBE 2025.07.02 release
(extract_alexandria_release_phases.py), on or off the hull.  Structure identity is pymatgen's StructureMatcher
with ltol 0.1, stol 0.2, angle_tol 4 deg.
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

sys.path.insert(0, 'src')

import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Composition, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from synthesizability.formula import sort_by_electronegativity
from synthesizability.paper_samples import OUTCOME_LABEL, OUTCOMES, load_classified_samples

HULL_CSV = Path('data/processed/hull_distances.csv')
ALEX_INDEX = Path('data/external/alexandria_hull_structures/index.csv')
ALEX_RELEASE = 'pbe_2025.07.02'
OUT_TEX = Path('results/publication_ready/database_table.tex')
OUT_CSV = Path('results/publication_ready/database_agreement.csv')
LABEL = 'tab:databases'

MATCHER = StructureMatcher(ltol=0.1, stol=0.2, angle_tol=4)

CAPTION = (
    'The {n} targets as seen by the three databases. For each database, '
    r'$\Delta E_{\mathrm{hull}}$ is the energy of its lowest entry at the target composition relative to the '
    'convex hull of all other compositions in that database (meV/atom; negative means the database predicts '
    'the compound to be stable), recomputed from the cached entries by the same construction for all three. '
    '``Same'' states whether the database\'s ground-state structure at that composition is the same structure '
    'as the OQMD prediction from which the target was selected (pymatgen StructureMatcher); the space group is '
    'given where it differs. Blank: no entry at the composition. '
    'Snapshots: OQMD v1.7 local dump; Materials Project API, April 2026; Alexandria PBE 2025.07.02 release.'
)


def spacegroup(struct):
    try:
        return SpacegroupAnalyzer(struct, symprec=0.05).get_space_group_symbol()
    except Exception:
        return '?'


def tex_sg(symbol):
    """'P6_3mc' -> '$P6_3mc$', 'Fm-3m' -> '$Fm\\bar{3}m$'."""
    s = re.sub(r'-(\d)', r'\\bar{\1}', symbol)
    return f'${s}$'


def mp_ground_state(formula):
    els = sorted(e.symbol for e in Composition(formula).elements)
    red = Composition(formula).reduced_formula
    best = None
    for cif in glob.glob(f'data/external/mp_ternary_phases/{"-".join(els)}/cifs/*.cif'):
        m = re.match(r'([A-Za-z0-9]+)_(mp-\d+)_stab([+-]\d+)meV', Path(cif).name)
        if not m or Composition(m.group(1)).reduced_formula != red:
            continue
        stab = int(m.group(3))
        if best is None or stab < best[0]:
            best = (stab, m.group(2), cif)
    return best


def agreement(pred, other_cif):
    other = Structure.from_file(other_cif)
    same = MATCHER.fit(pred, other)
    return same, spacegroup(other)


def main():
    df = load_classified_samples()
    hull = pd.read_csv(HULL_CSV).set_index('formula')
    alex = pd.read_csv(ALEX_INDEX) if ALEX_INDEX.exists() else pd.DataFrame(columns=['release', 'formula', 'cif', 'e_form_eV_atom'])
    alex = alex[alex['release'] == ALEX_RELEASE]

    rows = []
    for _, r in df.iterrows():
        formula = r['formula']
        pred = Structure.from_file(f'data/external/oqmd_structures/{formula}/{int(r["oqmd_entry_id"])}.cif')
        h = hull.loc[formula]
        row = {'sample': int(r['sample_number']), 'formula': formula, 'outcome': r['XRD Result'],
               'oqmd_hull_meV': h['oqmd_hull_meV'], 'oqmd_spacegroup': spacegroup(pred),
               'mp_hull_meV': h['mp_hull_meV'], 'mp_same': None, 'mp_spacegroup': None, 'mp_id': None,
               'alex_hull_meV': h['alex_pbe_hull_meV'], 'alex_same': None, 'alex_spacegroup': None, 'alex_id': None}
        mp = mp_ground_state(formula)
        if mp is not None:
            same, sg = agreement(pred, mp[2])
            row.update(mp_same=same, mp_spacegroup=sg, mp_id=mp[1])
        sub = alex[alex['formula'] == formula]
        if not sub.empty:
            best = sub.sort_values('e_form_eV_atom').iloc[0]
            same, sg = agreement(pred, best['cif'])
            row.update(alex_same=same, alex_spacegroup=sg, alex_id=best['mat_id'])
        rows.append(row)
    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)

    def hull_cell(v):
        return '' if pd.isna(v) else f'{v:.1f}'

    def same_cell(same, sg):
        if same is None or (isinstance(same, float) and pd.isna(same)):
            return ''
        return 'yes' if same else 'no (' + tex_sg(sg) + ')'

    lines = [
        r'\begin{table*}[!t]', r'\centering', r'\footnotesize',
        r'\caption{' + CAPTION.replace('{n}', str(len(df))) + '}', r'\label{' + LABEL + '}', r'\vspace{4pt}',
        r'\setlength{\tabcolsep}{5pt}',
        r'\begin{tabular}{rllrrlrl}', r'\toprule',
        r' & & \multicolumn{1}{c}{OQMD} & & \multicolumn{2}{c}{Materials Project} & \multicolumn{2}{c}{Alexandria (PBE)} \\',
        r'\cmidrule(lr){3-3}\cmidrule(lr){5-6}\cmidrule(lr){7-8}',
        r'\# & Formula & $\Delta E_{\mathrm{hull}}$ & Structure & $\Delta E_{\mathrm{hull}}$ & Same & $\Delta E_{\mathrm{hull}}$ & Same \\',
        r'\midrule',
    ]
    order = {'SS': ('disorder_probability', False), 'MP': ('oqmd_stability', False), 'P': ('oqmd_stability', False)}
    n = 0
    first = True
    for outcome in OUTCOMES:
        sub = df[df['XRD Result'] == outcome]
        if sub.empty:
            continue
        col, asc = order[outcome]
        sub = sub.sort_values([col, 'sample_number'], ascending=[asc, True])
        if not first:
            lines.append(r'\midrule')
        first = False
        lines.append(r'\multicolumn{8}{l}{\textit{' + OUTCOME_LABEL[outcome] + '}} \\\\')
        for _, r in sub.iterrows():
            o = out[out['sample'] == int(r['sample_number'])].iloc[0]
            n += 1
            lines.append(
                f"{n} & \\ch{{{sort_by_electronegativity(o['formula'])}}} & {hull_cell(o['oqmd_hull_meV'])} & "
                f"{tex_sg(o['oqmd_spacegroup'])} & {hull_cell(o['mp_hull_meV'])} & {same_cell(o['mp_same'], o['mp_spacegroup'])} & "
                f"{hull_cell(o['alex_hull_meV'])} & {same_cell(o['alex_same'], o['alex_spacegroup'])} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table*}']
    OUT_TEX.write_text('\n'.join(lines) + '\n')
    agree = out.dropna(subset=['mp_same']); alex_ok = out.dropna(subset=['alex_same'])
    print(f"MP: {int(agree['mp_same'].astype(bool).sum())}/{len(agree)} same structure; Alexandria: {int(alex_ok['alex_same'].astype(bool).sum())}/{len(alex_ok)} same")
    print(f'Saved {OUT_TEX} and {OUT_CSV}')


if __name__ == '__main__':
    main()
