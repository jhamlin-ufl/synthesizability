#!/usr/bin/env python3
"""Supplementary table: how OQMD, MP and Alexandria PBE see each target.

  in:  results/hull_energies.csv, databases/snapshot/structures/, tables/{paper_samples,outcomes}.csv
  out: results/tables/database_table.tex, results/database_agreement.csv

Per target and database: Delta H_f - C(x), and whether the database's lowest-energy
structure at the composition is the OQMD prediction the target was selected from
(pymatgen StructureMatcher, ltol 0.1, stol 0.2, angle_tol 4 deg; space group given
where it differs).  MP: lowest e_above_hull among the cached CIFs at the composition;
Alexandria: lowest formation energy among the 2025.07.02 entries at the composition.
"""
import re
import sys
import warnings
from pathlib import Path

sys.path.insert(0, 'workflow/src')
import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Composition, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from synthesizability.formula import ordered
from synthesizability.samples import CATEGORY_NAME, load_targets

STRUCT = Path('databases/snapshot/structures')
OUT_TEX = Path('results/tables/database_table.tex')
OUT_CSV = Path('results/database_agreement.csv')
warnings.filterwarnings('ignore', message='Issues encountered while parsing CIF')
MATCHER = StructureMatcher(ltol=0.1, stol=0.2, angle_tol=4)

CAPTION = (
    'The {n} candidates as seen by the three databases. For each database, '
    r'$\Delta H_{\mathrm{f}}-C(x)$ is the energy of its lowest entry at the nominal composition relative to the '
    'convex hull of all other compositions in that database (meV/atom; negative means the database predicts '
    'the compound to be stable), recomputed from the cached entries by the same construction for all three. '
    "``Same'' states whether the database's ground-state structure at that composition is the same structure "
    'as the OQMD prediction from which the candidate was selected (pymatgen StructureMatcher); the space group is '
    'given where it differs. Blank: no entry at the composition. '
    'Grouped and numbered as in Table~\\ref{tab:xrd_outcomes}. '
    'Snapshots: OQMD v1.7 local dump; Materials Project API, 7 April 2026; Alexandria PBE 2025.07.02 release.'
)


def spacegroup(s):
    try:
        return SpacegroupAnalyzer(s, symprec=0.05).get_space_group_symbol()
    except Exception:
        return '?'


def tex_sg(symbol):
    return '$' + re.sub(r'-(\d)', r'\\bar{\1}', symbol) + '$'


def mp_ground_state(formula):
    """(mp_id, cif) of the lowest-e_above_hull MP structure at the composition, or None."""
    red, best = Composition(formula).reduced_formula, None
    for cif in STRUCT.glob('mp/*.cif'):
        m = re.match(r'([A-Za-z0-9]+)_(mp-\d+)_stab([+-]\d+)meV', cif.name)
        if m and Composition(m.group(1)).reduced_formula == red and (best is None or int(m.group(3)) < best[0]):
            best = (int(m.group(3)), m.group(2), cif)
    return best and best[1:]


def alex_ground_state(alex, formula):
    sub = alex[alex['formula'] == formula]
    if sub.empty:
        return None
    best = sub.sort_values('e_form_eV_atom').iloc[0]
    return best['mat_id'], best['cif']


def _tidy(text):
    """Empty cells as '& &', as in the manuscript."""
    import re as _re
    while _re.search(r'&\s{2,}&', text):
        text = _re.sub(r'&\s{2,}&', '& &', text)
    return _re.sub(r'&\s{2,}\\\\', r'& \\\\', text)


def main():
    df = load_targets().merge(pd.read_csv('results/hull_energies.csv').drop(columns='formula'), on='sample_number')
    alex = pd.read_csv(STRUCT / 'alexandria_pbe_index.csv')
    rows = []
    for _, r in df.iterrows():
        pred = Structure.from_file(STRUCT / f'oqmd/{r.formula}_{r.oqmd_entry_id}.cif')
        row = dict(row=r.row, sample_number=r.sample_number, formula=r.formula, category=r.category,
                   oqmd_hull_meV=r.oqmd_hull_meV, oqmd_spacegroup=spacegroup(pred),
                   mp_hull_meV=r.mp_hull_meV, alex_hull_meV=r.alex_pbe_hull_meV)
        for db, gs in (('mp', mp_ground_state(r.formula)), ('alex', alex_ground_state(alex, r.formula))):
            if gs is not None:
                other = Structure.from_file(gs[1])
                row.update({f'{db}_id': gs[0], f'{db}_same': MATCHER.fit(pred, other),
                            f'{db}_spacegroup': spacegroup(other)})
        rows.append(row)
    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)

    def hull(v):
        return '' if pd.isna(v) else f'{v:.1f}'

    def same(r, db):
        v = r.get(f'{db}_same')
        return '' if pd.isna(v) else 'yes' if v else 'no (' + tex_sg(r[f'{db}_spacegroup']) + ')'

    lines = [r'\begin{table*}[!t]', r'\centering', r'\footnotesize',
             r'\caption{' + CAPTION.replace('{n}', str(len(out))) + '}', r'\label{tab:databases}',
             r'\vspace{4pt}', r'\setlength{\tabcolsep}{5pt}', r'\begin{tabular}{rllrrlrl}', r'\toprule',
             r' & & \multicolumn{1}{c}{OQMD} & & \multicolumn{2}{c}{Materials Project} & \multicolumn{2}{c}{Alexandria (PBE)} \\',
             r'\cmidrule(lr){3-3}\cmidrule(lr){5-6}\cmidrule(lr){7-8}',
             r'\# & Formula & $\Delta H_{\mathrm{f}}-C$ & Structure & $\Delta H_{\mathrm{f}}-C$ & Same & $\Delta H_{\mathrm{f}}-C$ & Same \\',
             r'\midrule']
    for cat, sub in out.groupby('category', sort=True):
        if cat != out['category'].min():
            lines.append(r'\midrule')
        lines.append(r'\multicolumn{8}{l}{\textit{' + CATEGORY_NAME[cat] + '}} \\\\')
        for _, r in sub.iterrows():
            lines.append(f"{r['row']} & \\ch{{{ordered(r['formula'])}}} & {hull(r['oqmd_hull_meV'])} & "
                         f"{tex_sg(r['oqmd_spacegroup'])} & {hull(r['mp_hull_meV'])} & {same(r, 'mp')} & "
                         f"{hull(r['alex_hull_meV'])} & {same(r, 'alex')} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table*}']
    OUT_TEX.parent.mkdir(parents=True, exist_ok=True)
    OUT_TEX.write_text(_tidy('\n'.join(lines) + '\n'))
    for db in ('mp', 'alex'):
        s = out[f'{db}_same'].dropna()
        print(f'{db}: {int(s.astype(bool).sum())}/{len(s)} same structure as the OQMD prediction')
    print(f'-> {OUT_TEX}, {OUT_CSV}')


if __name__ == '__main__':
    main()
