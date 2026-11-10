"""Table S1: the arc-melted buttons behind each candidate and the pattern used.

Reads tables/paper_samples.csv (candidates, Table 1 order, pattern used) and
tables/buttons.csv (every button, its role and mass loss). Buttons with
in_table_s1 = no (0618, a second TiPt5Si button made elsewhere) are left out.
Writes results/tables/buttons_table.tex (the tabular body of Table S1).
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results' / 'tables' / 'buttons_table.tex'


def loss(v):
    if pd.isna(v) or str(v).strip() == '':
        return '--'
    v = float(v)
    return '$<$0.1' if v < 0.05 else f'{v:.1f}'


def main():
    ps = pd.read_csv(ROOT / 'tables/paper_samples.csv', dtype=str)
    bt = pd.read_csv(ROOT / 'tables/buttons.csv', dtype=str)
    bt = bt[bt['in_table_s1'] != 'no']
    lines = [r'\# & Compound & Button & Earlier attempt & Remake & \multicolumn{2}{l}{Pattern used} \\', r'\midrule']
    cat = None
    for _, r in ps.sort_values('table_row', key=lambda c: c.astype(int)).iterrows():
        if r['category_name'] != cat:
            if cat is not None:
                lines.append(r'\midrule')
            cat = r['category_name']
            lines.append(rf'\multicolumn{{7}}{{l}}{{\textit{{{cat}}}}} \\')
        mine = bt[bt['target_sample'] == r['sample']]
        target = mine[mine['role'] == 'target']
        earlier = mine[mine['role'] == 'earlier attempt']
        later = mine[mine['role'].isin(['remake', 'second button'])]
        fmt = lambda d: ', '.join(f"{b['sample']} ({loss(b['mass_loss_pct'])})" + (r'$^{\dagger}$' if str(b['n_scans']) == '0' else '')
                                  for _, b in d.sort_values('sample').iterrows())
        state = 'annealed' if r['annealed'] == 'yes' else 'as cast'
        inst = 'PAN' if 'PANalytical' in r['instrument'] else 'D500'
        mount = 'ZB' if 'zero-background' in r['mount'] else ('PMMA' if 'PMMA' in r['mount'] else 'bulk')
        lines.append(rf"{r['table_row']} & \ch{{{r['formula']}}} & {fmt(target)} & {fmt(earlier)} & {fmt(later)} & {r['record_sample']} & "
                     rf"\parbox[t]{{0.25\textwidth}}{{\raggedright {state}, {inst}, {mount}\strut}} \\")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text('\n'.join(lines) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
