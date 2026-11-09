#!/usr/bin/env python3
"""Table 1 (outcomes) and Supplementary Table S (main phases).

  in:  tables/{paper_samples,outcomes}.csv, computed/p_disorder.csv, results/hull_energies.csv
  out: results/tables/sample_table.tex, results/tables/sample_phases_table.tex

One row per target, grouped by outcome category and ordered within a category
by decreasing P_disorder.  Both are complete RevTeX table* environments.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, 'workflow/src')
import pandas as pd

from synthesizability.formula import ordered
from synthesizability.samples import CATEGORY_NAME, load_targets

OUT = Path('results/tables')
REPORTED_MARK = {'yes': r'$^{*}$', 'composition': r'$^{\dagger}$', 'series': r'$^{\ddagger}$', 'no': ''}

CAPTION = (
    'Arc-melted targets grouped by XRD outcome category, ordered within a category by decreasing $P_{\\mathrm{disorder}}$. The outcome is from the best-characterized scan (after annealing where stated). List: candidate list the target came from (A, B or C; the lists are described in Methods). $\\Delta H_{\\mathrm{f}}-C(x)$: energy of the lowest entry at the target composition relative to the convex hull $C(x)$ of all other compositions~\\cite{Gautier2015} in OQMD, the Materials Project (MP) and Alexandria (PBE), in meV/atom; negative means predicted stable, blank means no entry. $P_{\\mathrm{disorder}}$: predicted probability of a chemically disordered structure~\\cite{Jakob2025}. A structure is on the hull if it is the lowest-energy structure in at least one database that places it below the hull. $^{*}$: structure previously reported at the target composition; $^{\\dagger}$: composition reported without a structure; $^{\\ddagger}$: predicted ordered member of a previously reported series whose known member has the same sites partly disordered.'
)
PHASES_CAPTION = (
    'Main phases of the arc-melted targets of the outcome table of the main text, in the same order, with '
    'weight fractions in percent from whole-pattern fitting of the as-cast scan unless stated. '
    'Ranges span the competitive fits; (known) marks a phase with an experimental record '
    '(ICSD, PDF or COD), (known type) a known structure type whose ternary composition has not '
    'been reported separately.'
)

# --- free text -> LaTeX -----------------------------------------------------
SPACE_GROUPS = {'Pnma', 'Cmmm', 'Cmcm', 'I4mm', 'P4/nmm', 'P4/mmm', 'P6_3/mmc', 'Pbam',
                'Amm2', 'I4/mmm', 'Fm-3m', 'Im-3m'}
STRUKTURBERICHT = {'B2': 'B2', 'A2': 'A2', 'A1': 'A1', 'C15': 'C15', 'L1_2': 'L1$_2$',
                   'D0_24': 'D0$_{24}$', 'D0_22': 'D0$_{22}$', 'C11b': 'C11$_b$'}
SYMBOLS = {'Rwp': r'$R_{wp}$', 'tau-boride': r'$\tau$-boride', 'eta': r'$\eta$',
           'deg': r'$^{\circ}$', '<=': r'$\le$', '>=': r'$\ge$'}
FORMULA = re.compile(r'^[A-Z][A-Za-z0-9().,~_]*$')


def is_formula(tok):
    return bool(FORMULA.match(tok)) and any(c.isdigit() for c in tok)


def tex_formula(tok):
    """'Y3Cu4Sn4' -> 'Y$_{3}$Cu$_{4}$Sn$_{4}$', 'GdFe~0.5Ge2' -> 'GdFe$_{\\sim0.5}$Ge$_{2}$'."""
    s = tok.replace('_x', '$_x$')
    s = re.sub(r'(?<=[A-Za-z)])(\d+(?:\.\d+)?)', r'$_{\1}$', s)
    s = re.sub(r'~(\d+(?:\.\d+)?)', r'$_{\\sim\1}$', s)
    return s.replace('~', r'$\sim$')


def tex_token(tok):
    lead = tok[:len(tok) - len(tok.lstrip('('))]
    tok = tok[len(lead):]
    core = tok.rstrip('(),;:+/')
    trail = tok[len(core):]
    if core in SPACE_GROUPS:
        core = '$' + re.sub(r'_(\d)', r'_{\1}', core.replace('-3', r'\bar{3}')) + '$'
    elif core in STRUKTURBERICHT:
        core = STRUKTURBERICHT[core]
    elif core.endswith('-type') and core[:-5] in STRUKTURBERICHT:
        core = STRUKTURBERICHT[core[:-5]] + '-type'
    elif core.endswith('-type') and is_formula(core[:-5]):
        core = tex_formula(core[:-5]) + '-type'
    elif is_formula(core):
        core = tex_formula(core)
    elif '-' in core and any(c.isdigit() for c in core) and all(FORMULA.match(p) for p in core.split('-')):
        core = '-'.join(tex_formula(p) for p in core.split('-'))       # 'Au2Al-Au2Mn'
    elif core in SYMBOLS:
        core = SYMBOLS[core]
    elif core.startswith('~'):
        core = r'$\sim$' + core[1:]
    return lead + core + trail


def tex_text(text):
    text = re.sub(r'(\d) A\b', r'\1~\\AA', text.replace('%', r'\%'))
    return ' '.join(tex_token(t) for t in text.split(' ')).replace('&', r'\&').replace('#', r'\#')


def cell(text, width):
    """One line, or two inside a makecell: broken at ' | ' if the text has one,
    else at the space nearest the middle when the text is longer than `width`."""
    if ' | ' in text:
        first, second = text.split(' | ', 1)
        return r'\makecell[tl]{' + tex_text(first) + r' \\ ' + tex_text(second) + '}'
    if len(text) > width:
        cut = max((i for i, c in enumerate(text) if c == ' ' and i <= len(text) // 2 + 10), default=None)
        if cut is not None:
            return r'\makecell[tl]{' + tex_text(text[:cut]) + r' \\ ' + tex_text(text[cut + 1:]) + '}'
    return tex_text(text)


def hull(v):
    return '' if pd.isna(v) else f'{v:.1f}'


def formula_cell(r):
    return r'\ch{' + ordered(r['formula']) + '}' + REPORTED_MARK[r['reported']]


def table(spec, header, caption, label, body, float_spec='!t', vspace='4pt', tabcolsep='5pt'):
    return '\n'.join([rf'\begin{{table*}}[{float_spec}]', r'\centering', r'\footnotesize',
                      r'\caption{' + caption + '}', r'\label{' + label + '}', rf'\vspace{{{vspace}}}',
                      rf'\setlength{{\tabcolsep}}{{{tabcolsep}}}', r'\renewcommand{\arraystretch}{1.0}',
                      rf'\begin{{tabular}}{{{spec}}}', r'\toprule', *header, r'\midrule', *body,
                      r'\bottomrule', r'\end{tabular}', r'\end{table*}']) + '\n'


def _tidy(text):
    """Empty cells as '& &', as in the manuscript."""
    import re as _re
    while _re.search(r'&\s{2,}&', text):
        text = _re.sub(r'&\s{2,}&', '& &', text)
    return _re.sub(r'&\s{2,}\\\\', r'& \\\\', text)


def main():
    df = load_targets().merge(pd.read_csv('results/hull_energies.csv').drop(columns='formula'), on='sample_number')
    main_rows, phase_rows = [], []
    for cat, sub in df.groupby('category', sort=True):
        if main_rows:
            main_rows.append(r'\midrule')
            phase_rows.append(r'\midrule')
        main_rows.append(rf'\multicolumn{{8}}{{l}}{{\textbf{{{CATEGORY_NAME[cat]}}}}} \\')
        for _, r in sub.iterrows():
            main_rows.append(f"{r['row']} & {formula_cell(r)} & {r['list']} & {hull(r['oqmd_hull_meV'])} & "
                             f"{hull(r['mp_hull_meV'])} & {hull(r['alex_pbe_hull_meV'])} & "
                             f"{r['disorder_probability']:.2f} & {r['outcome_tex'] if isinstance(r.get('outcome_tex'), str) and r['outcome_tex'] else cell(r['outcome'], 70)} \\\\")
            phase_rows.append(f"{r['row']} & {formula_cell(r)} & {cell(r['phases'], 80)} \\\\")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'sample_table.tex').write_text(_tidy(table(
        'rllrrrrl',
        [r' & & & \multicolumn{3}{c}{$\Delta H_{\mathrm{f}}-C(x)$ (meV/atom)} & & \\', r'\cmidrule(lr){4-6}',
         r'\# & Formula & List & OQMD & MP & Alex. & $P_{\mathrm{disorder}}$ & XRD outcome \\'],
        CAPTION, 'tab:xrd_outcomes', main_rows, float_spec='!tp', vspace='2pt', tabcolsep='4pt')))
    (OUT / 'sample_phases_table.tex').write_text(_tidy(table(
        'rll', [r'\# & Formula & Main phases (wt\%) \\'], PHASES_CAPTION, 'tab:xrd_phases', phase_rows)))
    print(f'{len(df)} targets -> {OUT}/sample_table.tex, sample_phases_table.tex')


if __name__ == '__main__':
    main()
