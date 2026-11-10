"""Supplementary tables of refined structural parameters of the three new compounds.

  in:  samples/<formula>/fits/<button>_structure/uiso_*_summary.json and work/fit.gpx
       (written by workflow/scripts/structure_table_refine.py, run from each main fit)
  out: results/tables/structure_table.tex   (crystal data and refinement, one column per compound)
       results/tables/structure_atoms_table.tex   (atomic coordinates and U_iso)
"""
import json
from pathlib import Path

ROWS = [  # formula, button, structure type, label of the main phase in the project, secondary-phase label for the table
    ('GdPt2Ge', '0480', r'YPd$_2$Si', r'\ch{Gd2O3}'),
    ('TiPt5Si', '0493', r'Pt$_5$InSi', 'fcc'),
    ('ZrPt5Si', '0561', r'Pt$_5$InSi', r'\ch{ZrPt3}'),
]
ATOM_NAME = {'Gd': 'Gd', 'Ge': 'Ge', 'Pt': 'Pt', 'Ti': 'Ti', 'Zr': 'Zr', 'Si': 'Si'}


def esd_fmt(val, esd, nd=None):
    """1.23456 with esd 0.00034 -> '1.2346(3)'; fixed values get nd decimals."""
    if esd is None or esd == 0:
        return f'{val:.{nd if nd is not None else 4}f}'
    d = 0
    while round(esd * 10 ** d) < 2 and d < 8:
        d += 1
    if round(esd * 10 ** d) >= 20:
        d = max(d - 1, 0)
    return f'{val:.{d}f}({round(esd * 10 ** d):d})'


def frac(v):
    for num, den in [(0, 1), (1, 2), (1, 4), (3, 4)]:
        if abs(v - num / den) < 1e-6:
            return '0' if num == 0 else f'{num}/{den}'
    return None


def wyckoff(sg, x, y, z):
    """Wyckoff letters for the two space groups of these tables (P4/mmm and Pnma, standard settings)."""
    fx, fy, fz = frac(x), frac(y), frac(z)
    if sg.replace(' ', '') == 'P4/mmm':
        table = {('0', '0', '0'): '1a', ('0', '0', '1/2'): '1b', ('1/2', '1/2', '0'): '1c',
                 ('1/2', '1/2', '1/2'): '1d', ('0', '1/2', '1/2'): '2e', ('1/2', '0', '1/2'): '2e',
                 ('0', '1/2', '0'): '2f', ('1/2', '0', '0'): '2f'}
        if (fx, fy, fz) in table:
            return table[(fx, fy, fz)]
        if fz is None:
            return {('0', '0'): '2g', ('1/2', '1/2'): '2h'}.get((fx, fy), '4i')
    if sg.replace(' ', '') == 'Pnma':
        return '4c' if fy in ('1/4', '3/4') else '8d'
    return '?'


def coord(v, e):
    f = frac(v)
    return f if (e is None and f is not None) else esd_fmt(v, e)


# Which refinement each table column shows: GdPt2Ge with refined U_iso; TiPt5Si and ZrPt5Si as the
# converged main fit with U_iso held (these data do not determine U_iso of the Pt5Si compounds reliably).
SHIFT_C = {'ZrPt5Si'}   # ZrPt5Si was refined with Zr at (1/2,1/2,0); tabulated in the TiPt5Si setting
STEP = {'GdPt2Ge': 'uiso_site', 'TiPt5Si': 'converged', 'ZrPt5Si': 'converged'}
# Weight fractions: main fit, with the spread among alternative refinements (second-phase and profile
# models) as the uncertainty, which is larger than the statistical one.
WT = {'GdPt2Ge': ('94(1)', '6(1)'), 'TiPt5Si': ('91(2)', '9(2)'), 'ZrPt5Si': ('97(1)', '3(1)')}


def load(formula, button):
    d = Path('samples') / formula / 'fits' / f'{button}_structure'
    return json.loads((d / f'{STEP[formula]}_summary.json').read_text()), d


def main():
    cols, atom_lines = [], []
    for formula, button, stype, second in ROWS:
        r, d = load(formula, button)
        k = r['Rvals']['GOF']          # standard uncertainties are scaled by GOF (misfit beyond counting statistics)
        conv = json.loads((d / 'converged_summary.json').read_text())
        kw = conv['Rvals']['GOF']
        ph = r['phases'][0]
        sc = lambda t: (t[0], t[1] * k if t[1] else t[1])
        for a in ph['atoms']:
            for key in ('x', 'y', 'z', 'uiso'):
                a[key] = sc(a[key])
            if formula in SHIFT_C:   # same origin as TiPt5Si: shift by c/2, then z -> 1 - z (mirror) to keep z < 1/2
                z0, ze = a['z']
                z1 = (z0 + 0.5) % 1.0
                a['z'] = ((1.0 - z1) % 1.0 if z1 > 0.5 + 1e-9 else z1, ze)
        c, e = ph['cell'], {key: v * k for key, v in ph['cell_esd'].items()}
        R, res = r['Rvals'], r['residuals']
        # weight fractions from the main fit (converged, U_iso as in the main fit), so they match the paper's text
        wts = {n: (v[0], v[1] * kw if v[1] else v[1]) for n, v in conv['weight_fractions'].items()}
        names = list(wts)
        wmain = wts[names[0]] if names else (None, None)
        wsec = wts[names[1]] if len(names) > 1 else (None, None)
        mass = {'GdPt2Ge': 620.04, 'TiPt5Si': 1051.44, 'ZrPt5Si': 1094.76}[formula]
        z = {'GdPt2Ge': 4, 'TiPt5Si': 1, 'ZrPt5Si': 1}[formula]
        rho = mass * z / (c['volume'] * 0.602214)
        tetr = ph['sg'].replace(' ', '') == 'P4/mmm'
        cols.append(dict(
            formula=rf'\ch{{{formula}}}', stype=stype, sg='$' + ph['sg'].replace(' ', '').replace('/', '/') + '$', z=str(z),
            a=esd_fmt(c['length_a'], e['length_a']), b='' if tetr else esd_fmt(c['length_b'], e['length_b']),
            c=esd_fmt(c['length_c'], e['length_c']), V=esd_fmt(c['volume'], e['volume']), rho=f'{rho:.2f}',
            button=button, nobs=f"{int(R['Nobs'])}", nref=f"{int(res.get('0:0:Nref', 0))}", nvar=str(r['nvars']),
            Rp=f"{res['R']:.2f}", Rwp=f"{res['wR']:.2f}", Rexp=f"{res['wRmin']:.2f}", chi2=f"{R['GOF'] ** 2:.1f}", gof=f"{R['GOF']:.2f}",
            RF=f"{res.get('0:0:Rf', float('nan')):.2f}",
            wmain=WT[formula][0], wsec=f'{second}, {WT[formula][1]}'))
        refined_types = {a['type'] for a in ph['atoms'] if a['uiso'][1]}
        for a in ph['atoms']:
            x, y, zc = a['x'], a['y'], a['z']
            u = a['uiso']
            if u[1]:
                ustr = esd_fmt(u[0], u[1])
            else:
                ustr = f'{u[0]:.3f}$^{{*}}$'
            atom_lines.append(rf"\ch{{{formula}}} & {ATOM_NAME.get(a['type'], a['type'])} & {wyckoff(ph['sg'], x[0], y[0], zc[0])} & "
                              rf"{coord(*x)} & {coord(*y)} & {coord(*zc)} & 1 & {ustr} \\")
        atom_lines.append(r'\addlinespace')
    atom_lines = atom_lines[:-1]

    def row(label, key):
        return label + ' & ' + ' & '.join(col[key] for col in cols) + r' \\'
    lines = [r'\begin{tabular}{l' + 'c' * len(cols) + '}', r'\toprule',
             row('', 'formula'), r'\midrule',
             row('Structure type', 'stype'), row('Space group', 'sg'), row('$Z$', 'z'),
             row(r'$a$ (\AA)', 'a'), row(r'$b$ (\AA)', 'b'), row(r'$c$ (\AA)', 'c'), row(r'$V$ (\AA$^3$)', 'V'),
             row(r'$\rho_\mathrm{calc}$ (g\,cm$^{-3}$)', 'rho'), r'\midrule',
             row('Button', 'button'), row('Data points', 'nobs'), row('Reflections (main phase)', 'nref'),
             row('Refined parameters', 'nvar'), row(r'$R_\mathrm{p}$ (\%)', 'Rp'), row(r'$R_\mathrm{wp}$ (\%)', 'Rwp'),
             row(r'$R_\mathrm{exp}$ (\%)', 'Rexp'), row(r'$\chi^2$', 'chi2'), row(r'GOF $= R_\mathrm{wp}/R_\mathrm{exp}$', 'gof'), row(r'$R_F$ (main phase, \%)', 'RF'),
             row(r'Main phase (wt\%)$^{a}$', 'wmain'), row(r'Second phase (wt\%)$^{a}$', 'wsec'), r'\bottomrule', r'\end{tabular}']
    out = Path('results/tables')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'structure_table.tex').write_text('\n'.join(lines) + '\n')
    alines = [r'\begin{tabular}{lllllllc}', r'\toprule',
              r'Compound & Atom & Site & $x$ & $y$ & $z$ & Occ. & $U_\mathrm{iso}$ (\AA$^2$) \\', r'\midrule',
              *atom_lines, r'\bottomrule', r'\end{tabular}']
    (out / 'structure_atoms_table.tex').write_text('\n'.join(alines) + '\n')
    print('wrote', out / 'structure_table.tex', out / 'structure_atoms_table.tex')


if __name__ == '__main__':
    main()
