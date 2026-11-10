#!/usr/bin/env python3
"""Write samples/<formula>/README.md for each of the 34 targets.

  in:  tables/paper_samples.csv, tables/outcomes.csv, tables/buttons.csv, tables/scans.csv,
       tables/si_fits.csv, samples/*/fits/*/summary.json, databases/snapshot/targets.csv,
       computed/p_disorder.csv, computed/phonons/dynamic_flags_5M_paper.csv
  out: samples/<formula>/README.md

Run from the repository root: python workflow/scripts/sample_readmes.py
"""
import csv
import json
from pathlib import Path

ROLE = {'FoR': 'main fit', 'C1': 'contrast', 'C2': 'contrast', 'C3': 'contrast', 'ascast': 'as-cast scan'}
REPORTED = {'yes': 'Structure previously reported at this composition.',
            'composition': 'Composition previously reported, without a structure.',
            'series': 'Ordered member of a reported series whose known member is partly disordered.',
            'no': 'Not previously reported.'}
# Fits whose results the Supplement text quotes but does not plot.
EXTRA = {
    '0494_pred_for_TiSiPt5_hs': 'predicted Ti2Pt8Si in place of the TiPt5Si-type phase',
    '0561_exchange_md': 'ordered ZrPt5Si with refined Zr/Si exchange (about 1 %), with texture',
    '0561_exchange_nomd': 'ordered ZrPt5Si with refined Zr/Si exchange, without texture',
    '0538_ascast_pnma': 'as-cast scan: Pnma YCoSi (known) + Y3Co2Si3 + YCo2Si2',
    '0538_ascast_pred_P4nmm': 'as-cast scan: predicted P4/nmm in place of Pnma, with the two minor phases',
}
DATABASES = ['oqmd_v1.7', 'mp', 'alexandria_pbe_2025.07.02', 'alexandria_pbesol_optimade']
SEARCH = ['oqmd_ternary_phases', 'mp_ternary_phases', 'alexandria_pbe_ternary_phases', 'cod']


def rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def rwp(fit_dir):
    s = json.loads((Path(fit_dir) / 'summary.json').read_text())
    rec = s.get('record') or s.get('rec')
    return f"{float(rec['wR']):.2f}"


def extra_label(name):
    if name in EXTRA:
        return EXTRA[name]
    return ''


def readme(t, out, buttons, scans, fits, system, oqmd_id, p_dis, flags):
    f, here = t['formula'], Path('samples') / t['formula']
    up = '../../'
    lines = [f'# {f}', '',
             f"Table 1 row {t['table_row']}, candidate list {t['list']}. "
             f"Discussed in the Supplement, Sec. S7, under {f}.", '',
             f"Outcome: {out['outcome'].replace(' | ', ' ')}", '',
             f"Category: {t['category_name']} ({t['confidence']})"
             + (f"; category {t['alt_category']} not excluded." if t['alt_category'] else '.')
             + f" {REPORTED[t['reported']]}", '',
             f"P_disorder = {float(p_dis):.2f}. Database structures below the hull that are dynamically "
             f"stable (MatterSim): {flags['n_stable_below_hull']} of {flags['n_below_hull']}.", '',
             '## Buttons', '',
             'Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) '
             'and the measured scans.', '',
             '| Button | Role | Mass loss (%) | Scans |', '|---|---|---|---|']
    for b in buttons:
        sc = [s for s in scans if s['sample'] == b['sample']]
        cell = '; '.join(f"`{Path(s['file']).name}` ({s['instrument'].split(' (')[0]}, {s['state']})" for s in sc) or 'none'
        lines.append(f"| {b['sample']} | {b['role']} | {b['mass_loss_pct'] or 'not recorded'} | {cell} |")
    rec = Path(t['record_scan']).relative_to(here)
    lines += ['', f"Scan used: `{rec}` ({t['instrument']}, {t['mount']}"
              + (', annealed 1000 C for one week).' if t['annealed'] == 'yes' else ', as cast).'), '']
    lines += ['## Fits', '',
              'Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), '
              '`summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it '
              '(see `../README.md`).', '',
              '| Folder | Role | Model | R_wp (%) |', '|---|---|---|---|']
    shown = set()
    for r in fits:
        if not r['fit_dir']:
            continue
        shown.add(r['fit_dir'])
        role = ROLE.get(r['role'], r['role'])
        if r['role'] == 'FoR' and str(t['category']) == '6':
            role = 'best fit found (not a satisfactory description; not shown in the Supplement)'
        lines.append(f"| `{Path(r['fit_dir']).name}` | {role} | {r['label']} | {rwp(r['fit_dir'])} |")
    for d in sorted((here / 'fits').glob('*')) if (here / 'fits').exists() else []:
        if str(d) not in shown:
            lines.append(f"| `{d.name}` | quoted in the Supplement text, not plotted | {extra_label(d.name)} | {rwp(d)} |")
    if not shown:
        lines.append('| | | no fit deposited | |')
    lines += ['', '## Database entries', '',
              f'Chemical system {system}. The entries are kept per system and shared with the other targets in it:', '',
              f"- hull energies: `{up}databases/snapshot/<database>/{system}.json` "
              f"({', '.join(d for d in DATABASES if (Path('databases/snapshot') / d / f'{system}.json').exists())}), "
              'with the binary and elemental subsystems in the same folders',
              f"- predicted structure (the OQMD entry the target was chosen from): "
              f"`{up}databases/snapshot/structures/oqmd/{oqmd_name(oqmd_id)}`",
              f"- search-match structures: `{up}databases/searchmatch/<source>/{system}/` "
              f"({', '.join(s for s in SEARCH if (Path('databases/searchmatch') / s / system).exists())})", '']
    return '\n'.join(lines)


def oqmd_name(oqmd_id):
    hits = sorted(Path('databases/snapshot/structures/oqmd').glob(f'*_{oqmd_id}.cif'))
    return hits[0].name if hits else f'(OQMD entry {oqmd_id})'


def main():
    out = {r['sample']: r for r in rows('tables/outcomes.csv')}
    targets = {int(r['sample_number']): r for r in rows('databases/snapshot/targets.csv')}
    p_dis = {int(r['sample_number']): r['disorder_probability'] for r in rows('computed/p_disorder.csv')}
    flags = {int(r['sample']): r for r in rows('computed/phonons/dynamic_flags_5M_paper.csv')}
    buttons, scans, fits = rows('tables/buttons.csv'), rows('tables/scans.csv'), rows('tables/si_fits.csv')
    for t in rows('tables/paper_samples.csv'):
        n = int(t['sample'])
        text = readme(t, out[t['sample']], [b for b in buttons if b['target_sample'] == t['sample']], scans,
                      [r for r in fits if r['formula'] == t['formula']], targets[n]['system'],
                      targets[n]['oqmd_entry_id'], p_dis[n], flags[n])
        (Path('samples') / t['formula'] / 'README.md').write_text(text)
    print('wrote', len(rows('tables/paper_samples.csv')), 'sample READMEs')


if __name__ == '__main__':
    main()
