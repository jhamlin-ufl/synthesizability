"""Converge a main fit and refine isotropic displacement parameters for the structure tables.

Starts from the GSAS-II project of a main fit (copied, never modified), then
  1. converge: repeat refinement cycles with the main fit's own refinement flags until the
     largest shift/esd is below 0.01 (or --max-rounds is reached);
  2. uiso: also refine U_iso of the main phase, one value per site (--uiso site), one value shared
     by the atoms of each element (--uiso element) or one value for all its atoms (--uiso all),
     and converge again.
Writes, for each step, summary.json (R values, cell and atoms with standard uncertainties, weight
fractions read from the GSAS-II listing, correlations of each U_iso with the other refined
variables above 0.8) and a CIF of every phase.

Usage (from the repository root, GSAS-II importable via GSAS2_DIR):
  python workflow/scripts/structure_table_refine.py --gpx samples/.../work/fit.gpx --out results/structure/<name> --uiso site
"""
import argparse, json, os, re, shutil, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, os.environ.get('GSAS2_DIR', '/home/jhamlin/repos/rietveld/external/gsas2/GSAS-II'))
from GSASII import GSASIIscriptable as G2sc  # noqa: E402


def converge(g, tag, max_rounds):
    """Refine until max shift/esd < 0.01; return the R values of the last round."""
    g.data['Controls']['data']['max cyc'] = 20
    rv = {}
    for k in range(max_rounds):
        g.refine()
        rv = dict(g['Covariance']['data'].get('Rvals', {}))
        msg = rv.get('Max shft/sig', np.nan)
        print(f'{tag} round {k + 1}: Rwp {rv.get("Rwp"):.3f}  max shift/esd {msg:.3g}', flush=True)
        if msg < 0.01:
            break
    return {k: float(v) for k, v in rv.items() if np.isscalar(v)}


def weight_fractions(lst):
    """Weight fractions (value, esd) per phase from the last 'Weight fraction' lines of the listing."""
    out = {}
    if not lst.exists():
        return out
    for name, val, esd in re.findall(r'Phase:\s*(\S+).*?\n.*?Weight fraction\s*:\s*([\d.]+)(?:,\s*sig\s*([\d.]+))?',
                                     lst.read_text(errors='ignore'), re.S):
        out[name] = (float(val), float(esd) if esd else None)
    return out


def record(g, out, tag):
    cov = g['Covariance']['data']
    vary = cov.get('varyList', [])
    sig = dict(zip(vary, cov.get('sig', [])))
    cm = cov.get('covMatrix')
    res = dict(step=tag, Rvals={k: float(v) for k, v in cov.get('Rvals', {}).items() if np.isscalar(v)},
               nvars=len(vary), phases=[])
    h = g.histograms()[0]
    res['residuals'] = {k: float(v) for k, v in h.residuals.items() if isinstance(v, (int, float, np.floating))}
    for ph in g.phases():
        cell, cesd = ph.get_cell_and_esd()
        atoms = []
        for i, a in enumerate(ph.data['Atoms']):
            pid = ph.id
            def v(name, idx):
                key = f'{pid}::{name}:{i}'
                return a[idx], sig.get(key)
            x, y, z = v('dAx', 3), v('dAy', 4), v('dAz', 5)   # GSAS-II refines position shifts dAx...
            u = v('AUiso', 10)
            occ = v('Afrac', 6)
            corr = {}
            ukey = f'{pid}::AUiso:{i}'
            if cm is not None and ukey in vary:
                s = np.sqrt(np.diag(cm))
                j = vary.index(ukey)
                c = cm[j] / (s[j] * s)
                corr = {vary[k]: round(float(c[k]), 3) for k in range(len(vary)) if k != j and abs(c[k]) > 0.8}
            atoms.append(dict(label=a[0], type=a[1], x=x, y=y, z=z, frac=a[6], occ=occ, uiso=u, flags=a[2],
                              site_sym=a[7] if len(a) > 7 else None, mult=a[8] if len(a) > 8 else None,
                              uiso_corr_gt_0p8=corr))
        res['phases'].append(dict(name=ph.name, sg=ph.data['General']['SGData']['SpGrp'],
                                  cell=cell, cell_esd=cesd, atoms=atoms))
        try:
            ph.export_CIF(str(out / f'{tag}_{ph.name}.cif'))
        except Exception as e:  # CIF export is a convenience; the summary is the record
            print(f'CIF export failed for {ph.name}: {e}', flush=True)
    res['weight_fractions'] = weight_fractions(Path(g.filename).with_suffix('.lst'))
    (out / f'{tag}_summary.json').write_text(json.dumps(res, indent=1, default=float))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--gpx', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--uiso', choices=['site', 'element', 'all', 'heavy'], default='site',
                    help='heavy: refine only the U_iso of the heaviest element (shared), others held')
    ap.add_argument('--max-rounds', type=int, default=6)
    ap.add_argument('--occ', action='store_true',
                    help='instead of U_iso: refine the occupancy of every site of the main phase except the '
                         'heaviest site with the largest multiplicity, which anchors the scale (U_iso held)')
    ap.add_argument('--improve', choices=['none', 'profile', 'profile_rough'], default='none',
                    help='before refining U_iso: profile = 20-term background, refined axial-divergence asymmetry '
                         'SH/L and generalized microstrain of the main phase; profile_rough also refines the '
                         'Suortti surface-roughness correction')
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    work = a.out / 'work'
    work.mkdir(exist_ok=True)
    gpx = work / 'fit.gpx'
    shutil.copy(a.gpx, gpx)
    g = G2sc.G2Project(str(gpx))

    converge(g, 'converge', a.max_rounds)
    record(g, a.out, 'converged')

    if a.improve != 'none':
        h = g.histograms()[0]
        refs = {'Background': {'type': 'chebyschev-1', 'no. coeffs': 20, 'refine': True},
                'Instrument Parameters': ['SH/L']}
        if a.improve == 'profile_rough':
            refs['Sample Parameters'] = ['SurfRoughA', 'SurfRoughB']
            h.data['Sample Parameters']['SurfRoughA'][0] = 0.1    # start off zero; Suortti model
            h.data['Sample Parameters']['SurfRoughB'][0] = 0.1
        h.set_refinements(refs)
        g.phases()[0].set_HAP_refinements({'Mustrain': {'type': 'generalized', 'refine': True}})
        converge(g, f'improve_{a.improve}', a.max_rounds)
        record(g, a.out, f'improved_{a.improve}')

    main_ph = g.phases()[0]
    atoms = main_ph.data['Atoms']
    pid = main_ph.id
    if a.occ:
        Z = {'Pt': 78, 'Gd': 64, 'Zr': 40, 'Ge': 32, 'Ti': 22, 'Si': 14}
        ref = max(range(len(atoms)), key=lambda i: (Z.get(atoms[i][1], 0), atoms[i][8]))
        for i, at in enumerate(atoms):
            at[2] = at[2].replace('F', '').replace('U', '') + ('' if i == ref else 'F')
        print('occupancy reference site:', atoms[ref][0], flush=True)
        converge(g, 'occ', a.max_rounds)
        record(g, a.out, 'occ')
        return
    # drop existing U_iso constraints of the main phase (the Pt5 fits tie the two Pt sites); the mode sets its own
    cons = g.data['Constraints']['data']
    cons['Phase'] = [c for c in cons.get('Phase', [])
                     if not any('AUiso' in str(term[1]) for term in c[:-3] if isinstance(term, list))]
    heavy = max({at[1] for at in atoms}, key=lambda el: {'Pt': 78, 'Gd': 64, 'Zr': 40, 'Ge': 32, 'Ti': 22, 'Si': 14}.get(el, 0))
    for at in atoms:
        at[2] = at[2].replace('U', '')
        if a.uiso != 'heavy' or at[1] == heavy:
            at[2] = at[2] + 'U'
    groups = {}
    if a.uiso == 'heavy':
        groups[heavy] = [f'{pid}::AUiso:{i}' for i, at in enumerate(atoms) if at[1] == heavy]
    if a.uiso == 'element':
        for i, at in enumerate(atoms):
            groups.setdefault(at[1], []).append(f'{pid}::AUiso:{i}')
    elif a.uiso == 'all':
        groups['all'] = [f'{pid}::AUiso:{i}' for i in range(len(atoms))]
    for members in groups.values():
        if len(members) > 1:
            g.add_EquivConstr(members)
    tag = f'uiso_{a.uiso}' + ('' if a.improve == 'none' else f'_{a.improve}')
    converge(g, tag, a.max_rounds)
    r = record(g, a.out, tag)
    for at in r['phases'][0]['atoms']:
        print(at['label'], 'Uiso', at['uiso'], 'corr>0.8', at['uiso_corr_gt_0p8'], flush=True)


if __name__ == '__main__':
    main()
