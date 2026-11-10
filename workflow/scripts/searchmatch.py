#!/usr/bin/env python
"""Automated multiphase XRD identification for one sample.

Stages (each writes into results/searchmatch/<sample>/):

1. library     symmetrised candidate phases of the chemical space from
               OQMD, Alexandria and MP, with disordered twins
2. mixtures    hull-geometry mixture hypotheses (per-database hulls merged)
3. searchmatch NNLS ranking of every hypothesis -> candidates.csv
4. refine      GSAS-II refinement of the top candidates, free and fixed
               lever-rule fractions -> refinements.csv, fits/*.txt

Usage:
  python workflow/scripts/searchmatch.py --sample 510 [--scan siemens|panalytical|<file>]
        [--top 50] [--max-excess 0.15] [--no-refine] [--time-budget SECONDS]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))

from pymatgen.core import Composition  # noqa: E402

from synthesizability.parsers.xrd import parse_xrd_file  # noqa: E402
from synthesizability.phaseid import searchmatch as sm  # noqa: E402
from synthesizability.phaseid.library import LIB_DIR, build_library, fractions  # noqa: E402
from synthesizability.phaseid.mixtures import atom_to_weight_fractions, enumerate_all  # noqa: E402

RESULTS = LIB_DIR.parent


def sample_formula(sample_number):
    """Nominal formula of a button, from tables/buttons.csv."""
    df = pd.read_csv(Path(__file__).resolve().parents[2] / 'tables' / 'buttons.csv', dtype={'sample': str})
    row = df[df['sample'].astype(int) == int(sample_number)]
    if row.empty:
        raise SystemExit(f'sample {sample_number} not in tables/buttons.csv')
    return row.iloc[0]['formula']


def _numpy_errors_default():
    """GSAS-II switches numpy to raise on floating-point underflow, which then
    aborts the next sample's pattern simulation; restore the defaults."""
    np.seterr(divide='warn', over='warn', under='ignore', invalid='warn')


def search(sample_number, formula, scan='siemens', max_excess=None, prior_scale=1.0,
           two_theta_range=sm.DEFAULT_RANGE, out_dir=None, verbose=True, structure_groups=True):
    _numpy_errors_default()
    out_dir = out_dir or RESULTS / f'{sample_number:04d}'
    out_dir.mkdir(parents=True, exist_ok=True)
    timing = {}
    t = time.time()
    elements, phases = build_library(formula, sample_number=sample_number, verbose=verbose)
    timing['library_s'] = time.time() - t
    t = time.time()
    groups = None
    if structure_groups:
        from synthesizability.phaseid.structure_types import framework_groups
        groups = framework_groups(phases, '-'.join(elements), verbose=verbose)
    timing['groups_s'] = time.time() - t

    t = time.time()
    c = fractions(Composition(formula), elements)
    X, members, mix = enumerate_all(phases, c, max_excess=max_excess)
    timing['mixtures_s'] = time.time() - t
    if verbose:
        print(f'{len(X)} composition points, {len(mix)} mixture hypotheses '
              f'(hulls: { {k: round(v * 1000, 1) for k, v in mix.attrs["hulls"].items()} } meV/atom)')

    path = sm.choose_scan(sample_number, scan)
    grid, y = sm.load_observed(path, two_theta_range)
    t = time.time()
    P, ptable = sm.phase_patterns(phases, grid, y)
    timing['phase_patterns_s'] = time.time() - t
    t = time.time()
    res = sm.score_mixtures(mix, members, phases, P, y, prior_scale=prior_scale, phase_table=ptable)
    timing['score_s'] = time.time() - t
    n_total = res.attrs.get('n_total', len(res))
    if verbose:
        print(f'{n_total} phase combinations scored ({len(res)} kept); '
              f'timing {json.dumps({k: round(v, 1) for k, v in timing.items()})}')

    ptable.to_csv(out_dir / 'phases.csv', index=False)
    res.drop(columns=['keys', 'fracs']).to_csv(out_dir / 'candidates.csv', index=False)
    meta = dict(sample=sample_number, formula=formula, scan=path.name, elements=elements,
                n_phases=len(phases), n_points=len(X), n_mixtures=len(mix), n_scored=int(n_total),
                two_theta_range=list(two_theta_range), timing=timing)
    (out_dir / 'search.json').write_text(json.dumps(meta, indent=1, default=str))
    return dict(elements=elements, phases=phases, members=members, mix=mix, grid=grid, y=y,
                P=P, ptable=ptable, res=res, scan=path, out_dir=out_dir, timing=timing, groups=groups)


MAX_REFL = 2500   # combinations with more reflections than this are too slow in GSAS-II


def unique_assemblages(res, groups=None):
    """Best-scoring row per assemblage.

    Without ``groups`` the key is the set of
    formulas: polymorphs of the same formula collapse onto one candidate, but
    a fully disordered twin counts as a different assemblage from its ordered
    parent (marked '*'), so both reach GSAS-II and the ordered-versus-
    disordered question is settled by refinement rather than by the quick
    ranking.

    With ``groups`` (``structure_types.framework_groups`` table)
    the key is formula plus framework group plus a disorder mark: the same
    structure from OQMD, MP and Alexandria collapses onto its best-scoring
    copy, while polymorphs of one formula stay apart.  The collapse is applied
    to the whole kept list before any cut, so no hypothesis is lost to a
    crowd of its own near-copies; the replaced rows are counted in
    ``variants``."""
    if groups is None:
        def key(row):
            parts = []
            for f, k in zip(row['formulas'].split(' + '), row['phases'].split(' + ')):
                parts.append(f + ('*' if (k.endswith('_dis') or '_proto_' in k) else ''))
            return ' + '.join(sorted(parts))
    else:
        g = groups.set_index('key') if 'key' in groups.columns else groups
        gid, dis = g['group'].to_dict(), (g['role'] == 'disordered member').to_dict()

        def key(row):
            parts = []
            for f, k in zip(row['formulas'].split(' + '), row['phases'].split(' + ')):
                parts.append(f'{f}@{gid[k]}' + ('*' if dis[k] else '') if k in gid else k)
            return ' + '.join(sorted(parts))
    out = res.assign(assemblage=res.apply(key, axis=1))
    out['variants'] = out.groupby('assemblage')['assemblage'].transform('size')
    return out.drop_duplicates('assemblage', keep='first')


def dominant_formula(row):
    """Formula of the phase carrying the largest NNLS share of the fit."""
    keys = row['phases'].split(' + ')
    shares = [float(x) for x in str(row['fit_share']).split(' / ')]
    return keys[int(np.argmax(shares))].split('_')[0]


REFINE_TIMEOUT = 1200   # s per refinement; the 20 to 50 min tail of the overnight runs was dense cells of poor fit


def refine_top(S, top=50, time_budget=None, instrument=None, verbose=True, max_refl=MAX_REFL,
               top_fixed=None, max_forced=6, diverse=False, cell_scan=True, pair_hypotheses=True,
               timeout=REFINE_TIMEOUT):
    """GSAS-II refinement of the ``top`` distinct assemblages.

    Free-fraction refinements run for the top ``top`` assemblages (one row per
    set of formulas); fixed lever-rule refinements run for the ``top_fixed``
    best of those by free wR (all of them when None).  Candidates whose phases
    together have more than ``max_refl`` reflections in range are skipped (a
    70-atom triclinic cell takes minutes per cycle).

    ``diverse=True`` keeps one assemblage per dominant phase when cutting at
    ``top`` (in a known-answer recall test, 0.78 instead of 0.68
    at top 15), and the best single disordered parent at the target
    composition is always included, because the ranking undervalues
    disordered cells (B2 answers at ranks 2 to 55, fcc at 2591).
    """
    from synthesizability.phaseid import refine as rf
    _numpy_errors_default()
    scan = S['scan']
    instrument = instrument or ('panalytical' if scan.suffix == '.xy' else 'siemens')
    instprm = rf.instprm_for(instrument)
    if not Path(instprm).exists():
        raise SystemExit(f'no instrument file {instprm}; run the calibration first')
    out_dir = S['out_dir']
    fits = out_dir / 'fits'
    fits.mkdir(exist_ok=True)
    d = parse_xrd_file(scan)
    xy = rf.write_xy(d['two_theta'], d['intensity'], out_dir / f'{scan.stem}.xy')
    phases = S['phases']
    rows = []
    t_start = time.time()
    res = unique_assemblages(S['res'], S.get('groups'))
    eligible = res[res['n_refl'] <= max_refl] if 'n_refl' in res else res
    skipped = len(res.head(top)) - len(eligible.head(top))
    if verbose and skipped:
        print(f'{skipped} of the top {top} assemblages skipped as too dense (> {max_refl} reflections)')
    if diverse and 'fit_share' in eligible:
        chosen = eligible.assign(_dom=eligible.apply(dominant_formula, axis=1)).drop_duplicates('_dom').head(top).drop(columns='_dom')
        # One assemblage per dominant phase starves samples with a single
        # dominant phase (0507: two hypotheses refined); fill up to ``top``
        # with the next-best assemblages in score order.
        if len(chosen) < top:
            rest = eligible.loc[~eligible.index.isin(chosen.index)].head(top - len(chosen))
            chosen = pd.concat([chosen, rest])
    else:
        chosen = eligible.head(top)
    # Hypotheses the ranking is known to undervalue are refined regardless of
    # rank: the best assemblage containing each disordered parent from
    # relaxed disordered cells ('_p2_') and the best single disordered parent of the target
    forced = []
    seen = set(chosen.index)
    single = eligible[(eligible['n_phases'] == 1)
                      & eligible['phases'].map(lambda k: k.endswith('_dis') or '_p2_' in k or '_proto_' in k)]
    if not single.empty and single.index[0] not in seen:
        forced.append(single.index[0])
        seen.add(single.index[0])
    for marker in ('_p2_',):
        for key in sorted({k for s in eligible['phases'] for k in s.split(' + ') if marker in k}):
            hit = eligible[eligible['phases'].map(lambda s: key in s.split(' + '))]
            if not hit.empty and hit.index[0] not in seen:
                forced.append(hit.index[0])
                seen.add(hit.index[0])
            if len(forced) >= max_forced:
                break
    if forced:
        chosen = pd.concat([chosen, eligible.loc[forced]])
        if verbose:
            print(f'{len(forced)} assemblages with relaxed disordered parents added to the refinement set (ranks {sorted(forced)})')
    key_index = {ph.key: i for i, ph in enumerate(phases)}

    def synthetic(rank, keys, note):
        """A refinement row for a hypothesis the ranking cannot express."""
        n = len(keys)
        return pd.Series(dict(phases=' + '.join(phases[k].key for k in keys), formulas=' + '.join(phases[k].formula for k in keys),
                              assemblage=note, n_phases=n, excess_meV=np.nan, rwp=np.nan, rwp_fixed=np.nan,
                              lever_wt=' / '.join(['%.2f' % (1 / n)] * n), fit_share='', keys=tuple(keys),
                              fracs=tuple([1 / n] * n), n_refl=sum(phases[k].d_hkl.size for k in keys)), name=rank)

    # Hypotheses forced in at the refinement stage, where the
    # strain model can judge them: two cells of the identified solution
    # framework (a cored solution), and the best ordered phase at the target
    # beside its own disordered parent (partial order).
    extra = []
    if pair_hypotheses:
        if not single.empty:
            k = key_index[single.iloc[0]['phases']]
            extra.append(synthetic(90001, [k, k], 'two cells of one solution framework'))
        for _, row in chosen.iterrows():
            ks = list(row['keys'])
            ordered = [k for k in ks if phases[k].at_target and not phases[k].disordered and phases[k].key + '_dis' in key_index]
            if ordered:
                k = ordered[0]
                partner = [q for q in ks if q != k][:1]
                extra.append(synthetic(90002, [k, key_index[phases[k].key + '_dis']] + partner, 'ordered phase + its disordered parent'))
                break
    if extra:
        chosen = pd.concat([chosen, pd.DataFrame(extra)])
        if verbose:
            print(f'{len(extra)} pair hypotheses forced into the refinement set: ' + '; '.join(e["assemblage"] for e in extra))

    # Starting cells: per-axis scan of every phase in the refinement set on
    # the stick pattern, so GSAS-II starts within a fraction of a
    # percent of the observed cell instead of at the DFT cell.
    cells_dir = fits / 'cells'
    cells_dir.mkdir(exist_ok=True)
    scanned = {}
    fixed_cells = set()
    if cell_scan:
        fwhm_ab = S['ptable'].attrs.get('fwhm_ab') if hasattr(S.get('ptable'), 'attrs') else None
        if fwhm_ab is None:
            fwhm_ab = sm.estimate_fwhm(S['grid'], S['y'])
        needed = sorted({k for _, row in chosen.iterrows() for k in row['keys']})
        t_scan = time.time()
        for k in needed:
            ph = phases[k]
            if ph.source in ('pdf', 'cod'):
                continue          # experimental cells are measured; never move them (0507: a 1.5 % move cost 7 in wR)
            try:
                s2, sc, r = sm.cell_scan(ph.cif, S['grid'], S['y'], fwhm_ab, ph.spg)
                from pymatgen.io.cif import CifWriter
                path = cells_dir / f'{ph.key}.cif'
                CifWriter(s2, symprec=0.01 if ph.spg > 2 else None).write_file(str(path))
                scanned[k] = (path, sc)
            except Exception as exc:
                if verbose:
                    print(f'  cell scan failed for {ph.key}: {exc}')
            if ph.spg <= 2:
                fixed_cells.add(ph.key[:40])
        if verbose:
            moved = [f"{phases[k].key}: {'/'.join('%+.1f%%' % (100 * (v - 1)) for v in sc)}" for k, (_, sc) in scanned.items() if any(abs(v - 1) > 0.004 for v in sc)]
            print(f'cell scan of {len(needed)} phases in {time.time() - t_scan:.0f} s; moved: ' + ('; '.join(moved) if moved else 'none'))
    todo = [(rank, row, 'free') for rank, row in chosen.iterrows()]
    free_results = {}

    def run_one(rank, row, mode):
        combo = row['keys']
        cifs = [scanned[k][0] if k in scanned else phases[k].cif for k in combo]
        names = [phases[k].key[:40] for k in combo]
        # the same phase twice (two cells of one framework): distinct names and
        # starting cells 0.7 % apart so the least squares can separate them
        if len(set(names)) < len(names):
            from pymatgen.core import Structure
            from pymatgen.io.cif import CifWriter
            new_cifs, new_names, seen_n = [], [], {}
            for c, nm in zip(cifs, names):
                seen_n[nm] = seen_n.get(nm, 0) + 1
                if seen_n[nm] > 1:
                    s = Structure.from_file(str(c))
                    s.scale_lattice(s.volume * 1.021)   # a 0.7 % larger cell for the second copy
                    c2 = cells_dir / f'{nm[:34]}_copy{seen_n[nm]}.cif'
                    CifWriter(s).write_file(str(c2))
                    new_cifs.append(c2); new_names.append(f'{nm[:34]}_copy{seen_n[nm]}')
                else:
                    new_cifs.append(c); new_names.append(nm)
            cifs, names = new_cifs, new_names
        comps = [phases[k].composition for k in combo]
        tag = f'{rank:04d}_{mode}'
        t = time.time()
        try:
            r = rf.Refinement(xy, instprm, cifs, names=names, workdir=fits / 'work' / tag,
                              limits=tuple(S['grid'][[0, -1]]))
            if mode == 'fixed':
                r.set_fixed_fractions(atom_to_weight_fractions(row['fracs'], comps))
            # Solution cells (disordered parents, prototype and relaxed random
            # cells) carry estimated lattice parameters, so they may move 5 %;
            # ordered database and experimental cells 2 % (per phase)
            with rf.time_limit(timeout):
                out = r.run(free_fractions=(mode == 'free'), max_cell_change=rf.cell_windows([phases[k].key for k in combo]),
                            fixed_cells=fixed_cells)
            r.export_pattern(fits / f'{tag}.txt')
            ticks = {k: v.tolist() for k, v in r.reflections().items()}
            (fits / f'{tag}_ticks.json').write_text(json.dumps(ticks))
            r.cleanup()
        except Exception as exc:  # keep going; record the failure
            out = dict(failed=True, error=str(exc), wR=None, gof=None, phases=[], seconds=time.time() - t)
        rec = dict(rank=rank, mode=mode, phases=row['phases'], formulas=row['formulas'],
                   assemblage=row['assemblage'],
                   n_phases=row['n_phases'], excess_meV=row['excess_meV'], nnls_rwp=row['rwp'],
                   nnls_rwp_fixed=row.get('rwp_fixed'), lever_wt=row['lever_wt'],
                   wR=out.get('wR'), Rp=out.get('Rp'), gof=out.get('gof'), n_var=out.get('n_var'),
                   failed=out.get('failed'), error=out.get('error'),
                   wt=' / '.join(f"{p['wt']:.2f}" for p in out.get('phases', [])),
                   rf=' / '.join(f"{p['rf']:.1f}" if p['rf'] is not None else 'nan' for p in out.get('phases', [])),
                   cell_change=' / '.join(
                       '{:+.2f}'.format(100 * max(p['cell_change'].values(), key=abs)) if p['cell_change'] else ''
                       for p in out.get('phases', [])),
                   seconds=out.get('seconds'), stage_s=' / '.join(str(x) for x in out.get('stage_s', [])),
                   start_scale=' | '.join('/'.join('%.3f' % v for v in scanned[k][1]) if k in scanned else '' for k in combo))
        rows.append(rec)
        if verbose:
            print(f"  {tag}: {row['formulas']}: wR={out.get('wR')} gof={out.get('gof')} "
                  f"wt={rec['wt']} {'FAILED ' + str(out.get('error')) if out.get('failed') else ''} "
                  f"({out.get('seconds', 0):.1f} s)", flush=True)
        pd.DataFrame(rows).to_csv(out_dir / 'refinements.csv', index=False)
        return rec

    for rank, row, mode in todo:
        if time_budget and time.time() - t_start > time_budget:
            if verbose:
                print(f'time budget reached before rank {rank}')
            break
        rec = run_one(rank, row, mode)
        if not rec['failed'] and rec['wR'] is not None:
            free_results[rank] = (rec['wR'], row)
    # Fixed lever-rule fractions for the best free fits with more than one phase
    best = sorted(free_results.items(), key=lambda kv: kv[1][0])
    n_fixed = len(best) if top_fixed is None else top_fixed
    for rank, (_, row) in best[:n_fixed]:
        if len(row['keys']) == 1 or rank >= 90000:
            continue
        if time_budget and time.time() - t_start > time_budget:
            break
        run_one(rank, row, 'fixed')
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample', type=int, required=True)
    ap.add_argument('--formula')
    ap.add_argument('--scan', default='siemens')
    ap.add_argument('--top', type=int, default=50, help='distinct assemblages refined with free fractions')
    ap.add_argument('--top-fixed', type=int, default=None, help='best free fits also refined with fixed lever-rule fractions')
    ap.add_argument('--max-excess', type=float, default=None, help='eV/atom; drop hypotheses above this')
    ap.add_argument('--prior-scale', type=float, default=1.0)
    ap.add_argument('--diverse', action='store_true', help='one assemblage per dominant phase in the refined set')
    ap.add_argument('--no-refine', action='store_true')
    ap.add_argument('--time-budget', type=float, default=None, help='seconds for the refinement stage')
    ap.add_argument('--range', type=float, nargs=2, default=sm.DEFAULT_RANGE)
    a = ap.parse_args()
    formula = a.formula or sample_formula(a.sample)
    print(f'{a.sample:04d} {formula}')
    S = search(a.sample, formula, a.scan, a.max_excess, a.prior_scale, tuple(a.range))
    print(unique_assemblages(S['res']).drop(columns=['keys', 'fracs']).head(15).to_string())
    if not a.no_refine:
        df = refine_top(S, a.top, a.time_budget, top_fixed=a.top_fixed, diverse=a.diverse)
        ok = df[df['failed'] != True]  # noqa: E712
        if not ok.empty:
            print(ok.sort_values('wR').head(10).to_string())


if __name__ == '__main__':
    main()
