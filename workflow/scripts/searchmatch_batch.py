#!/usr/bin/env python
"""Run the phase-identification pipeline over a list of samples.

Each sample is searched and refined with `searchmatch.py` machinery; a
failure in one sample is logged and the next sample starts.  A summary CSV
accumulates in results/searchmatch/batch_<name>.csv.

Usage:
  python workflow/scripts/searchmatch_batch.py --name A --samples 510 491 492 ... [--scan 472=panalytical 466=0466_HM_Hf2MoIr_annealed.txt]
         [--top 15] [--top-fixed 5]
"""
import argparse
import sys
import time
import traceback
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))

from searchmatch import RESULTS, refine_top, sample_formula, search, unique_assemblages  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', required=True)
    ap.add_argument('--samples', type=int, nargs='+', required=True)
    ap.add_argument('--scan', nargs='*', default=[], help='sample=scan overrides')
    ap.add_argument('--formula', nargs='*', default=[], help='sample=formula overrides for samples not in the dataframe')
    ap.add_argument('--suffix', nargs='*', default=[], help='sample=suffix: write to results/searchmatch/<sample>_<suffix> (second scan of one sample)')
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--top-fixed', type=int, default=5)
    ap.add_argument('--diverse', action='store_true', help='one assemblage per dominant phase in the refined set')
    ap.add_argument('--no-refine', action='store_true')
    ap.add_argument('--max-excess', type=float, default=0.2, help='eV/atom; hypotheses above this are dropped')
    ap.add_argument('--max-refl', type=int, default=2500, help='skip assemblages with more reflections than this')
    a = ap.parse_args()
    scans = dict(s.split('=', 1) for s in a.scan)
    formulas = dict(s.split('=', 1) for s in a.formula)
    suffixes = dict(s.split('=', 1) for s in a.suffix)
    summary = RESULTS / f'batch_{a.name}.csv'
    rows = []
    for s in a.samples:
        t0 = time.time()
        formula = formulas.get(str(s)) or sample_formula(s)
        scan = scans.get(str(s), 'siemens')
        out_dir = (RESULTS / f'{s:04d}_{suffixes[str(s)]}') if str(s) in suffixes else None
        print(f'\n===== {s:04d} {formula} ({scan}) {time.strftime("%H:%M")}', flush=True)
        rec = dict(sample=s, formula=formula, scan=scan, status='ok')
        try:
            S = search(s, formula, scan, max_excess=a.max_excess, out_dir=out_dir)
            top = unique_assemblages(S['res']).head(3)
            rec.update(search_s=round(sum(S['timing'].values())), n_phases=len(S['phases']),
                       n_scored=S['res'].attrs.get('n_total', len(S['res'])), top_search=' | '.join(top['formulas']))
            print(unique_assemblages(S['res']).drop(columns=['keys', 'fracs']).head(10).to_string(), flush=True)
            if not a.no_refine:
                df = refine_top(S, a.top, top_fixed=a.top_fixed, max_refl=a.max_refl, diverse=a.diverse)
                ok = df[(df['failed'] != True) & (df['mode'] == 'free')]  # noqa: E712
                if not ok.empty:
                    b = ok.sort_values('wR').iloc[0]
                    rec.update(best_free=b['formulas'], best_wR=round(b['wR'], 2), best_gof=round(b['gof'], 2),
                               best_wt=b['wt'], n_refined=len(df), n_failed=int(df['failed'].sum()))
                    print(ok.sort_values('wR').head(8)[['rank', 'formulas', 'wR', 'gof', 'wt', 'cell_change', 'seconds']].to_string(), flush=True)
        except Exception:
            rec['status'] = 'error: ' + traceback.format_exc().strip().splitlines()[-1]
            traceback.print_exc()
        rec['minutes'] = round((time.time() - t0) / 60, 1)
        rows.append(rec)
        pd.DataFrame(rows).to_csv(summary, index=False)
        print(f'----- {s:04d} done in {rec["minutes"]} min', flush=True)


if __name__ == '__main__':
    main()
