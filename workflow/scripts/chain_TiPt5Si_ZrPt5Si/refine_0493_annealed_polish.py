"""Residual-shape driven refinements of the 0493 annealed two-phase fit
(TiSiPt5 + disordered A1 Pt3Ti), one block at a time, each kept only if it
lowers Rwp by MIN_GAIN and stays physical:
  1. SH/L (axial-divergence asymmetry): the 41.95 deg (112) line has a
     low-angle tail and the zero-background slide does not sit where the
     quartz standard did
  2. sample transparency (one-sided low-angle tail growing as sin 2theta)
  3. Pt Uiso (Ti and Si held): intensity under-calculated more and more
     with angle (53.8, 78.1, 114.2 deg)
  4. surface roughness (Suortti): sprinkled powder on a thin film is rough
Output: results/refits/chain/0493_polish/ (0561: --sample 561,
results/refits/chain/0561_polish/)
"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from refine_0493_annealed_two_phase import Project, flags_main  # noqa: E402
from synthesizability.phaseid import refine as rf  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
CHAIN = ROOT / 'results/refits/chain'
# sample: (two-phase starting fit, output folder, single-phase summary that holds the limits)
SAMPLES = {
    '493': (CHAIN / '0493_two_phase/A1_Pt3Ti', CHAIN / '0493_polish',
            CHAIN / '0493_single/summary.json'),
    # 0561 ZrSiPt5 annealed, ZrSiPt5 + D0_24 ZrPt3, 8-150 deg; the residual has the
    # Uiso signature (high-angle lines under-, mid-angle lines over-calculated)
    '561': (CHAIN / '0561_two_phase/ZrPt3', CHAIN / '0561_polish',
            CHAIN / '0561_single/summary.json'),
}
MIN_GAIN = 0.05


def main():
    sample = sys.argv[sys.argv.index('--sample') + 1] if '--sample' in sys.argv else '493'
    BASE, OUT, single = SAMPLES[sample]
    base = json.loads((BASE / 'summary.json').read_text())
    limits = tuple(json.loads(single.read_text())['limits'])
    r = Project(BASE / 'work' / 'fit.gpx', OUT / 'work', limits)
    hname = r.hist.name
    log = []
    inst = lambda k: r.hist.data['Instrument Parameters'][0][k][1]   # noqa: E731
    samp = lambda k: r.hist.data['Sample Parameters'][k][0]           # noqa: E731

    def state():
        res = r.results()
        return dict(wt=[round(p['wt'], 4) for p in res['phases']], SHL=inst('SH/L'), transparency=samp('Transparency'),
                    rough=(samp('SurfRoughA'), samp('SurfRoughB')), shift=samp('Shift'),
                    uiso={a.label: round(float(a.uiso), 5) for a in r.phases[0].atoms()},
                    cell=[{k: round(v, 5) for k, v in p['cell'].items() if k in ('length_a', 'length_c')} for p in res['phases']])

    def stage(label, apply, accept=None, force=False):
        r.save()
        before = r.hist.get_wR()
        apply()
        try:
            wR = r._refine()
        except rf.RefinementFailed as exc:
            print(f'{label}: failed ({exc}), rolled back', flush=True)
            r._reload(); r.phases = r.gpx.phases(); log.append(dict(stage=label, wR=None, kept=False, note=str(exc))); return
        ok = force or wR < before - MIN_GAIN
        note = '' if ok else f'no gain over {before:.2f}'
        if ok and accept is not None and not accept():
            ok, note = False, 'unphysical value'
        st = state()
        print(f'{label}: Rwp {wR:.2f}' + ('' if ok else f' (rolled back: {note})') + f'   {st}', flush=True)
        log.append(dict(stage=label, wR=round(wR, 3), kept=bool(ok), note=note, state=st))
        if not ok:
            r._reload(); r.phases = r.gpx.phases()

    # The joint-cycle flags are still set in the saved project; hold everything
    # but scales and background while each new block is tried
    flags_main(r, False, dict(size=True))
    r.phases[1].clear_refinements({'Cell': True})
    stage('starting point (scales + background)', lambda: None, force=True)
    stage('SH/L (axial-divergence asymmetry)', lambda: r.hist.set_refinements({'Instrument Parameters': ['SH/L']}),
          accept=lambda: 0.0 < inst('SH/L') < 0.1)
    stage('sample transparency', lambda: r.hist.set_refinements({'Sample Parameters': ['Transparency']}),
          accept=lambda: samp('Transparency') >= 0)

    def pt_uiso():
        r.phases[0].set_refinements({'Atoms': {a.label: 'U' for a in r.phases[0].atoms() if a.element == 'Pt'}})
        r.gpx.add_EquivConstr(['0::AUiso:2', '0::AUiso:3'])   # Pt2 = Pt3
    stage('Pt Uiso (Pt2 = Pt3; Ti or Zr, and Si held 0.01)', pt_uiso, accept=lambda: all(a.uiso > 0.0005 for a in r.phases[0].atoms()))
    stage('surface roughness (Suortti)', lambda: r.hist.set_refinements({'Sample Parameters': ['SurfRoughA', 'SurfRoughB']}),
          accept=lambda: 0 <= samp('SurfRoughA') <= 1 and 0 <= samp('SurfRoughB') <= 1)

    def joint():
        flags_main(r, True, dict(size=True))
        r.phases[0].data['Histograms'][hname]['Size'][2][0] = True
    stage('joint cycle', joint, force=True)
    r.save()
    res = r.results()
    record = dict(wR=r.hist.get_wR(), gof=res.get('gof'), state=state(), md=r.phases[0].data['Histograms'][hname]['Pref.Ori.'][1],
                  start_wR=base['record']['wR'])
    r.export_pattern(OUT / 'fit.txt')
    (OUT / 'ticks.json').write_text(json.dumps({k: [float(v) for v in a] for k, a in r.reflections().items()}))
    (OUT / 'summary.json').write_text(json.dumps(dict(record=record, stages=log), indent=1, default=float))
    print(json.dumps(record, indent=1, default=float))


if __name__ == '__main__':
    main()
