"""0493 annealed scan: the single-phase TiSiPt5 main fit plus one
candidate fcc-type phase for the shifted humps on the fcc fundamental
reflections (residual maxima give a = 3.89-3.90 A).

Each candidate starts from the converged single-phase project
(annealed_single/fit_of_record/work/fit.gpx).  Sequence, a block kept only if
it lowers Rwp by MIN_GAIN and stays physical:
  1. candidate scale, TiSiPt5 parameters held
  2. one broadening term for the candidate: isotropic microstrain, or
     isotropic size if that fits better (Lorentzian, mixing held at 1).  A
     single broad hump cannot carry both; refined together they trade off
     into negative strain and nanometre sizes
  3. candidate cell, after the broadening (a sharp calculated peak cannot
     follow a broad hump, so the cell does not move before it)
  4. joint cycle: every kept TiSiPt5 parameter with the candidate scale and
     cell; the candidate broadening is held at its stage-2 value
Candidate Uiso stays at its CIF value (0.005 if none).
Output: results/phaseid/review/update_2026-10-22/annealed_two_phase/<name>/

Usage: python workflow/scripts/refine_0493_annealed_two_phase.py [--sample 493|561] [name ...]
"""
import json, shutil, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from synthesizability.phaseid import refine as rf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / 'results/phaseid/library/Pt-Si-Ti'
BASE = ROOT / 'results/phaseid/review/update_2026-10-22/annealed_single/fit_of_record/work/fit.gpx'
OUT = ROOT / 'results/phaseid/review/update_2026-10-22/annealed_two_phase'
CANDIDATES = {
    'A1_Pt3Ti': (LIB / 'TiPt3_cod_1539263_dis.cif', 'disordered A1 Pt0.75Ti0.25 (from TiPt3 COD 1539263)'),
    'L12_TiPt3': (LIB / 'TiPt3_cod_1539263.cif', 'ordered L1_2 TiPt3 (COD 1539263)'),
    'A1_Pt': (OUT / 'Pt_A1_a3.90.cif', 'fcc Pt (COD 1011103 atoms), start a = 3.90: Pt-rich solution'),
    'A1_TiSiPt5': (OUT / 'TiSiPt5_A1_a3.90.cif', 'disordered A1 at the TiSiPt5 composition, start a = 3.90'),
    # The fit-rule congruency scan names liquid + Ti3Pt7 as the peritectic
    # competitor; DFT-only C2/m cell (Alexandria agm003291108), lengths scaled
    # by 0.986 to the fcc-equivalent a = 3.90 A of the humps
    'Ti3Pt7': (OUT / 'Ti3Pt7_scaled_a3.90.cif', 'ordered C2/m Ti3Pt7 (Alexandria agm003291108), scaled to a_fcc 3.90'),
    'Ti3Pt7_dis': (OUT / 'Ti3Pt7_dis_scaled_a3.90.cif', 'Ti3Pt7 C2/m with mixed sites (agm003291108_dis), scaled to a_fcc 3.90'),
}
# 2026-10-23: the same test for the 0561 ZrSiPt5 annealed scan; the as-cast
# partners (D0_24 ZrPt3, Pt2Si, Pt-rich fcc ZrPt8) are the candidates
LIB_ZR = ROOT / 'results/phaseid/library/Pt-Si-Zr'
SAMPLES = {
    '493': (BASE, OUT, CANDIDATES),
    '561': (ROOT / 'results/phaseid/review/update_2026-10-23/0561_annealed_single/work/fit.gpx',
            ROOT / 'results/phaseid/review/update_2026-10-23/0561_annealed_two_phase',
            {'ZrPt3': (LIB_ZR / 'ZrPt3_oqmd_2017725.cif', 'D0_24 ZrPt3 (OQMD 2017725, the as-cast majority phase)'),
             'Pt2Si': (LIB_ZR / 'SiPt2_pdf_04-003-1631.cif', 'Pt2Si I4/mmm (PDF 04-003-1631, as-cast partner)'),
             'ZrPt8': (LIB_ZR / 'ZrPt8_oqmd_1374867.cif', 'Pt-rich fcc-derived ZrPt8 (OQMD 1374867, as-cast partner)')}),
}
MIN_GAIN = 0.05
MAX_CYCLES = 20


class Project(rf.Refinement):
    """A Refinement opened from an existing gpx instead of built from a scan."""

    def __init__(self, gpx_src, workdir, limits):
        G2sc = rf.g2sc()
        self.workdir = Path(workdir); self.workdir.mkdir(parents=True, exist_ok=True)
        self.gpx_path = self.workdir / 'fit.gpx'
        shutil.copy(gpx_src, self.gpx_path)
        with rf.quiet():
            self.gpx = G2sc.G2Project(str(self.gpx_path))
        self.hist = self.gpx.histograms()[0]
        self.phases = self.gpx.phases()
        self.limits = limits
        self.log, self.stage_s = [], []


def flags_main(r, on, kept):
    """Switch the kept TiSiPt5 parameters (cell, Pt2 z, size, strain, MD,
    Shift) on or off."""
    ph = r.phases[0]; hap = ph.data['Histograms'][r.hist.name]
    ph.data['General']['Cell'][0] = on
    for a in ph.atoms():
        f = ('X' if kept.get('pt2z', True) and a.label == 'Pt2' else '') + ('U' if kept.get('uiso', False) else '')
        a.refinement_flags = f if on else ''
    hap['Size'][2][0] = on and kept['size']
    hap['Mustrain'][2][0] = on
    if hap['Mustrain'][0] == 'generalized':
        hap['Mustrain'][5] = [on] * len(hap['Mustrain'][5])
    hap['Pref.Ori.'][2] = on
    r.hist.data['Sample Parameters']['Shift'][1] = on


def run(name, sample='493'):
    base, out_root, candidates = SAMPLES[sample]
    cif, title = candidates[name]
    base_summary = json.loads((base.parents[1] / 'summary.json').read_text())
    was_kept = lambda prefix: any(s['kept'] and s['stage'].startswith(prefix) for s in base_summary['stages'])  # noqa: E731
    kept = dict(size=any(s['kept'] and 'size' in s['stage'] for s in base_summary['stages']),
                pt2z=was_kept('Pt2 z'), uiso=was_kept('Uiso'))
    out_dir = out_root / name
    r = Project(base, out_dir / 'work', tuple(base_summary['limits']))
    r.gpx.set_Controls('cycles', MAX_CYCLES)
    with rf.quiet():
        r.gpx.add_phase(str(cif), phasename=name, histograms=[r.hist])
    r.phases = r.gpx.phases()
    hname = r.hist.name
    log = []
    cand = lambda: r.phases[1]                      # noqa: E731
    hap = lambda: cand().data['Histograms'][hname]  # noqa: E731
    for a in cand().atoms():
        if not a.uiso or a.uiso <= 0:
            a.uiso = 0.005

    def state():
        res = r.results()
        return dict(wt=[round(p['wt'], 4) for p in res['phases']],
                    cells=[{k: round(v, 5) for k, v in p['cell'].items() if k in ('length_a', 'length_c')} for p in res['phases']],
                    cand_size=hap()['Size'][1][0], cand_mustrain=hap()['Mustrain'][1][0],
                    main_size=r.phases[0].data['Histograms'][hname]['Size'][1][0],
                    main_mustrain=(r.phases[0].data['Histograms'][hname]['Mustrain'][4]
                                   if r.phases[0].data['Histograms'][hname]['Mustrain'][0] == 'generalized'
                                   else r.phases[0].data['Histograms'][hname]['Mustrain'][1][0]),
                    shift=r.hist.data['Sample Parameters']['Shift'][0])

    def stage(label, apply, force=False, accept=None):
        r.save()
        before = r.hist.get_wR() if log else None
        apply()
        try:
            wR = r._refine()
        except rf.RefinementFailed as exc:
            print(f'  {label}: failed ({exc}), rolled back', flush=True)
            r._reload(); r.phases = r.gpx.phases(); log.append(dict(stage=label, wR=None, kept=False, note=str(exc))); return False
        ok = force or before is None or wR < before - MIN_GAIN
        note = '' if ok else f'no gain over {before:.2f}'
        if ok and accept is not None and not accept():
            ok, note = False, 'unphysical value'
        st = state()
        print(f'  {label}: Rwp {wR:.2f}' + ('' if ok else f' (rolled back: {note})') + f'   {st}', flush=True)
        log.append(dict(stage=label, wR=round(wR, 3), kept=bool(ok), note=note, state=st))
        if not ok:
            r._reload(); r.phases = r.gpx.phases()
        return ok

    print(f'== {name}: {title}', flush=True)
    flags_main(r, False, kept)
    r.hist.set_refinements({'Sample Parameters': ['Scale']})

    def scale():
        cand().set_HAP_refinements({'Scale': True})
        hap()['Size'][1][2] = hap()['Mustrain'][1][2] = 1.0
        hap()['Size'][2][2] = hap()['Mustrain'][2][2] = False
    stage('candidate scale (TiSiPt5 held)', scale, force=True)

    def strain():
        r._bound(cand(), 'Mustrain;i', 0.0, 60000.0)
        cand().set_HAP_refinements({'Mustrain': {'type': 'isotropic', 'refine': True}})
    strain_ok = stage('candidate isotropic microstrain', strain, accept=lambda: 0 < hap()['Mustrain'][1][0] < 60000)
    w_strain = log[-1]['wR']

    def size_instead():
        hap()['Mustrain'][1][0] = 1000.0
        cand().clear_HAP_refinements({'Mustrain': {'type': 'isotropic', 'refine': True}})
        hap()['Mustrain'][2][0] = False
        hap()['Size'][1][0] = 0.1
        r._bound(cand(), 'Size;i', 0.005, 10.0)
        cand().set_HAP_refinements({'Size': {'type': 'isotropic', 'refine': True}})
    # Size in place of strain: kept only if it beats the strain stage
    r.save()
    before = r.hist.get_wR()
    if strain_ok and w_strain is not None:
        size_instead()
        try:
            w_size = r._refine()
        except rf.RefinementFailed:
            w_size = None
        good = w_size is not None and w_size < before - MIN_GAIN and 0.005 <= hap()['Size'][1][0] < 10.0
        print(f"  candidate isotropic size instead of strain: Rwp {w_size}" + ('' if good else ' (rolled back: strain kept)'), flush=True)
        log.append(dict(stage='candidate isotropic size instead of strain', wR=None if w_size is None else round(w_size, 3), kept=bool(good),
                        note='' if good else 'strain fits as well or better'))
        if not good:
            r._reload(); r.phases = r.gpx.phases()
    cell0 = cand().get_cell()
    stage('candidate cell', lambda: cand().set_refinements({'Cell': True}),
          accept=lambda: all(abs(cand().get_cell()[k] / cell0[k] - 1) < 0.03 for k in ('length_a', 'length_b', 'length_c')))
    def joint():
        flags_main(r, True, kept)
        hap()['Size'][2][0] = False; hap()['Mustrain'][2][0] = False
    stage('joint cycle (TiSiPt5 + candidate scale and cell)', joint, force=True)

    r.save()
    res = r.results()
    record = dict(name=name, title=title, wR=r.hist.get_wR(), gof=res.get('gof'), state=state(),
                  main_atoms={a.label: round(float(a.coordinates[2]), 5) for a in r.phases[0].atoms()},
                  md=r.phases[0].data['Histograms'][hname]['Pref.Ori.'][1])
    r.export_pattern(out_dir / 'fit.txt')
    (out_dir / 'ticks.json').write_text(json.dumps({k: [float(v) for v in a] for k, a in r.reflections().items()}))
    (out_dir / 'summary.json').write_text(json.dumps(dict(record=record, stages=log), indent=1, default=float))
    print(f"  RESULT {name}: Rwp {record['wR']:.2f}  {record['state']}", flush=True)


def main():
    args = sys.argv[1:]
    sample = '493'
    if args[:1] == ['--sample']:
        sample, args = args[1], args[2:]
    for name in (args or list(SAMPLES[sample][2])):
        run(name, sample)


if __name__ == '__main__':
    main()
