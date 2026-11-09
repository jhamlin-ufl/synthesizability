"""Single-phase Rietveld fit of the predicted P4/mmm TiSiPt5 cell to the 0493
annealed scan (1000 C 7 d, zero-background slide, 8 to 120 deg).

Sequence (one block at a time; a block is kept only if it lowers Rwp by
MIN_GAIN and its values are physical, otherwise the previous state is
restored):
  1. scale + 8-term Chebyshev background
  2. cell + sample displacement (Shift); Zero held at the quartz calibration
     (Zero and Shift are nearly collinear and are not refined together)
  3. sample broadening on the fixed quartz instrument profile:
     isotropic microstrain (tan theta), then isotropic size (1/cos theta),
     Lorentz/Gauss mixing held at 1 unless refining it stays inside [0, 1]
  4. March-Dollase along (001), the platelet axis of the layered cell
  5. Pt2 z, then Uiso grouped by element type (Ti = Si, Pt2 = Pt3)
  6. background 12 terms if it pays
Diagnostics, always rolled back and only reported: Ti and Si site
occupancies with Uiso held, and a free Zero.  Kalpha2/Kalpha1 is held at 0.5
(no monochromator).
Output: results/phaseid/review/update_2026-10-22/annealed_single/fit_of_record/

The same sequence serves other Pt5InSi-type scans (2026-10-23: 0561 ZrSiPt5
annealed): --scan, --cif, --out and --name override the 0493 defaults.  The
non-Pt sites (Ti or Zr, and Si) are the ones held at U_LIGHT and grouped.
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from synthesizability.parsers.xrd import parse_xrd_file  # noqa: E402
from synthesizability.phaseid import refine as rf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/phaseid/review/update_2026-10-22/annealed_single/fit_of_record'
SCAN = ROOT / 'samples/TiPt5Si/buttons/0493/0493_HM_TiPt5Si_annealed_ZB.txt'
CIF = ROOT / 'results/phaseid/library/Pt-Si-Ti/TiSiPt5_alex_agm027752216.cif'
LIMITS = (8.0, 119.98)
MIN_GAIN = 0.05     # Rwp points a block must gain to be kept
MAX_CYCLES = 20     # least-squares cycles per GSAS-II call (stops earlier on convergence)


def main():
    global OUT, SCAN, CIF, LIMITS
    ap = argparse.ArgumentParser()
    ap.add_argument('--scan', type=Path, default=SCAN)
    ap.add_argument('--cif', type=Path, default=CIF)
    ap.add_argument('--out', type=Path, default=OUT)
    ap.add_argument('--name', default='TiSiPt5', help='phase name in the project and the tick file')
    ap.add_argument('--hi', type=float, default=LIMITS[1],
                    help='upper 2theta limit: the scan end, or a background point just below it (0561 runs to 150 deg)')
    a = ap.parse_args()
    OUT, SCAN, CIF = a.out.resolve(), a.scan.resolve(), a.cif.resolve()
    LIMITS = (LIMITS[0], a.hi)
    OUT.mkdir(parents=True, exist_ok=True)
    dat = parse_xrd_file(SCAN)
    xy = rf.write_xy(dat['two_theta'], dat['intensity'], OUT / 'scan.xy')
    r = rf.Refinement(xy, rf.instprm_for('siemens'), [str(CIF)], names=[a.name],
                      workdir=OUT / 'work', limits=LIMITS, n_bkg=8)
    # The pipeline caps GSAS-II at 3 least-squares cycles per call (screening
    # speed); a careful fit needs each block converged
    r.gpx.set_Controls('cycles', MAX_CYCLES)
    hname = r.hist.name
    log = []
    ph = lambda: r.phases[0]                       # noqa: E731  (re-read after every rollback)
    hap = lambda: ph().data['Histograms'][hname]   # noqa: E731

    def snapshot():
        res = r.results()
        return dict(cell={k: round(v, 5) for k, v in res['phases'][0]['cell'].items() if k.startswith('length')},
                    shift=r.hist.data['Sample Parameters']['Shift'][0], zero=r.hist.data['Instrument Parameters'][0]['Zero'][1],
                    size=hap()['Size'][1][:3], mustrain=hap()['Mustrain'][1][:3], md=hap()['Pref.Ori.'][1] if hap()['Pref.Ori.'][2] else None,
                    atoms={a.label: dict(z=round(float(a.coordinates[2]), 5), uiso=round(float(a.uiso), 5), occ=round(float(a.occupancy), 3)) for a in ph().atoms()})

    def stage(name, apply, force=False, accept=None, diagnostic=False):
        r.save()
        before = r.hist.get_wR() if log else None
        apply()
        try:
            wR = r._refine()
        except rf.RefinementFailed as exc:
            print(f'{name}: failed ({exc}), rolled back', flush=True)
            r._reload(); log.append(dict(stage=name, wR=None, kept=False, note=str(exc))); return False
        kept = force or before is None or wR < before - MIN_GAIN
        note = '' if kept else f'no gain over {before:.2f}'
        if kept and accept is not None and not accept():
            kept, note = False, 'unphysical value'
        if diagnostic:
            kept, note = False, 'diagnostic'
        snap = snapshot()
        print(f'{name}: Rwp {wR:.2f}' + ('' if kept else f' (rolled back: {note})'), flush=True)
        print(f'    {snap}', flush=True)
        log.append(dict(stage=name, wR=round(wR, 3), kept=bool(kept), note=note, state=snap))
        if not kept:
            r._reload()
        return kept

    def mix(value, refine):
        h = hap()
        h['Size'][1][2] = h['Mustrain'][1][2] = value
        h['Size'][2][2] = h['Mustrain'][2][2] = refine

    stage('scale + background (8 Chebyshev)', lambda: r.hist.set_refinements({'Sample Parameters': ['Scale']}), force=True)
    stage('cell + displacement (Shift)', lambda: (ph().set_refinements({'Cell': True}),
                                                 r.hist.set_refinements({'Sample Parameters': ['Shift']})), force=True)

    def strain():
        mix(1.0, False)
        r._bound(ph(), 'Mustrain;i', 0.0, 60000.0)
        ph().set_HAP_refinements({'Mustrain': {'type': 'isotropic', 'refine': True}})
    stage('isotropic microstrain (Lorentzian)', strain, accept=lambda: 0 < hap()['Mustrain'][1][0] < 60000)

    def size():
        r._bound(ph(), 'Size;i', 0.02, 10.0)
        ph().set_HAP_refinements({'Size': {'type': 'isotropic', 'refine': True}})
    stage('isotropic size (Lorentzian)', size, accept=lambda: 0.02 <= hap()['Size'][1][0] < 10.0)

    def free_mix():
        r._bound(ph(), 'Size;mx', 0.0, 1.0); r._bound(ph(), 'Mustrain;mx', 0.0, 1.0)
        mix(1.0, True)
    stage('Lorentz/Gauss mixing', free_mix,
          accept=lambda: 0 <= hap()['Size'][1][2] <= 1 and 0 <= hap()['Mustrain'][1][2] <= 1)

    def po():
        p = hap()['Pref.Ori.']
        p[0] = 'MD'; p[1] = 1.0; p[2] = True; p[3] = [0, 0, 1]
    stage('March-Dollase (001)', po, accept=lambda: 0.7 <= hap()['Pref.Ori.'][1] <= 1.3)
    stage('Pt2 z', lambda: ph().set_refinements({'Atoms': {'Pt2': 'X'}}))

    def uiso():
        for a in ph().atoms():
            a.uiso = 0.005
        ph().set_refinements({'Atoms': {a.label: 'U' for a in ph().atoms()}})
        r.gpx.add_EquivConstr(['0::AUiso:0', '0::AUiso:1'])   # Ti = Si
        r.gpx.add_EquivConstr(['0::AUiso:2', '0::AUiso:3'])   # Pt2 = Pt3
    stage('Uiso grouped (Ti = Si, Pt2 = Pt3)', uiso, accept=lambda: all(0 < a.uiso < 0.05 for a in ph().atoms()))
    # Anisotropic strain broadening: grinding a ductile fcc-derived
    # intermetallic leaves hkl-dependent strain and faulting.  Uniaxial (c
    # against the basal plane) first; then the Stephens generalized model,
    # started from the isotropic value as the GSAS-II GUI does
    def uniaxial_strain():
        mu = hap()['Mustrain'][1][0]
        hap()['Mustrain'][1][:2] = [mu, mu]
        ph().set_HAP_refinements({'Mustrain': {'type': 'uniaxial', 'direction': [0, 0, 1], 'refine': True}})
    stage('uniaxial microstrain along (001)', uniaxial_strain, accept=lambda: min(hap()['Mustrain'][1][:2]) > 0)

    def generalized_strain():
        from GSASII import GSASIIspc as G2spc
        g = ph().data['General']
        mu = max(hap()['Mustrain'][1][:2]) if hap()['Mustrain'][0] == 'uniaxial' else hap()['Mustrain'][1][0]
        hap()['Mustrain'][4] = list(G2spc.Muiso2Shkl(mu, g['SGData'], g['Cell'][1:7]))
        ph().set_HAP_refinements({'Mustrain': {'type': 'generalized', 'refine': True}})
    stage('generalized (Stephens) microstrain', generalized_strain)

    def size_again():
        hap()['Size'][1][0] = 0.2
        r._bound(ph(), 'Size;i', 0.02, 10.0)
        ph().set_HAP_refinements({'Size': {'type': 'isotropic', 'refine': True}})
    stage('isotropic size, after strain model settled', size_again, accept=lambda: 0.02 <= hap()['Size'][1][0] < 10.0)
    stage('background 12 Chebyshev', lambda: r.hist.set_refinements({'Background': {'type': 'chebyschev-1', 'no. coeffs': 12, 'refine': True}}))
    stage('final cycle, all kept parameters', lambda: None, force=True)

    r.save()
    record = dict(wR=r.hist.get_wR(), gof=r.results().get('gof'), state=snapshot())
    r.export_pattern(OUT / 'fit.txt')
    (OUT / 'ticks.json').write_text(json.dumps({k: [float(v) for v in a] for k, a in r.reflections().items()}))

    def occupancies():
        for a in ph().atoms():
            a.refinement_flags = ''
        ph().set_refinements({'Atoms': {at.label: 'F' for at in ph().atoms() if at.element != 'Pt'}})
    stage('DIAGNOSTIC Ti and Si occupancy (Uiso held)', occupancies, diagnostic=True)
    stage('DIAGNOSTIC free Zero with Shift', lambda: r.hist.set_refinements({'Instrument Parameters': ['Zero']}), diagnostic=True)

    summary = dict(record=record, limits=LIMITS, stages=log)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(record, indent=1, default=float))


if __name__ == '__main__':
    main()
