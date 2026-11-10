"""0615 Al3Cu2Y annealed scans: is the CaCu5-type main phase the predicted
ordered YAl3Cu2 (Al on 3g, Cu on 2c), or an Al/Cu-random Y(Cu,Al)5?

Al/Cu order in the P6/mmm cell adds no reflections; it moves intensity
between lines that both models index: (001) at 21.7 deg is strong when
ordered and weak when random, (101) at 29.1 and (200) at 38.6 the other way.
Texture along c also raises (001), but raises (002) with it, whereas order
lowers (002).  The variants therefore cross the occupancy model with a
March-Dollase (001) term:
  --occ pred     3g Al 1, 2c Cu 1 (held)
  --occ random   both sites Al 0.6 / Cu 0.4 (held, the YAl3Cu2 composition)
  --occ free     3g Al and 2c Cu fractions refined, each site summing to 1
  --md           March-Dollase (001) on the main phase, ratio free
Partners: the two phases of the 0540 fit (rank 90313), the
Y3Al2Cu7_alex_agm005614889_dis and Y3Al8Cu3_alex_agm005688824 cells.

Sequence (one block at a time; a block is kept only if it lowers Rwp by
MIN_GAIN and its values are physical, otherwise the previous state is
restored):
  1. scale + 8-term Chebyshev background
  2. partner phase fractions
  3. cells (main 2 %, ordered partner 2 %, _dis partner 5 %) + Shift; Zero
     held at the quartz calibration
  4. main isotropic microstrain, then size; partner microstrain
  5. March-Dollase (001) on the main phase (--md only)
  6. site fractions (--occ free only)
  7. one Uiso for the main phase
  8. background 12 terms if it pays, final cycle
Output: <out>/fit.txt, ticks.json, summary.json (stage log, occupancies,
line heights obs/calc).

Usage: python workflow/scripts/refine_0615_ordering.py --scan <file> --occ pred --out <dir> [--md]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from synthesizability.parsers.xrd import parse_xrd_file  # noqa: E402
from synthesizability.phaseid import refine as rf  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
LIB = None   # folder holding the starting CIFs, set from --space
# Default partners: the two of the 0540 fit; --partners replaces them with
# library keys (the annealed scans' own search picks them)
PARTNERS = ['Y3Al2Cu7_alex_agm005614889_dis', 'Y3Al8Cu3_alex_agm005688824']
START_CELL = (5.2707, 4.0903)     # 0540 rank 90313 refined cell
OCC = {'pred': (1.0, 1.0), 'random': (0.6, 0.4), 'free': (0.9, 0.9),   # (3g Al, 2c Cu)
       'order': (0.92, 0.88)}   # S = 0.8 start; composition held at YAl3Cu2
# --occ order: one long-range order parameter S at the YAl3Cu2 composition,
# Al(3g) = 0.6 + 0.4 S, Cu(2c) = 0.4 + 0.6 S (S = 1 predicted, S = 0 random);
# four fractions, three constraints (two site sums and 3 Al(3g) + 2 Al(2c) = 3)
LINES = [(0, 0, 1), (1, 0, 1), (1, 1, 0), (2, 0, 0), (1, 1, 1), (0, 0, 2), (2, 0, 1), (0, 0, 3)]
WAVELENGTH = 1.5406
LIMITS = (8.0, 149.98)            # zero-background slide: nothing usable below 8 deg
MIN_GAIN = 0.05
MAX_CYCLES = 20


def others_all(r):
    return range(1, len(r.phases))


def exchange_cif(src, dst, lab_a, lab_b, x0=0.95):
    """Copy a CIF and split two single-element sites: site A gets element B
    and site B gets element A, at fraction 1 - x0 (start eta = 2 x0 - 1)."""
    lines = Path(src).read_text().splitlines()
    rows = {}
    for i, l in enumerate(lines):
        t = l.split()
        if len(t) >= 7 and t[1] in (lab_a, lab_b):
            rows[t[1]] = (i, t)
    assert set(rows) == {lab_a, lab_b}, f'--exchange sites {lab_a}, {lab_b} not both in {src}'
    (ia, ta), (ib, tb) = rows[lab_a], rows[lab_b]
    assert ta[2] == tb[2], '--exchange needs two sites of equal multiplicity'
    def row(t, el, lab, occ):
        return '  ' + '  '.join([el, lab] + t[2:6] + [f'{occ:.4f}'] + t[7:])
    lines[ia] = row(ta, ta[0], ta[1], x0)
    lines[ib] = row(tb, tb[0], tb[1], x0)
    extra = [row(ta, tb[0], tb[0] + 'x' + ta[1], 1 - x0), row(tb, ta[0], ta[0] + 'x' + tb[1], 1 - x0)]
    last = max(ia, ib)
    lines = lines[:last + 1] + extra + lines[last + 1:]
    Path(dst).write_text('\n'.join(lines) + '\n')
    return dst


def main_cif(path, al3g, cu2c):
    """CaCu5-type YAl3Cu2 with Al/Cu on both sites (labels Al3g, Cu3g, Cu2c, Al2c)."""
    a, c = START_CELL
    atoms = [('Y1a', 'Y', 0, 0, 0, 1.0),
             ('Al3g', 'Al', 0.5, 0, 0.5, al3g), ('Cu3g', 'Cu', 0.5, 0, 0.5, 1 - al3g),
             ('Cu2c', 'Cu', 1 / 3, 2 / 3, 0, cu2c), ('Al2c', 'Al', 1 / 3, 2 / 3, 0, 1 - cu2c)]
    lines = ['data_YAl3Cu2_CaCu5', "_symmetry_space_group_name_H-M   'P 6/m m m'", '_symmetry_Int_Tables_number   191',
             f'_cell_length_a   {a}', f'_cell_length_b   {a}', f'_cell_length_c   {c}',
             '_cell_angle_alpha   90', '_cell_angle_beta   90', '_cell_angle_gamma   120',
             'loop_', '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x', '_atom_site_fract_y',
             '_atom_site_fract_z', '_atom_site_occupancy', '_atom_site_U_iso_or_equiv']
    lines += [f'{lab} {el} {x:.6f} {y:.6f} {z:.6f} {occ:.4f} 0.005' for lab, el, x, y, z, occ in atoms if occ > 0]
    Path(path).write_text('\n'.join(lines) + '\n')
    return path


def line_heights(r):
    """Per cent of the maximum of (obs - bkg): observed and calculated
    maximum within 0.15 deg of each main-phase line, at the position of the
    refined hexagonal cell (Kalpha1).  Calculated includes the partners, so a
    line they overlap is marked."""
    x = r.hist.getdata('x'); yo = r.hist.getdata('yobs'); yc = r.hist.getdata('ycalc'); bg = r.hist.getdata('background')
    m = (x >= LIMITS[0]) & (x <= LIMITS[1])
    top = np.max((yo - bg)[m])
    cell = r.phases[0].get_cell(); a, c = cell['length_a'], cell['length_c']
    others = np.concatenate([v for k, v in r.reflections().items() if k != r.phases[0].name] or [np.zeros(0)])
    out = {}
    for h, k, l in LINES:
        d = 1 / np.sqrt(4 / 3 * (h * h + h * k + k * k) / a ** 2 + l * l / c ** 2)
        t = float(np.degrees(2 * np.arcsin(WAVELENGTH / (2 * d))))
        w = m & (abs(x - t) < 0.15)
        out[f'({h}{k}{l})'] = dict(two_theta=round(t, 2), obs=round(100 * float(np.max((yo - bg)[w])) / top, 1),
                                   calc=round(100 * float(np.max((yc - bg)[w])) / top, 1),
                                   overlap=bool(np.any(abs(others - t) < 0.15)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scan', type=Path, required=True)
    ap.add_argument('--occ', choices=list(OCC), required=True)
    ap.add_argument('--md', action='store_true')
    ap.add_argument('--md-hold', type=float, default=None, metavar='RATIO',
                    help='March-Dollase on the main phase with the ratio held at RATIO (not refined); implies --md. '
                         'Default off: --md alone refines the ratio as before')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--uiso', type=float, default=None,
                    help='hold one Uiso for the main phase at this value instead of refining it')
    ap.add_argument('--uiso-max', type=float, default=0.04,
                    help='upper limit for the refined Uiso (0.08: an effective term that also takes up the thin-layer intensity loss at high angle)')
    ap.add_argument('--main-key', default=None,
                    help='library key used as the main phase in place of the CaCu5 cell (--occ pred only), e.g. the '
                         'hull-stable Pmmm Y2(AlCu)5 ordering of the same lattice')
    ap.add_argument('--space', required=True, help='folder holding the starting CIFs <key>.cif (each frozen fit keeps its own in start/), '
                    'or the name of a search-match library space under results/searchmatch/library')
    ap.add_argument('--md-axis', type=int, nargs=3, default=[0, 0, 1], help='March-Dollase axis (default 0 0 1)')
    ap.add_argument('--hold-cells', type=int, nargs='*', default=[],
                    help='partner numbers (1-based) whose cells are held (a sharp minor line the cell would walk away from)')
    ap.add_argument('--gen-strain', action='store_true',
                    help='Stephens generalized microstrain on the main phase (after isotropic strain/size)')
    ap.add_argument('--uni-size', type=int, nargs=3, default=None,
                    help='uniaxial crystallite size on the main phase along this axis (after the strain model)')
    ap.add_argument('--tie-partner-broadening', type=int, nargs='*', default=[], metavar='N',
                    help='partner N takes the main phase\'s refined isotropic size and microstrain (held, not refined): '
                         'a matched contrast in which an added phase cannot absorb misfit through its own broadening')
    ap.add_argument('--partner-mustrain', nargs=2, action='append', default=[], metavar=('N', 'VALUE'),
                    help='starting isotropic microstrain for partner N (a broad twin of the main phase)')
    ap.add_argument('--inst-refine', nargs='*', default=[], choices=['U', 'V', 'W', 'X', 'Y', 'SH/L'],
                    help='instrument profile terms refined after the strain model (the quartz profile carries '
                         'transparency broadening a strongly absorbing sample lacks)')
    ap.add_argument('--inst-uv0', action='store_true',
                    help='with --inst-refine W: set U = V = 0 first (constant instrument Gaussian)')
    ap.add_argument('--free-mix', action='store_true', help='refine the main-phase Lorentz/Gauss mixing (bounded 0-1)')
    ap.add_argument('--twin', type=float, default=None, metavar='MUSTRAIN',
                    help='bimodal model after convergence: add a copy of the main phase at its refined cell and '
                         'positions with this isotropic microstrain, then refine its scale, its width, both widths '
                         'and scales, and its cell (one block at a time)')
    ap.add_argument('--positions', action='store_true', help='refine the main-phase atom coordinates (one block, after Uiso)')
    ap.add_argument('--instprm', type=Path, default=None,
                    help='instrument profile file (default instruments/siemens_d500.instprm, from the quartz standard)')
    ap.add_argument('--limits', type=float, nargs=2, default=None,
                    help='2theta window (default 8 to 149.98; the 0540 circular-holder scan ends at 120)')
    ap.add_argument('--thin', action='store_true',
                    help='finite-thickness (thin layer) intensity correction: Suortti form with SurfRoughA held at '
                         '1000, so f ~ 1 - exp(-SurfRoughB / sin theta), SurfRoughB = 2 mu t refined')
    ap.add_argument('--partner-sites', type=int, default=None,
                    help='partner number (1-based) whose Al/Cu split sites are refined, each pair summing to 1')
    ap.add_argument('--partners', nargs='*', default=PARTNERS, help='library keys in results/phaseid/library/Al-Cu-Y')
    ap.add_argument('--partner-scales', nargs='*', default=[],
                    help="a/b/c start factors per partner, e.g. 1.010/0.980/1.000 (the search's start_scale column)")
    ap.add_argument('--zero', action='store_true', help='DIAG: refine the Zero offset with the main cell instead of Shift (pipeline protocol)')
    ap.add_argument('--zero-and-shift', action='store_true', help='DIAG: refine Zero and Shift with the main cell')
    ap.add_argument('--hold-shift', action='store_true', help='hold sample displacement at 0 (sample at the reference height); cells take up any line shift')
    ap.add_argument('--strict-final', action='store_true', help='DIAG: final cycle kept only if it does not raise Rwp; re-refine after a clamp')
    ap.add_argument('--main-window', type=float, default=None, help='DIAG: main cell window (default 0.02; pipeline used 0.05 for _proto_/_dis keys)')
    ap.add_argument('--no-uiso', action='store_true', help='DIAG: skip the Uiso stage (pipeline protocol)')
    ap.add_argument('--no-bkg12', action='store_true', help='DIAG: skip the 12-term background stage')
    ap.add_argument('--accept-bound', action='store_true', help='DIAG: accept a microstrain clamped at its 60000 bound (pipeline behaviour)')
    ap.add_argument('--bkg-peak', type=float, nargs=2, default=None, metavar=('POS', 'FWHM'), help='DIAG: one refinable pseudo-Voigt background peak (amorphous hump of the slide/holder)')
    ap.add_argument('--hold-main-cell', action='store_true', help='DIAG: main cell held during the broadening stages (pipeline), freed again before Uiso')
    ap.add_argument('--exchange', nargs=2, default=None, metavar=('SITE_A', 'SITE_B'),
                    help='one exchange parameter between two single-element main-phase sites of equal multiplicity '
                         '(e.g. Zr0 Si1): each site gets the other element, composition held; eta = 2x - 1, 1 ordered, 0 random')
    ap.add_argument('--joint-strain', action='store_true', help='DIAG: one microstrain stage for every phase >= 1 wt%% (pipeline stage 4), before the per-phase stages')
    a = ap.parse_args()
    global LIMITS, LIB
    LIB = Path(a.space) if Path(a.space).is_dir() else ROOT / 'results/searchmatch/library' / a.space
    if a.limits:
        LIMITS = tuple(a.limits)
    partners = [(k, LIB / f'{k}.cif', 0.05 if rf.is_solution_key(k) else 0.02) for k in a.partners]
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    dat = parse_xrd_file(a.scan.resolve())
    xy = rf.write_xy(dat['two_theta'], dat['intensity'], out / 'scan.xy')
    al3g, cu2c = OCC[a.occ]
    if a.main_key:
        assert a.occ == 'pred', '--main-key takes the cell as given: --occ pred'
        cif, main_name = LIB / f'{a.main_key}.cif', a.main_key
    else:
        cif, main_name = main_cif(out / 'YAl3Cu2_start.cif', al3g, cu2c), 'YAl3Cu2'
    if a.exchange:
        cif = exchange_cif(cif, out / 'main_exchange.cif', *a.exchange)
    r = rf.Refinement(xy, a.instprm.resolve() if a.instprm else rf.instprm_for('siemens'), [str(cif)] + [str(p[1]) for p in partners],
                      names=[main_name] + [p[0] for p in partners], workdir=out / 'work', limits=LIMITS, n_bkg=8)
    r.gpx.set_Controls('cycles', MAX_CYCLES)
    # Start partner cells where the search's cell scan put them
    for i in others_all(r):
        for at in r.phases[i].atoms():
            if not at.uiso or at.uiso <= 0:
                at.uiso = 0.005             # COD CIFs can carry U = 0
    for n_, val in a.partner_mustrain:
        r.phases[int(n_)].data['Histograms'][r.hist.name]['Mustrain'][1][0] = float(val)
    for i, fac in enumerate(a.partner_scales, start=1):
        f = [float(v) for v in fac.split('/')]
        cell = r.phases[i].data['General']['Cell']
        for j in range(3):
            cell[1 + j] *= f[j]
        cell[7] *= f[0] * f[1] * f[2]
    if a.bkg_peak:
        pos, fw = a.bkg_peak
        bk = r.hist.data['Background'][1]
        bk['nPeaks'] = 1
        bk['peaksList'] = [[pos, True, 2.0e5, True, (100 * fw) ** 2 / (8 * np.log(2)), True, 0.1, False]]
    r.save()
    hname = r.hist.name
    windows = [a.main_window or 0.02] + [p[2] for p in partners]
    others = range(1, len(r.phases))
    cells0 = [ph.get_cell() for ph in r.phases]
    log = []
    main = lambda: r.phases[0]                        # noqa: E731  (re-read after every rollback)
    hap = lambda i=0: r.phases[i].data['Histograms'][hname]   # noqa: E731
    idx = {at.label: n for n, at in enumerate(main().atoms())}

    def snapshot():
        res = r.results(cells0, windows)
        return dict(wt=[round(p['wt'], 4) for p in res['phases']],
                    cells=[{k: round(v, 5) for k, v in p['cell'].items() if k in ('length_a', 'length_b', 'length_c')} for p in res['phases']],
                    shift=r.hist.data['Sample Parameters']['Shift'][0],
                    size=hap()['Size'][1][0], mustrain=hap()['Mustrain'][1][0],
                    md=hap()['Pref.Ori.'][1] if (hap()['Pref.Ori.'][2] or a.md_hold is not None) else None,
                    atoms={at.label: dict(occ=round(float(at.occupancy), 4), uiso=round(float(at.uiso), 5)) for at in main().atoms()},
                    cell_ok=not any(p['diverged'] for p in res['phases']))

    def stage(name, apply, force=False, accept=None):
        r.save()
        before = r.hist.get_wR() if log else None
        apply()
        try:
            wR = r._refine()
        except rf.RefinementFailed as exc:
            print(f'{name}: failed ({exc}), rolled back', flush=True)
            r._reload(); log.append(dict(stage=name, wR=None, kept=False, note=str(exc))); return False
        if r._clamp_broadening() and a.strict_final:
            try:
                wR = r._refine()
            except rf.RefinementFailed as exc:
                r._reload(); log.append(dict(stage=name, wR=None, kept=False, note='after clamp: ' + str(exc))); return False
            r._clamp_broadening()
        kept = force or before is None or wR < before - MIN_GAIN
        if a.strict_final and force and before is not None and wR > before + MIN_GAIN and name.startswith('final'):
            kept = False
        note = '' if kept else f'no gain over {before:.2f}'
        snap = snapshot()
        if kept and not snap['cell_ok']:
            res = r.results(cells0, windows)
            kept, note = False, 'cell outside its window: ' + ', '.join(p['name'] for p in res['phases'] if p['diverged'])
        elif kept and accept is not None and not accept():
            kept, note = False, 'unphysical value'
        print(f'{name}: Rwp {wR:.2f}' + ('' if kept else f' (rolled back: {note})'), flush=True)
        print(f'    {snap}', flush=True)
        log.append(dict(stage=name, wR=round(wR, 3), kept=bool(kept), note=note, state=snap))
        if not kept:
            r._reload()
        return kept

    def mix1(i):
        h = hap(i)
        h['Size'][1][2] = h['Mustrain'][1][2] = 1.0
        h['Size'][2][2] = h['Mustrain'][2][2] = False

    stage('scale + background (8 Chebyshev)', lambda: r.hist.set_refinements({'Sample Parameters': ['Scale']}), force=True)
    stage('partner phase fractions', lambda: [r.phases[i].set_HAP_refinements({'Scale': True}) for i in others], force=True)
    # One cell at a time: all three together broke the R-3m partner's
    # metric tensor in an early run
    def cell_disp():
        main().set_refinements({'Cell': True})
        if a.zero or a.zero_and_shift:
            r.hist.set_refinements({'Instrument Parameters': ['Zero']})
        if not a.zero and not a.hold_shift:
            r.hist.set_refinements({'Sample Parameters': ['Shift']})
    stage('main cell + displacement (Shift/Zero)', cell_disp, force=True)
    for i in others:
        if r.weight_fractions()[i] < 0.01 or i in a.hold_cells:   # a phase at zero scale has no defined cell
            continue
        stage(f'{r.phases[i].name} cell', lambda i=i: r.phases[i].set_refinements({'Cell': True}))
        # then held: a minor partner's cell drifted out of its window in
        # later stages on the annealed_2x scan
        r.phases[i].data['General']['Cell'][0] = False
        r.save()

    def strain(i):
        def f():
            mix1(i)
            r._bound(r.phases[i], 'Mustrain;i', 0.0, 60000.0)
            r.phases[i].set_HAP_refinements({'Mustrain': {'type': 'isotropic', 'refine': True}})
        return f
    if a.hold_main_cell:
        main().data['General']['Cell'][0] = False; r.save()
    if a.joint_strain:
        def joint():
            for i in range(len(r.phases)):
                if i == 0 or r.weight_fractions()[i] >= 0.01:
                    strain(i)()
        stage('joint isotropic microstrain, all phases', joint,
              accept=lambda: all(0 < hap(i)['Mustrain'][1][0] <= 60000 for i in range(len(r.phases))))
    stage('main isotropic microstrain', strain(0), accept=lambda: 0 < hap()['Mustrain'][1][0] < 60000 or (a.accept_bound and hap()['Mustrain'][1][0] == 60000))

    def size():
        hap()['Size'][1][0] = 0.5
        r._bound(main(), 'Size;i', 0.02, 10.0)
        main().set_HAP_refinements({'Size': {'type': 'isotropic', 'refine': True}})
    stage('main isotropic size', size, accept=lambda: 0.025 <= hap()['Size'][1][0] < 10.0)   # at the 0.02 bound = not a fit
    for i in others:
        if i in a.tie_partner_broadening:
            # matched contrast: copy the main phase's broadening and hold it
            hp = r.phases[i].data['Histograms'][r.hist.name]
            hp['Mustrain'][1][0] = hap()['Mustrain'][1][0]
            hp['Size'][1][0] = hap()['Size'][1][0]
            hp['Mustrain'][2][0] = False
            hp['Size'][2][0] = False
            r.save()
            continue
        if r.weight_fractions()[i] >= 0.03:
            stage(f'{r.phases[i].name} microstrain', strain(i), accept=lambda i=i: 0 < hap(i)['Mustrain'][1][0] < 60000 or (a.accept_bound and hap(i)['Mustrain'][1][0] == 60000))

    if a.hold_main_cell:
        main().data['General']['Cell'][0] = True; r.save()
    if a.gen_strain:
        def generalized_strain():
            from GSASII import GSASIIspc as G2spc
            g = main().data['General']
            mu = hap()['Mustrain'][1][0]
            hap()['Mustrain'][4] = list(G2spc.Muiso2Shkl(mu, g['SGData'], g['Cell'][1:7]))
            main().set_HAP_refinements({'Mustrain': {'type': 'generalized', 'refine': True}})
        stage('generalized (Stephens) microstrain, main phase', generalized_strain)
    if a.uni_size:
        def uniaxial_size():
            sz = hap()['Size'][1][0]
            hap()['Size'][1][:2] = [sz, sz]
            r._bound(main(), 'Size;a', 0.02, 10.0)
            main().set_HAP_refinements({'Size': {'type': 'uniaxial', 'direction': list(a.uni_size), 'refine': True}})
        stage(f'uniaxial size along {tuple(a.uni_size)}, main phase', uniaxial_size,
              accept=lambda: all(0.02 <= v < 10 for v in hap()['Size'][1][:2]))
    if a.inst_refine:
        def gauss_ok():
            ip = r.hist.data['Instrument Parameters'][0]
            th = np.radians(np.linspace(LIMITS[0], LIMITS[1], 50) / 2)
            sig2 = ip['U'][1] * np.tan(th) ** 2 + ip['V'][1] * np.tan(th) + ip['W'][1]
            return bool(np.all(sig2 > 0))
        def inst():
            if a.inst_uv0:
                ip = r.hist.data['Instrument Parameters'][0]
                th = np.radians(40.0 / 2)      # keep the 40 deg Gaussian width as the starting W
                ip['W'][1] = ip['U'][1] * np.tan(th) ** 2 + ip['V'][1] * np.tan(th) + ip['W'][1]
                ip['U'][1] = ip['V'][1] = 0.0
            r.hist.set_refinements({'Instrument Parameters': list(a.inst_refine)})
        stage(f"instrument profile {' '.join(a.inst_refine)}" + (' (U = V = 0)' if a.inst_uv0 else ''), inst, accept=gauss_ok)
    if a.free_mix:
        def free_mix():
            r._bound(main(), 'Size;mx', 0.0, 1.0); r._bound(main(), 'Mustrain;mx', 0.0, 1.0)
            hap()['Size'][2][2] = hap()['Mustrain'][2][2] = True
        stage('Lorentz/Gauss mixing, main phase', free_mix,
              accept=lambda: 0 <= hap()['Size'][1][2] <= 1 and 0 <= hap()['Mustrain'][1][2] <= 1)
    if a.md_hold is not None:
        # Ratio held (0605: free MD ran to a worse minimum; the verdict fit held (001) at 0.70)
        def po_held():
            p = hap()['Pref.Ori.']
            p[0] = 'MD'; p[1] = a.md_hold; p[2] = False; p[3] = list(a.md_axis)
        stage(f'March-Dollase {tuple(a.md_axis)} held at {a.md_hold}, main phase', po_held, force=True)
    elif a.md:
        def po():
            p = hap()['Pref.Ori.']
            p[0] = 'MD'; p[1] = 1.0; p[2] = True; p[3] = list(a.md_axis)
        # Kept even without gain: the variant asks what the occupancies do
        # with texture free, so the term stays in
        stage(f'March-Dollase {tuple(a.md_axis)}, main phase', po, force=True, accept=lambda: 0.3 <= hap()['Pref.Ori.'][1] <= 3.0)

    def site_fractions(sites):
        """Refine the fractions of the listed sites; each split pair sums to 1
        and shares one Uiso."""
        def f():
            pairs = {'3g': ('Al3g', 'Cu3g'), '2c': ('Cu2c', 'Al2c')}
            for s in sites:
                a1, a2 = pairs[s]
                main().set_refinements({'Atoms': {a1: 'F', a2: 'F'}})
                r.gpx.add_EqnConstr(1.0, [f'0::Afrac:{idx[a1]}', f'0::Afrac:{idx[a2]}'], [1, 1])
        return f
    in_range = lambda: all(-0.02 <= at.occupancy <= 1.02 for at in main().atoms())   # noqa: E731
    if a.occ == 'order':
        def order_parameter():
            main().set_refinements({'Atoms': {l: 'F' for l in ('Al3g', 'Cu3g', 'Cu2c', 'Al2c')}})
            v = {l: f'0::Afrac:{idx[l]}' for l in ('Al3g', 'Cu3g', 'Cu2c', 'Al2c')}
            r.gpx.add_EqnConstr(1.0, [v['Al3g'], v['Cu3g']], [1, 1])
            r.gpx.add_EqnConstr(1.0, [v['Cu2c'], v['Al2c']], [1, 1])
            r.gpx.add_EqnConstr(3.0, [v['Al3g'], v['Al2c']], [3, 2])
        stage('order parameter S (composition held)', order_parameter, force=True, accept=in_range)
    if a.occ == 'free':
        if not stage('site fractions 3g + 2c', site_fractions(['3g', '2c']), accept=in_range):
            # GSAS-II limits do not act on constrained fractions; hold the
            # site that ran out of range at its bound and refine the other
            for held, other, (lab, val) in (('3g', '2c', ('Al3g', 1.0)), ('2c', '3g', ('Cu2c', 1.0))):
                def one(other=other, lab=lab, val=val):
                    partner = {'Al3g': 'Cu3g', 'Cu2c': 'Al2c'}[lab]
                    for at in main().atoms():
                        if at.label == lab: at.occupancy = val
                        if at.label == partner: at.occupancy = 1 - val
                    site_fractions([other])()
                if stage(f'site fractions: {held} held ({lab} {val}), {other} free', one, accept=in_range):
                    break

    if a.exchange:
        la, lb = a.exchange
        ta = next(at.type for at in main().atoms() if at.label == la)
        tb = next(at.type for at in main().atoms() if at.label == lb)
        xa, xb = tb + 'x' + la, ta + 'x' + lb       # element B on site A, element A on site B
        def exchange_parameter():
            main().set_refinements({'Atoms': {l: 'F' for l in (la, lb, xa, xb)}})
            v = {l: f'0::Afrac:{idx[l]}' for l in (la, lb, xa, xb)}
            r.gpx.add_EqnConstr(1.0, [v[la], v[xa]], [1, 1])     # site A full
            r.gpx.add_EqnConstr(1.0, [v[lb], v[xb]], [1, 1])     # site B full
            r.gpx.add_EqnConstr(0.0, [v[la], v[lb]], [1, -1])    # composition held (equal multiplicity)
        stage(f'{ta}/{tb} exchange parameter ({la}, {lb})', exchange_parameter, force=True, accept=in_range)

    def uiso():
        main().set_refinements({'Atoms': {at.label: 'U' for at in main().atoms()}})
        labs = [at.label for at in main().atoms()]
        if len(labs) > 1:  # a one-site cell (e.g. fcc) needs no equivalence
            r.gpx.add_EquivConstr([f'0::AUiso:{idx[l]}' for l in labs])
    if a.thin:
        # GSAS-II has no thin-sample term for Bragg-Brentano; its Suortti
        # roughness f = A + (1 - A) exp(-B / sin theta) tends to the thin-layer
        # factor 1 - exp(-2 mu t / sin theta) (up to scale) as A grows
        def thin():
            sp = r.hist.data['Sample Parameters']
            sp['SurfRoughA'] = [1000.0, False]
            sp['SurfRoughB'] = [1.0, True]
        stage('thin-layer correction (SurfRoughB = 2 mu t)', thin, force=True,
              accept=lambda: 0.02 < r.hist.data['Sample Parameters']['SurfRoughB'][0] < 10)
    if a.no_uiso:
        pass
    elif a.uiso is not None:
        def hold_uiso():
            for at in main().atoms():
                at.uiso = a.uiso
        stage(f'Uiso held at {a.uiso}', hold_uiso, force=True)
    else:
        stage('one Uiso, main phase', uiso, accept=lambda: all(0 < at.uiso < a.uiso_max for at in main().atoms()))
    if a.positions:
        # free coordinates only (GSAS-II holds special-position components by symmetry)
        stage('main-phase atom positions', lambda: main().set_refinements({'Atoms': {at.label: 'X' for at in main().atoms()}}),
              accept=lambda: True)
    if a.partner_sites:
        pi = a.partner_sites
        def partner_fractions():
            ph = r.phases[pi]
            atoms = ph.atoms()
            pos = {}
            for n, at in enumerate(atoms):
                pos.setdefault(tuple(round(float(c), 4) for c in at.coordinates), []).append((n, at))
            for group in pos.values():
                if len(group) == 2 and {g[1].element for g in group} == {'Al', 'Cu'}:
                    ph.set_refinements({'Atoms': {g[1].label: 'F' for g in group}})
                    r.gpx.add_EqnConstr(1.0, [f'{ph.id}::Afrac:{g[0]}' for g in group], [1, 1])
        stage(f'{r.phases[pi].name} Al/Cu site fractions', partner_fractions,
              accept=lambda: all(-0.02 <= at.occupancy <= 1.02 for at in r.phases[pi].atoms()))
    if not a.no_bkg12: stage('background 12 Chebyshev', lambda: r.hist.set_refinements({'Background': {'type': 'chebyschev-1', 'no. coeffs': 12, 'refine': True}}))
    stage('final cycle, all kept parameters', lambda: None, force=True)
    if a.twin:
        # Bimodal microstructure (sharp cores + broad tails): the usual Rietveld
        # treatment is two copies of one structure with different broadening
        r.save()
        with rf.quiet():
            r.gpx.add_phase(str(cif), phasename=main_name + '_broad', histograms=[r.hist])
        r.phases = r.gpx.phases()
        tw = len(r.phases) - 1
        m0, t0 = r.phases[0], r.phases[tw]
        t0.data['General']['Cell'] = list(m0.data['General']['Cell'])
        t0.data['General']['Cell'][0] = False
        for am, at in zip(m0.atoms(), t0.atoms()):
            at.coordinates = list(am.coordinates); at.uiso = am.uiso
        th = t0.data['Histograms'][hname]
        th['Scale'] = [0.1 * m0.data['Histograms'][hname]['Scale'][0], False]
        th['Mustrain'][1][0] = a.twin
        th['Mustrain'][1][2] = th['Size'][1][2] = 1.0
        windows.append(0.01); cells0.append(t0.get_cell())
        r.save()
        stage('broad copy: scale', lambda: r.phases[tw].set_HAP_refinements({'Scale': True}), force=True)

        def twin_width():
            r._bound(r.phases[tw], 'Mustrain;i', 0.0, 60000.0)
            r.phases[tw].set_HAP_refinements({'Mustrain': {'type': 'isotropic', 'refine': True}})
        stage('broad copy: isotropic microstrain', twin_width, accept=lambda: 0 < hap(tw)['Mustrain'][1][0] < 60000)
        # the main phase's kept width terms are still flagged, so this cycle
        # refines both widths with the copy's scale (the main scale is the
        # histogram scale; its own HAP scale would be collinear)
        both = lambda: min(r.weight_fractions()[[0, tw]]) > 0.05   # noqa: E731
        stage('joint cycle: both widths, copy scale', lambda: None, force=True, accept=both)
        stage('broad copy: cell (1 %)', lambda: r.phases[tw].set_refinements({'Cell': True}), accept=both)
        stage('final cycle with the broad copy', lambda: None, force=True, accept=both)

    r.save()
    snap = snapshot()
    occ = {k: v['occ'] for k, v in snap['atoms'].items()}
    x_al = 3 * occ.get('Al3g', 0) + 2 * occ.get('Al2c', 0)
    res = r.results(cells0, windows)
    lines = {} if a.main_key else line_heights(r)
    if a.main_key:
        x_al = None
    partner_occ = ({at.label: round(float(at.occupancy), 3) for at in r.phases[a.partner_sites].atoms()}
                   if a.partner_sites else None)
    S = None
    if a.occ == 'order':
        S = round((next(at.occupancy for at in main().atoms() if at.label == 'Al3g') - 0.6) / 0.4, 3)
    record = dict(zero=r.hist.data['Instrument Parameters'][0]['Zero'][1], instprm=str(a.instprm) if a.instprm else 'siemens_d500.instprm', inst={k: r.hist.data['Instrument Parameters'][0][k][1] for k in 'UVWXY'},
                  mix=(hap()['Size'][1][2], hap()['Mustrain'][1][2]), main=main_name, space=a.space, md_axis=a.md_axis if (a.md or a.md_hold is not None) else None, partner_occ=partner_occ, S=S,
                  two_mu_t=r.hist.data['Sample Parameters']['SurfRoughB'][0] if a.thin else None, scan=a.scan.name, occ=a.occ, md=a.md, uiso_held=a.uiso, uiso_max=a.uiso_max, partners=a.partners, partner_scales=a.partner_scales, wR=r.hist.get_wR(), gof=res.get('gof'),
                  x_al=None if x_al is None else round(x_al, 3), lines=lines, state=snap)
    if a.md_hold is not None:
        record['md_hold'] = a.md_hold
    r.export_pattern(out / 'fit.txt')
    (out / 'ticks.json').write_text(json.dumps({k: [float(v) for v in arr] for k, arr in r.reflections().items()}))
    (out / 'summary.json').write_text(json.dumps(dict(record=record, limits=LIMITS, stages=log), indent=1, default=float))
    print(json.dumps(record, indent=1, default=float))


if __name__ == '__main__':
    main()
