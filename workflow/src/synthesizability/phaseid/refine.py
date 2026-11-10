"""GSAS-II refinement of a mixture hypothesis.

Uses the GSAS-II checkout and compiled binaries that the ``rietveld`` repo
set up (``GSAS2_DIR``; binaries in ``~/.GSASII/bin``).  Every candidate is
refined with the same fixed sequence so that the fit statistics are
comparable:

1. histogram scale and Chebyshev background
2. plus zero offset and unit cells
3. plus phase fractions (free mode) or nothing (fixed lever-rule mode)
4. plus isotropic size and microstrain per phase

Instrument profile parameters are never refined on candidates; they come
from an ``.instprm`` fitted once per diffractometer (``calibrate_profile``).

Fixed lever-rule mode uses the GSAS-II convention that the weight fraction
of phase *i* is ``S_i M_i / sum_j S_j M_j`` with ``M_i`` the unit-cell mass:
phase scales are set to ``w_i / M_i`` and not refined, and only the
histogram scale floats.
"""
from __future__ import annotations

import contextlib
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

# GSAS-II: either importable as a package, or point GSAS2_DIR at its GSAS-II source folder
GSAS2_DIR = Path(os.environ['GSAS2_DIR']) if os.environ.get('GSAS2_DIR') else None
if GSAS2_DIR is not None and str(GSAS2_DIR) not in sys.path:
    sys.path.insert(0, str(GSAS2_DIR))

from ..paths import ROOT, INSTPRM_DIR
PANALYTICAL_INSTPRM = INSTPRM_DIR / 'panalytical_xpert_lab6.instprm'

_G2sc = None


class RefinementTimeout(RuntimeError):
    pass


@contextlib.contextmanager
def time_limit(seconds):
    """Raise RefinementTimeout if the block runs longer than ``seconds``
    (SIGALRM; main thread only, which is how the batch workers run)."""
    import signal
    if not seconds:
        yield
        return

    def handler(signum, frame):
        raise RefinementTimeout(f'refinement exceeded {seconds} s')
    old = signal.signal(signal.SIGALRM, handler)
    signal.alarm(int(seconds))
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


class _QuickStop(Exception):
    """Internal: leave ``Refinement.run`` after the screening stages."""


class RefinementFailed(RuntimeError):
    """GSAS-II reported a failed least-squares cycle."""


def g2sc():
    """Import GSASIIscriptable lazily and quietly."""
    global _G2sc
    if _G2sc is None:
        with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
            from GSASII import GSASIIscriptable as G2sc
        _G2sc = G2sc
    return _G2sc


@contextlib.contextmanager
def quiet():
    with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull), \
            contextlib.redirect_stderr(devnull):
        yield


def write_xy(two_theta, intensity, path):
    np.savetxt(path, np.column_stack([two_theta, intensity]), fmt='%.4f %.1f')
    return path


def is_solution_key(key):
    """Library keys of cells whose lattice parameters are estimates rather
    than relaxed or measured: disordered twins (``_dis``), Vegard prototype
    cells (``_proto_``), MatterSim-relaxed random cells (``_p2_``) and
    pair-merged interpolations (``_pm_``)."""
    key = str(key)
    return key.endswith('_dis') or '_proto_' in key or '_p2_' in key or '_pm_' in key


def cell_windows(keys, ordered=0.02, solution=0.05):
    """Per-phase cell-move windows for ``Refinement.run``: solution cells may
    move 5 %, ordered database and experimental cells 2 %.  A single window
    per refinement (5 % for every phase as soon as one solution cell is in the
    set) would let ordered cells drift 2 to 5 % without the fit being rejected."""
    return [solution if is_solution_key(k) else ordered for k in keys]


def instprm_for(instrument):
    """Instrument parameter file for 'siemens' or 'panalytical'."""
    if instrument == 'panalytical':
        p = INSTPRM_DIR / 'panalytical.instprm'
        return p if p.exists() else PANALYTICAL_INSTPRM
    return INSTPRM_DIR / 'siemens_d500.instprm'


class Refinement:
    """One GSAS-II project: a scan plus a list of phases."""

    def __init__(self, xy_file, instprm, cifs, names=None, workdir=None,
                 limits=(18.0, 118.0), n_bkg=8):
        G2sc = g2sc()
        self.workdir = Path(workdir or tempfile.mkdtemp(prefix='phaseid_'))
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.gpx_path = self.workdir / 'fit.gpx'
        with quiet():
            self.gpx = G2sc.G2Project(newgpx=str(self.gpx_path))
            self.hist = self.gpx.add_powder_histogram(str(xy_file), str(instprm))
            self.phases = []
            for n, cif in enumerate(cifs):
                name = names[n] if names else Path(cif).stem
                self.phases.append(self.gpx.add_phase(str(cif), phasename=name, histograms=[self.hist]))
        self.hist.set_refinements({'Limits': list(limits)})
        self.hist.set_refinements({'Background': {'type': 'chebyschev-1', 'no. coeffs': n_bkg, 'refine': True}})
        self.limits = limits
        self.log = []
        self.stage_s = []

    def _refine(self):
        """One GSAS-II refinement; raises RefinementFailed on the messages
        GSAS-II prints (without raising) when the least squares breaks.
        Wall time of every call is appended to ``self.stage_s``."""
        import io
        t0 = time.time()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            try:
                self.gpx.do_refinements([{}])
            except Exception as exc:
                raise RefinementFailed(str(exc)) from exc
        text = buf.getvalue()
        if 'Refinement failed' in text or 'Invalid metric tensor' in text:
            line = next((ln for ln in text.splitlines() if 'Error' in ln or 'ERROR' in ln), 'refinement failed')
            raise RefinementFailed(line.strip())
        self.log.append(text)
        self.stage_s.append(round(time.time() - t0, 1))
        wR = self.hist.get_wR()
        if wR is None or not np.isfinite(wR) or wR > 95:
            raise RefinementFailed(f'diverged (wR = {wR})')
        return wR

    # Mustrain 60000 is a Lorentzian width of about 1.3 deg at 2theta = 40 deg,
    # beyond anything in this study; Size 0.02 um is the lower limit of what
    # lab data can distinguish from strain.
    BOUNDS = {'Size': (0.02, 5.0), 'Mustrain': (0.0, 60000.0), 'mix': (0.0, 1.0)}

    def _bound(self, phase, parm, lo, hi):
        """Parameter limits for a histogram-phase variable such as 'Mustrain;i'.

        GSAS-II freezes a variable that leaves its limits but writes the phase
        data before clamping, so the out-of-range value survives; see
        ``_clamp_broadening`` for the fix applied after each stage."""
        var = f'{phase.id}:0:{parm}'
        self.gpx.set_Controls('parmMin', lo, variable=var)
        self.gpx.set_Controls('parmMax', hi, variable=var)

    def _clamp_broadening(self):
        """Reset size, strain and mixing values that left their bounds to the
        bound and stop refining them.  Returns the number clamped."""
        hname = self.hist.name
        n = 0
        for ph in self.phases:
            hap = ph.data['Histograms'][hname]
            for key in ('Size', 'Mustrain'):
                lo, hi = self.BOUNDS[key]
                vals, flags = hap[key][1], hap[key][2]
                if not lo <= vals[0] <= hi:
                    vals[0] = min(max(vals[0], lo), hi)
                    flags[0] = False
                    n += 1
                mlo, mhi = self.BOUNDS['mix']
                if not mlo <= vals[2] <= mhi:
                    vals[2] = min(max(vals[2], mlo), mhi)
                    flags[2] = False
                    n += 1
        return n

    def _reload(self):
        """Reload the state saved by the last ``save()``.  GSAS-II writes the
        gpx again after a refinement that has run to the end even when it has
        diverged, so the project file itself cannot be trusted: ``save()``
        keeps a snapshot copy and this restores it (without it, every late-stage
        divergence came back as wR 100)."""
        if self.snapshot.exists():
            shutil.copy(self.snapshot, self.gpx_path)
        self.gpx = type(self.gpx)(str(self.gpx_path))
        self.hist = self.gpx.histograms()[0]
        self.phases = self.gpx.phases()

    def weight_fractions(self):
        hname = self.hist.name
        sm = np.array([ph.data['Histograms'][hname]['Scale'][0] * ph.data['General']['Mass']
                       for ph in self.phases])
        return sm / sm.sum() if sm.sum() > 0 else np.zeros(len(sm))

    def set_fixed_fractions(self, weight_fractions):
        """Phase scales proportional to w_i / M_i, not refined."""
        for ph, w in zip(self.phases, weight_fractions):
            mass = ph.data['General']['Mass']
            for h in ph.data['Histograms']:
                ph.data['Histograms'][h]['Scale'] = [float(w) / mass, False]

    REFINE_MIXING = False   # the size + Lorentz/Gauss mixing stage was 70 % of a refinement for 0.7 in wR;
    MIXING_VALUE = 1.0      # GSAS-II tutorials leave the sample mixing at its default (1 = Lorentzian) and never refine it for lab data

    def run(self, free_fractions=True, size_strain=True, max_cell_change=0.02,
            min_wt_for_cell=0.05, quick=False, fixed_cells=(), refine_mixing=None, mixing_value=None):
        """Fixed refinement sequence; returns a results dict.

        ``max_cell_change`` is the cell-move window as a fraction (one number
        for every phase, or one per phase in the order of ``cifs``, see
        ``cell_windows``): a stage after which a phase's cell moved further
        is rolled back, and the final result is marked failed if a phase above
        ``min_wt_for_cell`` ends outside its window.

        ``quick=True`` stops after stage 2 (scale, background and phase
        fractions with cells and profile fixed): the screening stage of the
        candidate ladder, a few seconds per assemblage.

        1. histogram scale + background
        2. phase fractions (free mode only)
        3. zero offset + cells of phases above ``min_wt_for_cell``
        4. isotropic microstrain of those phases, bounded (the instrument
           profile is fixed, so sample broadening is absorbed here)
        5. isotropic size, bounded, plus the Lorentz/Gauss mixing of both
           terms (GSAS-II defaults to pure Lorentzian broadening, which cannot
           absorb Gaussian-like sample broadening)
        Stages 4 and 5 fall back to the previous converged state if they
        misbehave.  A GSAS-II failure elsewhere ends the sequence; the dict then
        has ``failed=True`` and ``error`` set, with the statistics of the last
        stage that converged.
        """
        t0 = time.time()
        cells0 = [ph.get_cell() for ph in self.phases]
        windows = self._per_phase(max_cell_change)
        out = {'stage_wR': [], 'error': None}
        try:
            self.hist.set_refinements({'Sample Parameters': ['Scale']})
            out['stage_wR'].append(self._refine())
            if free_fractions and len(self.phases) > 1:
                for ph in self.phases[1:]:
                    ph.set_HAP_refinements({'Scale': True})
                out['stage_wR'].append(self._refine())
            if quick:
                raise _QuickStop()
            wt = self.weight_fractions()
            # ``fixed_cells``: phases whose cell is held (P1 random-decoration
            # cells, which would otherwise refine six free parameters; they
            # start from the scanned cell instead)
            live = [ph for ph, w in zip(self.phases, wt) if w >= min_wt_for_cell and ph.name not in fixed_cells]
            self.hist.set_refinements({'Instrument Parameters': ['Zero']})
            for ph in live:
                ph.set_refinements({'Cell': True})
            self.save()
            try:
                out['stage_wR'].append(self._refine())
                moved = [ph.name for ph, c0, mc in zip(self.phases, cells0, windows)
                         if any(abs(ph.get_cell()[k] - c0[k]) / c0[k] > mc
                                for k in ('length_a', 'length_b', 'length_c'))]
                if moved:
                    raise RefinementFailed('cell ran away: ' + ', '.join(moved))
            except RefinementFailed as exc:
                # A minor phase's cell ran away; reload the pre-cell state and
                # free only the dominant phase's cell
                out['cell_retry'] = str(exc)
                self._reload()
                wt = self.weight_fractions()
                live = [self.phases[int(np.argmax(wt))]]
                for ph in self.phases:
                    ph.clear_refinements({'Cell': True})
                if live[0].name not in fixed_cells:
                    live[0].set_refinements({'Cell': True})
                out['stage_wR'].append(self._refine())
            if size_strain:
                live_names = {ph.name for ph, w in zip(self.phases, self.weight_fractions()) if w >= min_wt_for_cell}
                # Cells stay fixed while the broadening terms settle; they are
                # polished again at the end for the dominant phase only
                for ph in self.phases:
                    ph.clear_refinements({'Cell': True})
                # Stage 4: microstrain alone, bounded
                self.save()
                for ph in self.phases:
                    if ph.name in live_names:
                        self._bound(ph, 'Mustrain;i', *self.BOUNDS['Mustrain'])
                        ph.set_HAP_refinements({'Mustrain': {'type': 'isotropic', 'refine': True}})
                try:
                    out['stage_wR'].append(self._refine())
                    if self._clamp_broadening():
                        out['stage_wR'].append(self._refine())
                except RefinementFailed as exc:
                    out['mustrain_skipped'] = str(exc)
                    self._reload()
                # Stage 5: size and the Lorentz/Gauss mixing of both terms, bounded
                self.save()
                hname = self.hist.name
                if refine_mixing is None:
                    refine_mixing = self.REFINE_MIXING
                for ph in self.phases:
                    if ph.name not in live_names:
                        continue
                    self._bound(ph, 'Size;i', 0.02, 5.0)
                    ph.set_HAP_refinements({'Size': {'type': 'isotropic', 'refine': True}})
                    hap = ph.data['Histograms'][hname]
                    if refine_mixing:
                        self._bound(ph, 'Size;mx', 0.0, 1.0)
                        self._bound(ph, 'Mustrain;mx', 0.0, 1.0)
                        hap['Size'][2][2] = True
                        hap['Mustrain'][2][2] = True
                    else:
                        # half Lorentzian, half Gaussian, held: absorbs Gaussian-like
                        # sample broadening without the slow mixing refinement
                        mv = self.MIXING_VALUE if mixing_value is None else mixing_value
                        hap['Size'][1][2] = mv
                        hap['Mustrain'][1][2] = mv
                        hap['Size'][2][2] = False
                        hap['Mustrain'][2][2] = False
                try:
                    out['stage_wR'].append(self._refine())
                    clamped = self._clamp_broadening()
                    if clamped:
                        out['clamped'] = clamped
                        out['stage_wR'].append(self._refine())
                except RefinementFailed as exc:
                    out['size_lgmix_skipped'] = str(exc)
                    self._reload()
                # Stage 6: cell of the dominant phase once more, with the
                # broadening terms fixed
                self.save()
                hname = self.hist.name
                for ph in self.phases:
                    hap = ph.data['Histograms'][hname]
                    hap['Size'][2][0] = hap['Size'][2][2] = False
                    hap['Mustrain'][2][0] = hap['Mustrain'][2][2] = False
                wt = self.weight_fractions()
                if self.phases[int(np.argmax(wt))].name not in fixed_cells:
                    self.phases[int(np.argmax(wt))].set_refinements({'Cell': True})
                try:
                    out['stage_wR'].append(self._refine())
                except RefinementFailed as exc:
                    out['cell_polish_skipped'] = str(exc)
                    self._reload()
                # Values that left their bounds in the last cycle are reported
                # at the bound (the pattern is not recomputed for this)
                self._clamp_broadening()
        except _QuickStop:
            pass
        except RefinementFailed as exc:
            out['error'] = str(exc)
        stage_error = out['error']
        out.update(self.results(cells0, windows, min_wt=min_wt_for_cell))
        out['error'] = stage_error or out['error']
        out['failed'] = out['failed'] or out['error'] is not None
        out['seconds'] = time.time() - t0
        out['stage_s'] = list(self.stage_s)
        return out

    def _per_phase(self, value):
        """``value`` as one window per phase (a scalar applies to all)."""
        if np.isscalar(value):
            return [float(value)] * len(self.phases)
        value = [float(v) for v in value]
        if len(value) != len(self.phases):
            raise ValueError(f'{len(value)} cell windows for {len(self.phases)} phases')
        return value

    def results(self, cells0=None, max_cell_change=0.02, min_wt=0.0):
        """Fit statistics, weight fractions and cells.  A cell that moved by
        more than ``max_cell_change`` (scalar or per phase) marks the fit as
        failed, except for phases below ``min_wt`` (their cells are never
        refined on purpose and any drift is irrelevant)."""
        windows = self._per_phase(max_cell_change)
        res = dict(self.hist.residuals)
        hname = self.hist.name
        scales, masses, fracs = [], [], []
        for ph in self.phases:
            hap = ph.data['Histograms'][hname]
            scales.append(hap['Scale'][0])
            masses.append(ph.data['General']['Mass'])
        sm = np.array(scales) * np.array(masses)
        wt = sm / sm.sum() if sm.sum() > 0 else np.zeros(len(sm))
        cells = [ph.get_cell() for ph in self.phases]
        flags = []
        for n, (ph, cell) in enumerate(zip(self.phases, cells)):
            d = {}
            if cells0 is not None:
                c0 = cells0[n]
                for k in ('length_a', 'length_b', 'length_c'):
                    d[k[-1]] = (cell[k] - c0[k]) / c0[k]
            flags.append(dict(name=ph.name, wt=float(wt[n]), scale=float(scales[n]),
                              cell={k: float(v) for k, v in cell.items()},
                              cell_change=d,
                              diverged=(any(abs(v) > windows[n] for v in d.values())
                                        and float(wt[n]) >= min_wt) if d else False,
                              rf=res.get(f'{n}:0:Rf'), rf2=res.get(f'{n}:0:Rf^2')))
        n_obs = int(np.sum((self.hist.getdata('x') >= self.limits[0]) & (self.hist.getdata('x') <= self.limits[1])))
        n_var = self._n_refined()
        wR, wRmin = res.get('wR'), res.get('wRmin')
        gof = (wR / wRmin) if (wR and wRmin) else None
        return dict(wR=wR, Rp=res.get('R'), wRmin=wRmin, gof=gof, chi2=(gof ** 2 if gof else None),
                    n_obs=n_obs, n_var=n_var, zero=self.hist.data['Instrument Parameters'][0]['Zero'][1],
                    phases=flags, failed=bool(wR is None or not np.isfinite(wR) or wR > 95 or any(f['diverged'] for f in flags)),
                    error=(None if wR is not None and np.isfinite(wR) and wR <= 95 and not any(f['diverged'] for f in flags)
                           else 'no wR' if wR is None or not np.isfinite(wR)
                           else f'diverged (wR = {wR:.1f})' if wR > 95
                           else 'cell moved beyond its window: ' + ', '.join(
                               f"{f['name']} (> {100 * windows[n]:.0f}%)" for n, f in enumerate(flags) if f['diverged'])))

    def _n_refined(self):
        try:
            return len(self.gpx.data['Covariance']['data']['varyList'])
        except Exception:
            return None

    def export_pattern(self, path):
        """Columns: two_theta, obs, calc, background, diff (within limits)."""
        x = self.hist.getdata('x'); yo = self.hist.getdata('yobs'); yc = self.hist.getdata('ycalc')
        bg = self.hist.getdata('background')
        m = (x >= self.limits[0]) & (x <= self.limits[1])
        arr = np.column_stack([x[m], yo[m], yc[m], bg[m], (yo - yc)[m]])
        np.savetxt(path, arr, fmt='%.4f', header='two_theta obs calc background diff')
        return path

    def reflections(self):
        """{phase name: array of 2theta positions} for tick marks."""
        out = {}
        for ph in self.phases:
            try:
                refl = self.hist.reflections()[ph.name]['RefList']
                out[ph.name] = np.asarray(refl[:, 5], float)
            except Exception:
                out[ph.name] = np.zeros(0)
        return out

    def save(self):
        with quiet():
            self.gpx.save()
        shutil.copy(self.gpx_path, self.snapshot)

    @property
    def snapshot(self):
        return Path(self.gpx_path).with_name('snapshot.gpx')

    def cleanup(self):
        shutil.rmtree(self.workdir, ignore_errors=True)


def calibrate_profile(xy_file, base_instprm, cifs, out_instprm, limits=(18.0, 118.0),
                      names=None, workdir=None):
    """Fit U V W X Y SH/L and Zero on a well-identified scan; save an .instprm.

    The phases' cells, scales, fractions and background are refined first so
    the profile terms absorb only the instrumental broadening.
    """
    r = Refinement(xy_file, base_instprm, cifs, names=names, workdir=workdir, limits=limits)
    stages = []
    r.hist.set_refinements({'Sample Parameters': ['Scale']})
    stages.append(('scale+bkg', r._refine()))
    for ph in r.phases[1:]:
        ph.set_HAP_refinements({'Scale': True})
    stages.append(('fractions', r._refine()))
    # Gaussian width through W only (U and V fitted freely on one sample go
    # negative at high angle); the Lorentzian X, Y terms carry the angle
    # dependence.
    inst = r.hist.data['Instrument Parameters'][0]
    inst['U'][1] = inst['U'][0] = 0.0
    inst['V'][1] = inst['V'][0] = 0.0
    r.hist.set_refinements({'Instrument Parameters': ['W']})
    stages.append(('W', r._refine()))
    r.hist.set_refinements({'Instrument Parameters': ['X', 'Y']})
    stages.append(('XY', r._refine()))
    r.hist.set_refinements({'Instrument Parameters': ['Zero']})
    for ph in r.phases:
        ph.set_refinements({'Cell': True})
    stages.append(('zero+cells', r._refine()))
    r.hist.set_refinements({'Instrument Parameters': ['SH/L']})
    wR = r._refine()
    stages.append(('SH/L', wR))
    out_instprm = Path(out_instprm)
    out_instprm.parent.mkdir(parents=True, exist_ok=True)
    r.hist.SaveProfile(str(out_instprm))
    inst = r.hist.data['Instrument Parameters'][0]
    summary = {k: inst[k][1] for k in ('U', 'V', 'W', 'X', 'Y', 'SH/L', 'Zero')}
    summary['wR'] = wR
    summary['stages'] = stages
    summary['results'] = r.results()
    r.save()
    return summary, r
