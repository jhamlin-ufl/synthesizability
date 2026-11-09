"""Fast whole-pattern ranking of mixture hypotheses.

The observed pattern (background removed by SNIP) is treated as a
non-negative combination of simulated phase patterns.  Simulated lattice
parameters from DFT are off by up to about a percent, so each phase is first
matched on its own at a few isotropic cell scalings and the best scaling is
kept.  Each mixture hypothesis is then scored by non-negative least squares
with its (two or three) phases as columns; the fit statistic is an unweighted
Rwp-like residual on square-root intensities, which stops the strongest line
from deciding everything.  Hull excess enters as a soft prior.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from ..parsers.xrd import parse_xrd_file
from ..xrd_background import subtract_background
from ..paths import button_dir
from ..parsers.xrd import is_xrd_file


def find_xrd_files(sample_number):
    """All parseable measured patterns of one button (samples/<formula>/buttons/<NNNN>/)."""
    d = button_dir(sample_number)
    return [f for f in sorted(d.iterdir()) if f.suffix in ('.xy', '.txt') and is_xrd_file(f)] if d else []
from .library import CU_KA1
from .mixtures import atom_to_weight_fractions

CU_KA2 = 1.544426
KA2_RATIO = 0.5
STEP = 0.02
DEFAULT_RANGE = (18.0, 118.0)
BG_WINDOW = 3.0
SCALES = np.linspace(0.98, 1.02, 17)      # isotropic d-spacing scalings (DFT cells are off by up to 2 %)
FWHM = (0.0, 0.0, 0.03)                   # Caglioti U, V, W (deg^2)


def choose_scan(sample_number, which='siemens'):
    """Path of the scan to fit: Siemens powder .txt by default."""
    files = find_xrd_files(sample_number)
    if not files:
        raise FileNotFoundError(f'no XRD scans for {sample_number:04d}')
    if which == 'siemens':
        pick = [f for f in files if f.suffix == '.txt']
    elif which == 'panalytical':
        pick = [f for f in files if f.suffix == '.xy']
    else:
        pick = [f for f in files if f.name == which]
    if not pick:
        raise FileNotFoundError(f'{sample_number:04d}: no {which} scan among {[f.name for f in files]}')
    # Annealed scans are kept only when asked for by name
    plain = [f for f in pick if 'anneal' not in f.name.lower()]
    return (plain or pick)[0]


HUMP_LEVEL = 0.05     # lower envelope above this fraction of the maximum marks the low-angle hump
HUMP_WINDOW = 1.5     # degrees, running-minimum window for the envelope
HUMP_MAX_CUT = 32.0   # never trim beyond this angle


def trim_low_angle(grid, y):
    """Start of the usable range: the air-scatter and holder hump at the
    low-angle edge of many Siemens scans survives the SNIP background and is
    then fitted by dense low-symmetry cells.  The lower envelope (running
    minimum over HUMP_WINDOW) exceeds HUMP_LEVEL only where the hump
    dominates, so the range starts after the last such point."""
    from scipy.ndimage import minimum_filter1d
    n = max(3, int(HUMP_WINDOW / STEP))
    env = minimum_filter1d(y, n)
    cand = np.where((env > HUMP_LEVEL) & (grid < HUMP_MAX_CUT))[0]
    if cand.size == 0:
        return 0
    return int(cand[-1] + 1)


def load_observed(path, two_theta_range=DEFAULT_RANGE, trim=True):
    """(grid, background-subtracted intensity normalised to max 1), the
    range trimmed past the low-angle hump when ``trim`` is set."""
    d = parse_xrd_file(Path(path))
    tt, ii = np.asarray(d['two_theta'], float), np.asarray(d['intensity'], float)
    ii, _ = subtract_background(tt, ii, window_deg=BG_WINDOW)
    lo, hi = two_theta_range
    # never extend past the data: GSAS-II cannot evaluate a histogram whose
    # limits lie outside the measured range (0526, scan to 100 deg)
    lo, hi = max(lo, float(tt.min()) + STEP), min(hi, float(tt.max()) - STEP)
    grid = np.arange(lo, hi + STEP / 2, STEP)
    y = np.interp(grid, tt, ii, left=0, right=0)
    y = np.clip(y, 0, None)
    y = y / y.max()
    if trim:
        i0 = trim_low_angle(grid, y)
        if i0 > 0:
            grid, y = grid[i0:], y[i0:]
            y = y / y.max()
    return grid, y


def caglioti(two_theta, uvw=FWHM):
    u, v, w = uvw
    t = np.tan(np.radians(two_theta / 2))
    return np.sqrt(np.maximum(u * t ** 2 + v * t + w, 1e-4))


def estimate_fwhm(grid, y, n_peaks=8, default=(0.12, 0.0025)):
    """Linear FWHM(2theta) = a + b * 2theta fitted to the strongest observed
    peaks (Gaussian fits), so simulated widths follow each sample's own
    broadening.  Falls back to ``default`` when too few peaks fit."""
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks

    def g(x, a, c, s, b):
        return a * np.exp(-0.5 * ((x - c) / s) ** 2) + b

    # Below 65 deg the K-alpha doublet is not resolved, so single Gaussians
    # measure the real width; above it they would overestimate it
    pk, props = find_peaks(y, height=0.08, distance=int(0.5 / STEP))
    low = grid[pk] < 65
    pk, heights = pk[low], props['peak_heights'][low]
    order = np.argsort(heights)[::-1][:n_peaks]
    tt, fw = [], []
    for i in pk[order]:
        c = grid[i]
        w = (grid > c - 0.7) & (grid < c + 0.7)
        try:
            popt, _ = curve_fit(g, grid[w], y[w], p0=[y[i], c, 0.08, 0.0], maxfev=2000)
        except Exception:
            continue
        f = abs(popt[2]) * 2.3548
        if 0.05 < f < 1.5:
            tt.append(c)
            fw.append(f)
    if len(fw) < 3:
        return default
    tt, fw = np.array(tt), np.array(fw)
    b, a = np.polyfit(tt, fw, 1) if len(fw) >= 4 else (default[1], np.median(fw) - default[1] * np.median(tt))
    b = float(np.clip(b, 0.0, 0.003))
    a = float(np.clip(np.median(fw - b * tt), 0.08, 1.0))
    return a, b


def profile(d_hkl, intensity, grid, scale=1.0, uvw=FWHM, fwhm_ab=None):
    """K-alpha doublet Gaussian profile on ``grid`` for scaled d-spacings.

    Peak width is either the Caglioti form ``uvw`` or, when ``fwhm_ab`` is
    given, the linear ``a + b * 2theta`` estimated from the observed scan.
    Vectorised over reflections (2026-10-06): every peak's +-5 sigma window is
    evaluated in one array and accumulated with bincount, which is 20 to 50
    times faster than the former per-peak loop and numerically identical.
    """
    y = np.zeros_like(grid)
    d = np.asarray(d_hkl, float) * scale
    lo, hi = grid[0], grid[-1]
    n = grid.size
    for lam, ratio in ((CU_KA1, 1.0), (CU_KA2, KA2_RATIO)):
        s = lam / (2 * d)
        ok = s < 1
        if not ok.any():
            continue
        tt = 2 * np.degrees(np.arcsin(s[ok]))
        inten = np.asarray(intensity, float)[ok] * ratio
        fw = (fwhm_ab[0] + fwhm_ab[1] * tt) if fwhm_ab is not None else caglioti(tt, uvw)
        sig = fw / 2.3548
        keep = (tt > lo - 5 * fw) & (tt < hi + 5 * fw)
        if not keep.any():
            continue
        c, h, sg = tt[keep], inten[keep], sig[keep]
        half = int(np.ceil(5 * sg.max() / STEP)) + 1
        centre = np.round((c - lo) / STEP).astype(int)
        idx = centre[:, None] + np.arange(-half, half + 1)[None, :]
        x = lo + idx * STEP
        val = h[:, None] * np.exp(-0.5 * ((x - c[:, None]) / sg[:, None]) ** 2)
        m = (idx >= 0) & (idx < n) & (np.abs(x - c[:, None]) <= 5 * sg[:, None])
        y += np.bincount(idx[m], weights=val[m], minlength=n)[:n]
    return y


def profiles_scaled(d_hkl, intensity, grid, scales, fwhm_ab):
    """``profile`` for every scale at once: (n_scales, n_grid), one bincount."""
    scales = np.asarray(scales, float)
    n, ns = grid.size, scales.size
    lo, hi = grid[0], grid[-1]
    d0 = np.asarray(d_hkl, float)
    I0 = np.asarray(intensity, float)
    out = np.zeros(ns * n)
    for lam, ratio in ((CU_KA1, 1.0), (CU_KA2, KA2_RATIO)):
        s = lam / (2 * d0[None, :] * scales[:, None])          # (ns, nrefl)
        ok = s < 1
        if not ok.any():
            continue
        tt = np.where(ok, 2 * np.degrees(np.arcsin(np.where(ok, s, 0.5))), np.nan)
        inten = np.broadcast_to(I0[None, :] * ratio, tt.shape)
        fw = fwhm_ab[0] + fwhm_ab[1] * np.nan_to_num(tt, nan=lo)
        keep = ok & (tt > lo - 5 * fw) & (tt < hi + 5 * fw)
        if not keep.any():
            continue
        c, h, sg = tt[keep], inten[keep], (fw / 2.3548)[keep]
        row = np.broadcast_to(np.arange(ns)[:, None], tt.shape)[keep]
        half = int(np.ceil(5 * sg.max() / STEP)) + 1
        centre = np.round((c - lo) / STEP).astype(int)
        idx = centre[:, None] + np.arange(-half, half + 1)[None, :]
        x = lo + idx * STEP
        val = h[:, None] * np.exp(-0.5 * ((x - c[:, None]) / sg[:, None]) ** 2)
        m = (idx >= 0) & (idx < n) & (np.abs(x - c[:, None]) <= 5 * sg[:, None])
        flat = (row[:, None] * n + idx)[m]
        out += np.bincount(flat, weights=val[m], minlength=ns * n)[:ns * n]
    return out.reshape(ns, n)


# ---------------------------------------------------------------------------
# Per-axis cell scan: the starting cell handed to GSAS-II (2026-10-06)
# ---------------------------------------------------------------------------

def _lattice_system(spg):
    if spg >= 195:
        return 'cubic'
    if spg >= 143:
        return 'hexagonal'      # trigonal and hexagonal: a = b, c
    if spg >= 75:
        return 'tetragonal'
    if spg >= 16:
        return 'orthorhombic'
    return 'low'                # monoclinic, triclinic: isotropic only


def cell_scan(cif, grid, y_obs, fwhm_ab, spg, span=0.03, step=0.005, two_theta_range=(5, 135), min_gain=0.02):
    """Best per-axis scaling of a phase's cell against the observed pattern.

    Reflection indices and intensities come from one pymatgen pattern
    calculation; d-spacings for each candidate cell follow from the scaled
    reciprocal metric, so the scan costs one profile per candidate cell.
    Returns (scaled Structure, (sa, sb, sc), correlation).  Cubic and low-
    symmetry cells are scanned isotropically; tetragonal and hexagonal in
    (a, c); orthorhombic in (a, b, c) with a coarser step.
    """
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Structure
    s = Structure.from_file(str(cif))
    pat = XRDCalculator(wavelength=CU_KA1).get_pattern(s, two_theta_range=two_theta_range)
    if not len(pat.y):
        return s, (1.0, 1.0, 1.0), 0.0
    hkl = np.array([h[0]['hkl'][:3] for h in pat.hkls], float)
    inten = np.asarray(pat.y, float)
    inten = inten / inten.max() * 100
    yo = np.sqrt(y_obs)
    yo_n = yo / (np.linalg.norm(yo) or 1.0)
    system = _lattice_system(spg)
    f = np.arange(1 - span, 1 + span + 1e-9, step)
    if system == 'cubic' or system == 'low':
        cand = [(a, a, a) for a in f]
    elif system in ('tetragonal', 'hexagonal'):
        cand = [(a, a, c) for a in f for c in f]
    else:
        f3 = np.arange(1 - span, 1 + span + 1e-9, 2 * step)
        cand = [(a, b, c) for a in f3 for b in f3 for c in f3]
    Ginv0 = s.lattice.reciprocal_lattice_crystallographic.metric_tensor
    best = (-1.0, (1.0, 1.0, 1.0))
    r0 = None
    for sa, sb, sc in [(1.0, 1.0, 1.0)] + [c for c in cand if c != (1.0, 1.0, 1.0)]:
        S = np.diag([1 / sa, 1 / sb, 1 / sc])        # scaling real axes by s scales the reciprocal metric by 1/s
        Ginv = S @ Ginv0 @ S
        d = 1 / np.sqrt(np.einsum('ni,ij,nj->n', hkl, Ginv, hkl))
        prof = np.sqrt(profile(d, inten, grid, 1.0, fwhm_ab=fwhm_ab))
        nrm = np.linalg.norm(prof)
        r = float(prof @ yo_n / nrm) if nrm > 0 else 0.0
        if r0 is None:
            r0 = r                       # the unscaled cell
        if r > best[0]:
            best = (r, (sa, sb, sc))
    r, (sa, sb, sc) = best
    # A phase that is absent or minor latches onto other phases' lines and
    # drifts to the edge of the window.  Keep the scanned cell only when it
    # clearly beats the unscaled one and stays within 2 %; otherwise start
    # GSAS-II from the DFT cell as before.
    if r0 is not None and (r - r0 < min_gain or max(abs(v - 1) for v in (sa, sb, sc)) > 0.02 + 1e-9):
        r, (sa, sb, sc) = r0, (1.0, 1.0, 1.0)
    lat = s.lattice
    new = lat.__class__.from_parameters(lat.a * sa, lat.b * sb, lat.c * sc, lat.alpha, lat.beta, lat.gamma)
    scaled = Structure(new, s.species_and_occu, s.frac_coords)
    return scaled, (sa, sb, sc), r


def _corr(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0


# Per-phase width factors on the sample width.  Tested 2026-10-06 on 0470 and
# 0512 with (1.0, 1.6, 2.5): the correlation criterion prefers the broadest
# profile for almost every phase (broad profiles correlate better with any
# pattern), the single bcc solution of 0470 only improved from 0.58 to 0.51
# while the bcc-derived mixtures kept winning, and 0512 did not change.  Kept
# at a single width until the choice is made on the NNLS residual rather than
# on the correlation.
WIDTHS = (1.0,)


def phase_patterns(phases, grid, y_obs, scales=SCALES, fwhm_ab=None, widths=WIDTHS):
    """Best-scaled sqrt profile of every phase and its single-phase score.

    Each phase is tried at every cell scaling and, at its best scaling, at
    every width factor in ``widths`` (2026-10-05: one width per sample made
    the ranking blind to a broadened solid solution next to sharp lines, 0470
    and 0512).  Returns (P (n_phases, n_grid) sqrt-intensity profiles
    normalised to unit norm, DataFrame with ``scale``, ``width`` and ``corr``
    per phase).
    """
    if fwhm_ab is None:
        fwhm_ab = estimate_fwhm(grid, y_obs)
    yo = np.sqrt(y_obs)
    P = np.zeros((len(phases), grid.size))
    rows = []
    scales = np.asarray(scales, float)
    yo_n = yo / (np.linalg.norm(yo) or 1.0)
    for n, ph in enumerate(phases):
        if ph.d_hkl.size == 0:
            rows.append(dict(key=ph.key, formula=ph.formula, source=ph.source, spg=ph.spg,
                             disordered=ph.disordered, at_target=ph.at_target,
                             e_form=ph.e_form, stab=ph.stab, scale=1.0, width=1.0, corr=0.0))
            continue
        # coarse-to-fine: every fourth scale first, then the neighbours of the best
        coarse = scales[::4] if scales.size > 8 else scales
        profs = np.sqrt(profiles_scaled(ph.d_hkl, ph.intensity, grid, coarse, fwhm_ab))
        norms = np.linalg.norm(profs, axis=1)
        corrs = np.where(norms > 0, profs @ yo_n / np.where(norms > 0, norms, 1), 0.0)
        ib = int(np.argmax(corrs))
        r, sc, prof = float(corrs[ib]), float(coarse[ib]), profs[ib]
        if scales.size > 8:
            step = scales[1] - scales[0]
            fine = scales[(np.abs(scales - sc) < 4 * step - 1e-9) & (np.abs(scales - sc) > 1e-9)]
            if fine.size:
                profs2 = np.sqrt(profiles_scaled(ph.d_hkl, ph.intensity, grid, fine, fwhm_ab))
                norms2 = np.linalg.norm(profs2, axis=1)
                corrs2 = np.where(norms2 > 0, profs2 @ yo_n / np.where(norms2 > 0, norms2, 1), 0.0)
                jb = int(np.argmax(corrs2))
                if corrs2[jb] > r:
                    r, sc, prof = float(corrs2[jb]), float(fine[jb]), profs2[jb]
        wf = 1.0
        for w in widths:
            if w == 1.0:
                continue
            wide = np.sqrt(profile(ph.d_hkl, ph.intensity, grid, sc, fwhm_ab=(fwhm_ab[0] * w, fwhm_ab[1] * w)))
            rw = _corr(wide, yo)
            if rw > r + 0.01:
                r, prof, wf = rw, wide, w
        if prof is not None and prof.max() > 0:
            P[n] = prof / np.linalg.norm(prof)
        rows.append(dict(key=ph.key, formula=ph.formula, source=ph.source, spg=ph.spg,
                         disordered=ph.disordered, at_target=ph.at_target,
                         e_form=ph.e_form, stab=ph.stab, scale=sc, width=wf, corr=r))
    table = pd.DataFrame(rows)
    table.attrs['fwhm_ab'] = fwhm_ab
    return P, table


def batched_nnls(G, b, yy, combos):
    """Non-negative least squares of y on up to three unit-norm columns, for
    many column triples at once.

    ``G`` is the Gram matrix of all profiles, ``b`` their projections on the
    observed pattern and ``yy`` its squared norm.  ``combos`` is an (n, 3)
    integer array (-1 pads shorter combinations).  With at most three columns
    the KKT conditions can be checked by solving the normal equations on every
    non-empty subset of columns and keeping the feasible (all coefficients
    non-negative) subset with the smallest residual; the unconstrained
    solution of the active set is optimal when the remaining gradients are
    non-negative, and the best feasible subset over all subsets satisfies
    that.  Returns (coefficients (n, 3), squared residuals (n,)).
    """
    n = len(combos)
    valid = combos >= 0
    idx = np.where(valid, combos, 0)
    Gc = G[idx[:, :, None], idx[:, None, :]]            # (n, 3, 3)
    bc = b[idx]                                          # (n, 3)
    best_res = np.full(n, np.inf)
    best_coef = np.zeros((n, 3))
    for mask in range(1, 8):
        cols = [t for t in range(3) if mask >> t & 1]
        ok = valid[:, cols].all(axis=1)
        if not ok.any():
            continue
        sub = np.ix_(np.where(ok)[0], cols, cols)
        A = Gc[sub]
        rhs = bc[np.where(ok)[0]][:, cols]
        A = A + 1e-10 * np.eye(len(cols))
        try:
            x = np.linalg.solve(A, rhs[..., None])[..., 0]
        except np.linalg.LinAlgError:
            continue
        feas = (x >= -1e-12).all(axis=1)
        if not feas.any():
            continue
        res = yy - np.einsum('ij,ij->i', x, rhs)         # ||y||^2 - x.b for LS optimum
        rows = np.where(ok)[0][feas]
        better = res[feas] < best_res[rows]
        rows = rows[better]
        best_res[rows] = res[feas][better]
        coef = np.zeros((len(rows), 3))
        coef[:, cols] = x[feas][better]
        best_coef[rows] = coef
    return best_coef, np.maximum(best_res, 0.0)


TRACE_SHARE = 0.03        # a fitted phase below this share of the intensity is "trace"
TRACE_PENALTY = 0.05      # per trace phase
PHASE_PENALTY = 0.02      # per phase beyond the first (Occam)
OFF_NOMINAL_PENALTY = 0.1  # times the barycentric overshoot
REFL_REF = 150            # reflections in range above which a phase is "dense"
REFL_PENALTY = 0.03       # per decade of reflections above REFL_REF (Occam)


PAIR_TOP = 5                                  # polymorphs of the target paired in the ranking (2026-10-06)
PAIR_RATIOS = np.linspace(0.1, 0.9, 9)        # fixed-mode split of the target fraction between the two


def score_mixtures(mix, members, phases, P, y_obs, prior_scale=1.0, top_polymorphs=2,
                   phase_table=None, chunk=200_000, keep=100_000, pairs=True):
    """NNLS fit of every mixture hypothesis, vectorised over hypotheses.

    Each composition point of a mixture may be realised by any of its
    polymorphs; the ``top_polymorphs`` best single-phase matches at that point
    are tried, except at the target composition, where every polymorph is
    tried (the predicted structure and its competitors are the point of the
    exercise).  Returns one row per phase combination with ``rwp`` (free
    NNLS, on sqrt intensities), ``rwp_fixed`` (one scale factor, phases in
    their lever-rule weight ratio), intensity shares, the lever-rule
    fractions and the composite score

        rwp + prior_scale * excess_eV + 0.02 (n_phases - 1)
            + 0.05 n_trace + 0.1 off_nominal
    """
    import itertools
    yo = np.sqrt(y_obs)
    yy = float(yo @ yo)
    corr = phase_table['corr'].to_numpy() if phase_table is not None else None

    def choices(point):
        idx = list(members[point])
        if corr is not None:
            idx.sort(key=lambda t: -corr[t])
        if any(phases[t].at_target for t in idx):
            return idx
        return idx[:top_polymorphs]

    # Expand mixtures into phase combinations (int32 arrays, not tuples:
    # large spaces reach 10^7 combinations)
    combo_blocks, row_blocks = [], []
    for n, row in enumerate(mix.itertuples(index=False)):
        pts = [int(t) for t in (row.i, row.j, row.k) if t >= 0]
        ch = [choices(p) for p in pts]
        block = np.array(list(itertools.product(*ch)), dtype=np.int32).reshape(-1, len(pts))
        if block.shape[1] < 3:
            block = np.hstack([block, np.full((len(block), 3 - block.shape[1]), -1, dtype=np.int32)])
        combo_blocks.append(block)
        row_blocks.append(np.full(len(block), n, dtype=np.int32))
    # Polymorph pairs at the target composition (2026-10-06): two polymorphs
    # of the target in one hypothesis, alone or with one partner, so a
    # partially transformed or two-polymorph sample can be proposed.  The free
    # NNLS solves their ratio; for the fixed mode the target fraction is split
    # over PAIR_RATIOS.  Each pair row is appended to ``mix`` so that the
    # fraction, excess and prior bookkeeping below needs no special case.
    if pairs:
        tgt = [p for p in range(len(members)) if any(phases[t].at_target for t in members[p])]
        if tgt:
            p_t = tgt[0]
            idx = list(members[p_t])
            if corr is not None:
                idx.sort(key=lambda t: -corr[t])
            prs = list(itertools.combinations(idx[:PAIR_TOP], 2))
            new_rows, new_combos = [], []
            cols = list(mix.columns)
            for n, row in enumerate(mix.itertuples(index=False)):
                pts = [int(t) for t in (row.i, row.j, row.k) if t >= 0]
                if p_t not in pts or len(pts) > 2:
                    continue
                fr_row = [row.f1, row.f2, row.f3]
                f_t = fr_row[pts.index(p_t)]
                others = [p for p in pts if p != p_t]
                partner_choices = [choices(others[0])[:top_polymorphs]] if others else [[-1]]
                f_o = fr_row[pts.index(others[0])] if others else 0.0
                base = row._asdict()
                for i, j in prs:
                    for r in PAIR_RATIOS:
                        for t_o in partner_choices[0]:
                            rec = dict(base)
                            rec.update(f1=r * f_t, f2=(1 - r) * f_t, f3=f_o)
                            new_rows.append([rec.get(c) for c in cols])
                            new_combos.append([i, j, t_o])
            if new_rows:
                mix = pd.concat([mix, pd.DataFrame(new_rows, columns=cols)], ignore_index=True)
                combo_blocks.append(np.array(new_combos, dtype=np.int32))
                row_blocks.append(np.arange(len(mix) - len(new_rows), len(mix), dtype=np.int32))
    combos_all = np.vstack(combo_blocks)
    mix_row_all = np.concatenate(row_blocks)
    del combo_blocks, row_blocks

    G = P @ P.T
    b = P @ yo
    masses = np.array([p.composition.weight / p.composition.num_atoms for p in phases])
    n_refl = np.array([p.d_hkl.size for p in phases])
    fr_all = mix[['f1', 'f2', 'f3']].to_numpy()
    excess_all = mix['excess'].to_numpy()
    off_all = mix['off_nominal'].to_numpy() if 'off_nominal' in mix else np.zeros(len(mix))

    # Score in chunks and keep only the best ``keep`` rows: the (n, 3, 3)
    # temporaries of the batched solve are what exhausts memory otherwise
    kept = None
    for start in range(0, len(combos_all), chunk):
        combos = combos_all[start:start + chunk].astype(int)
        mix_row = mix_row_all[start:start + chunk]
        coef, res2 = batched_nnls(G, b, yy, combos)
        rwp = np.sqrt(res2) / np.sqrt(yy)
        valid = combos >= 0
        cidx = np.where(valid, combos, 0)
        fr = fr_all[mix_row]
        w = np.where(valid, fr * masses[cidx], 0.0)
        w = w / w.sum(axis=1, keepdims=True)
        share = coef ** 2
        ssum = share.sum(1, keepdims=True)
        share = np.where(valid, share / np.where(ssum > 0, ssum, 1), 0.0)
        # Fixed-ratio fit: one column, phases weighted by lever-rule weight
        # fraction (intensity taken proportional to weight), one free scale
        Gc = G[cidx[:, :, None], cidx[:, None, :]]
        vv = np.einsum('ni,nij,nj->n', w, Gc, w)
        vb = np.einsum('ni,ni->n', w, b[cidx])
        del Gc
        res_fixed = yy - np.where((vb > 0) & (vv > 0), vb ** 2 / np.where(vv > 0, vv, 1), 0.0)
        rwp_fixed = np.sqrt(np.maximum(res_fixed, 0)) / np.sqrt(yy)
        n_trace = ((share < TRACE_SHARE) & valid).sum(1)
        n_ph = valid.sum(1)
        excess = excess_all[mix_row]
        off = off_all[mix_row]
        # Dense low-symmetry cells (hundreds of lines) fit anything; penalise them
        dense = np.where(valid, np.log10(np.maximum(n_refl[cidx], 1) / REFL_REF).clip(0, None), 0.0).sum(1)
        score = (rwp + prior_scale * excess + PHASE_PENALTY * (n_ph - 1)
                 + TRACE_PENALTY * n_trace + OFF_NOMINAL_PENALTY * off + REFL_PENALTY * dense)
        block = dict(combos=combos, mix_row=mix_row, coef=coef, rwp=rwp, rwp_fixed=rwp_fixed,
                     share=share, w=w, n_trace=n_trace, score=score)
        if kept is None:
            kept = block
        else:
            kept = {k: np.concatenate([kept[k], block[k]]) for k in block}
        if len(kept['score']) > keep:
            # the top ``keep`` by score, plus every single-phase hypothesis
            # (there are only as many as phases at the target composition,
            # and a solid-solution cell must never fall out of the kept set)
            sel = np.argpartition(kept['score'], keep)[:keep]
            single = np.where((kept['combos'] >= 0).sum(1) == 1)[0]
            sel = np.union1d(sel, single)
            kept = {k: v[sel] for k, v in kept.items()}
    n_total = len(combos_all)
    del combos_all, mix_row_all

    combos, mix_row, coef = kept['combos'], kept['mix_row'], kept['coef']
    rwp, rwp_fixed, share, w = kept['rwp'], kept['rwp_fixed'], kept['share'], kept['w']
    n_trace, score = kept['n_trace'], kept['score']
    valid = combos >= 0
    cidx = np.where(valid, combos, 0)
    n_ph = valid.sum(1)
    fr = fr_all[mix_row]
    excess = excess_all[mix_row]
    off = off_all[mix_row]
    m = mix.iloc[mix_row].reset_index(drop=True)
    keys = np.array([p.key for p in phases])
    formulas = np.array([p.formula for p in phases])

    def join(arr, fmt=None):
        out = []
        for row, v in zip(arr, valid):
            vals = row[v]
            out.append(' / '.join(fmt % t for t in vals) if fmt else ' + '.join(vals))
        return out

    df = pd.DataFrame(dict(
        phases=join(keys[cidx].astype(object)),
        formulas=join(formulas[cidx].astype(object)),
        n_phases=n_ph, excess_meV=excess * 1000, rwp=rwp, rwp_fixed=rwp_fixed,
        score=score, off_nominal=off, n_trace=n_trace,
        n_refl=np.where(valid, n_refl[cidx], 0).sum(1),
        lever_wt=join(w, '%.2f'), fit_share=join(share, '%.2f'),
        zero_phases=((coef <= 1e-6) & valid).sum(1),
        prior_source=m['prior_source'].to_numpy() if 'prior_source' in m else '',
        keys=[tuple(int(t) for t in c if t >= 0) for c in combos],
        fracs=[tuple(float(t) for t in f[:k]) for f, k in zip(fr, n_ph)]))
    df = df.sort_values('score').reset_index(drop=True)
    df.attrs['n_total'] = n_total
    return df


def unexplained_peaks(grid, y_obs, P, combo, min_height=0.03):
    """Residual maxima after the NNLS fit of ``combo``: (two_theta, height)."""
    yo = np.sqrt(y_obs)
    coef, _ = nnls(P[list(combo)].T, yo)
    resid = yo - P[list(combo)].T @ coef
    from scipy.signal import find_peaks
    pk, props = find_peaks(resid, height=min_height * yo.max(), distance=int(0.3 / STEP))
    return list(zip(grid[pk].round(2), (props['peak_heights'] ** 2).round(3)))
