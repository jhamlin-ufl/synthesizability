"""Shared machinery for the measured-vs-predicted XRD figures.

Two scripts draw these panels -- ``scripts/plot_xrd_representative.py`` (the
six-panel figure in the body of the paper) and
``scripts/plot_xrd_comparison_grid.py`` (the per-outcome grids in the
appendix) -- and they must agree on every choice that affects what a panel
shows: which measurement file is plotted, how the background is removed, how
the predicted structure is simulated, which structural solution is overlaid.
All of that lives here so it is decided once.

Each panel carries, from bottom to top:

  red    the pattern simulated from the most thermodynamically stable OQMD
         structure at that composition -- i.e. what was predicted
  blue   the measured pattern, after background subtraction
  green  the pattern simulated from our structural solution, where we have one
         (``SOLUTION`` below)

Page geometry is here too, so that the body figure and the appendix figures
have the same width and the same per-row panel height.
"""
from __future__ import annotations

import string
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MultipleLocator
from pymatgen.core import Composition, Lattice, Structure

from synthesizability.formula import formula_to_mathtext
from synthesizability.parsers.xrd import is_xrd_file, parse_xrd_file
from synthesizability.xrd_background import subtract_background
from synthesizability.xrd_simulator import simulate_pattern

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
RAW_DIR = Path('data/raw')
OQMD_DIR = Path('data/external/oqmd_structures')

# Plotted 2-theta window.  The lower limit matters: the clearest single piece of
# evidence for disorder in this set is a predicted superlattice reflection in the
# 20-30 deg region that is simply absent from the measurement, so the window has
# to reach down there.  Below ~20 deg the Siemens patterns are almost entirely
# the amorphous hump from the slide and the grease.
TWO_THETA_RANGE = (20.0, 120.0)

FWHM = 0.4              # deg, broadening applied to the simulated patterns
BG_WINDOW = 4.0         # deg, SNIP clipping half-width (see xrd_background)

# The background is estimated over TWO_THETA_RANGE widened by this much, then
# the result is clipped back.  The padding keeps SNIP's edge behaviour outside
# the plotted window.  It is deliberately small: on the Siemens scans the
# amorphous hump from the slide and grease peaks near 12 deg and is wide enough
# that SNIP reads it as a peak rather than as background, so including it in the
# fit range leaves its whole flank behind after subtraction.
BG_FIT_PAD = 2.0        # deg

# Samples with more than one usable measurement.  The default rule picks the
# file covering the most of TWO_THETA_RANGE, breaking ties on point count, which
# in practice means the 5-120 deg Siemens scan over a shorter Panalytical one.
# List a filename here to overrule that.  Sample 0483's Panalytical file is
# mislabelled "La5Al3Cu2" but is the La5Al3Co2 measurement.
# 0466 Hf2MoIr, 0473 Sc3CuPd2, 0474 Sc2CuPd and 0539 CoSi3Y were annealed (1000 C,
# one week, sealed under Ar) and measured again afterwards.  Per the 2026-09
# review the paper shows the annealed scan for 0466, 0473 and 0539 (cleaner, and
# the low-angle argument holds in both states) and the as-cast bulk scan for
# 0474 (the annealed scan, on ground material, lost most of its peaks).
FILE_OVERRIDE = {
    466: '0466_HM_Hf2MoIr_annealed.txt',
    473: '0473_HM_Sc3CuPd2_annealed.txt',
    474: '0474_HM_Sc2CuPd.xy',
    539: '0539_HM_CoSi3Y_annealed.txt',
}


def b2(a, site_a, site_b):
    """CsCl-type (B2) cell: ``site_a`` at the origin, ``site_b`` at the body
    centre, each a species or an occupancy dict such as {'Fe': 0.25, 'Co': 0.75}.
    With the same occupancy on both sites it is a disordered bcc (A2) cell."""
    return Structure(Lattice.cubic(a), [site_a, site_b],
                     [[0, 0, 0], [0.5, 0.5, 0.5]])


def disordered_bcc(formula, a):
    """A2 (bcc, one mixed-occupancy site) structure with lattice parameter ``a``."""
    occupancy = dict(Composition(formula).fractional_composition.as_dict())
    return b2(a, occupancy, occupancy)


def _interp_occ(c, full, eta):
    return {e: v for e, v in {e: c.get(e, 0) + eta * (full.get(e, 0) - c.get(e, 0))
                             for e in c}.items() if v > 1e-9}


def _complement(c, site, weight):
    return {e: v for e, v in {e: c[e] - weight * (site.get(e, 0) - c[e])
                             for e in c}.items() if v > 1e-9}


def b2_partial(a, formula, ordered_element, site_fraction):
    """B2 with ``ordered_element`` occupying ``site_fraction`` of its own site
    and the rest of the composition spread in nominal ratio (anti-site
    disorder).  site_fraction = 1 is full B2 order, 0.5 is A2."""
    c = dict(Composition(formula).fractional_composition.as_dict())
    eta = (site_fraction - c[ordered_element]) / (1 - c[ordered_element])
    A = _interp_occ(c, {ordered_element: 1.0}, eta)
    return b2(a, A, _complement(c, A, 1.0))


def fcc(a, occupancy):
    """Disordered fcc with one mixed site."""
    return Structure(Lattice.cubic(a), [occupancy] * 4,
                     [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])


def l12_partial(a, formula, corner_element, corner_fraction):
    """Cu3Au-type (L1_2) with ``corner_element`` occupying ``corner_fraction``
    of the corner site; the faces take the remainder in nominal ratio."""
    c = dict(Composition(formula).fractional_composition.as_dict())
    eta = (corner_fraction - c[corner_element]) / (1 - c[corner_element])
    corner = _interp_occ(c, {corner_element: 1.0}, eta)
    face = _complement(c, corner, 1 / 3)
    return Structure(Lattice.cubic(a), [corner, face, face, face],
                     [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])


def c11b_au2_almn(a=3.333, c=8.897):
    """MoSi2-type Au2(Al2/3 Mn1/3), on the Au2Al cell of PDF 04-003-0365."""
    template = Structure.from_file(RAW_DIR / '0479_HM_MnAl2Au6' / 'PDF Card - 04-003-0365.cif')
    out = Structure(Lattice.tetragonal(a, c), [], [])
    for site in template:
        out.append({'Al': 2 / 3, 'Mn': 1 / 3} if site.specie.symbol == 'Al' else 'Au',
                   site.frac_coords)
    return out


class Jade:
    """Use the Jade whole-pattern fit (calculated minus background) as the
    solution trace.  ``filename`` is the ``*_XRD_fit.txt`` export in the
    sample's raw folder; None means ``<sample>_XRD_fit.txt``."""

    def __init__(self, filename=None):
        self.filename = filename


# The phase or mixture we believe actually formed, keyed by sample number, as
# (spec, label).  ``spec`` is a pymatgen Structure (a model we simulate), a CIF
# path, or ``Jade(...)`` (the students' whole-pattern fit, drawn as exported).
# ``label`` is printed on the panel.  Samples absent from this table get no
# green trace: no structural solution has been identified for them.
#
# Site fractions in the disordered models come from the integrated-intensity
# analysis of 2026-09-30 (data/temp/occupancy_goodness.py); lattice constants
# from least-squares fits to the measured peak positions of the scan shown.
SOLUTION: dict[int, tuple] = {
    # --- solid solutions: simulated models --------------------------------
    # 0466 Hf2MoIr, annealed: fully disordered bcc; the L2_1 (111) at 23.5 deg
    # and the B2 (100) at 27.4 deg are both absent.  (A little Mo is also
    # present; not included in the trace.)
    466: (disordered_bcc('Hf2MoIr', 3.258), 'A2 bcc (Hf,Mo,Ir)'),
    # 0473 Sc3CuPd2, annealed: B2, Sc site 0.94 Sc.  The predicted P4/mmm
    # line at 8.9 deg is absent.
    473: (b2_partial(3.277, 'Sc3CuPd2', 'Sc', 0.94), 'B2 Sc(Cu,Pd), Sc site 0.94'),
    # 0474 Sc2CuPd, as-cast bulk: B2, Sc site 0.89; the L2_1 (111) at 23.6 deg
    # (predicted 5.5 %) is absent below 0.7 %.
    474: (b2_partial(3.271, 'Sc2CuPd', 'Sc', 0.89), 'B2 Sc(Cu,Pd), Sc site 0.89'),
    # 0477 ScPd5Au2: fcc solid solution; orderings indistinguishable on the
    # 30-130 deg scan.  Interim until the remake.
    477: (fcc(3.968, {'Sc': 1 / 8, 'Pd': 5 / 8, 'Au': 2 / 8}), 'fcc (Sc,Pd,Au)'),
    # 0478 Sc2Pd5Au: L1_2, corner site 0.77 Sc, faces 0.19 Sc.
    478: (l12_partial(3.980, 'Sc2Pd5Au', 'Sc', 0.77), r'L1$_2$ Sc(Pd,Au)$_3$, corner 0.77 Sc'),
    # 0479 MnAl2Au6: MoSi2-type Au2(Al,Mn); the predicted line at 19.4 deg is
    # absent.  Peaks are broad; site order not refined.
    479: (c11b_au2_almn(), r'MoSi$_2$-type Au$_2$(Al,Mn)'),
    # 0481 NiGe2Pt: fcc (Ni,Pt,Ge); no evidence for Ge ordering.  A Pt
    # germanide second phase is not included.
    481: (fcc(3.587, {'Ni': .25, 'Pt': .25, 'Ge': .5}), 'fcc (Ni,Pt,Ge)'),
    # 0535, 0545: B2 Al(Fe,Co), Al site 0.92 Al.
    535: (b2_partial(2.864, 'Al8FeCo7', 'Al', 0.92), 'B2 Al(Fe,Co), Al site 0.92'),
    545: (b2_partial(2.868, 'Al4FeCo3', 'Al', 0.92), 'B2 Al(Fe,Co), Al site 0.92'),
    # 0556 Gd3CoSn6: the known CeNiSi2-type GdCo0.27Sn2 (ordered child =
    # the prediction) plus CoSn3, as fitted.
    556: (Jade(), r'Whole-pattern fit: GdCo$_{0.27}$Sn$_2$ + CoSn$_3$'),
    # --- multiphase: the Jade whole-pattern fits --------------------------
    483: (Jade(), r'Whole-pattern fit: La$_4$Co$_3$Al$_3$ only'),
    494: (Jade(), 'Whole-pattern fit: TiSi + Pt-rich fcc'),
    514: (Jade(), 'Whole-pattern fit: four La-Ni-Si phases',
          'Whole-pattern fit:\nLa$_3$Ni$_3$Si$_7$ 52 wt%\nLaNiSi$_2$ 31 wt%\nSi 12 wt%\nLaNi$_6$Si$_6$ 6 wt%'),
    537: (Jade(), 'Whole-pattern fit: four Gd-Co-Sn phases'),
    538: (Jade(), r'Whole-pattern fit: $Pnma$ YCoSi + YCo$_2$Si$_2$ + Si',
          'Whole-pattern fit:\n$Pnma$ YCoSi 89 wt%\nYCo$_2$Si$_2$ 9 wt%\nSi 2 wt%'),
    540: (Jade(), 'Whole-pattern fit: three Y-Al-Cu phases'),
    553: (Jade(), 'Whole-pattern fit: two Gd-Fe-Ge phases + Ge'),
    554: (Jade(), 'Whole-pattern fit: three Gd-Fe-Ge phases + Ge'),
    555: (Jade(), 'Whole-pattern fit: three Gd-Fe-Ge phases'),
    # --- predicted structure found ----------------------------------------
    507: (Jade(), r'Whole-pattern fit: YCuSn + Y$_3$Cu$_4$Sn$_4$',
          'Whole-pattern fit:\nYCuSn (PDF 04-014-1199) 91 wt%\nY$_3$Cu$_4$Sn$_4$ 9 wt%'),
    # 0539 annealed: Cmmm YCoSi3 (+8 meV) ~80 wt%; the second phase of that
    # fit is not physical and the 28.1 deg line is unexplained.
    539: (Jade('0539_annealed_XRD_fit.txt'), r'Whole-pattern fit: $Cmmm$ YCoSi$_3$ (+8 meV), 80 wt%'),
}

# ---------------------------------------------------------------------------
# Page geometry
# ---------------------------------------------------------------------------
# Figures are sized to the A4 text block, less the space a caption will occupy
# underneath.  Springer Nature sets npj figure legends in 8 pt sans;
# CAPTION_LEADING is the line spacing.  Nothing is reserved *inside* the figure
# -- the caption goes below it in the manuscript -- so savefig must not be given
# bbox_inches='tight', which would change the size computed here.
CM = 1 / 2.54                 # cm -> inch
A4_CM = (21.0, 29.7)
MARGIN_CM = 1.4
CAPTION_PT = 8.0
CAPTION_LEADING = 1.15
CAPTION_GAP_CM = 0.4          # clear space between the axes block and the caption

FIG_W_CM = A4_CM[0] - 2 * MARGIN_CM


def page_height_cm(caption_lines):
    """Height left on an A4 page for the figure above a caption of N lines."""
    caption_cm = caption_lines * CAPTION_PT * CAPTION_LEADING / 72 / CM
    return A4_CM[1] - 2 * MARGIN_CM - caption_cm - CAPTION_GAP_CM


# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------
COLOR_MEASURED = '#1565C0'
COLOR_PREDICTED = '#D32F2F'
COLOR_SOLUTION = '#2E7D32'

LW_MEASURED = 0.55
LW_SIMULATED = 0.55

GAP = 0.12              # normalised units of clear space between stacked traces

XLABEL = r'2$\theta$ (deg)'
YLABEL = r'$I$ (arb. units)'
XTICK_MAJOR = 20              # deg; the edge ticks (20 and 120) are drawn but not
                              # labelled, so neighbouring columns' labels never collide
XTICK_MINOR = 10              # deg, default minor spacing (large panels use 5)

LEGEND_PREDICTED = 'OQMD predicted'
LEGEND_MEASURED = 'Measured (background subtracted)'
LEGEND_SOLUTION = 'Model or whole-pattern fit'

# ---------------------------------------------------------------------------
# Data location
# ---------------------------------------------------------------------------


def find_xrd_files(sample_number):
    """All parseable measured patterns for one sample."""
    dirs = [d for d in RAW_DIR.iterdir()
            if d.is_dir() and d.name.startswith(f'{sample_number:04d}_')]
    if not dirs:
        return []
    return [f for f in sorted(dirs[0].iterdir())
            if f.suffix in ('.xy', '.txt') and is_xrd_file(f)]


def choose_xrd_file(sample_number):
    """Pick the measurement to plot: the widest coverage of TWO_THETA_RANGE."""
    candidates = find_xrd_files(sample_number)
    if not candidates:
        return None

    override = FILE_OVERRIDE.get(sample_number)
    if override is not None:
        match = [f for f in candidates if f.name == override]
        if not match:
            raise FileNotFoundError(
                f'FILE_OVERRIDE for {sample_number:04d} names {override!r}, '
                f'which is not among {[f.name for f in candidates]}')
        return match[0]

    lo, hi = TWO_THETA_RANGE
    best, best_key = None, None
    for path in candidates:
        try:
            d = parse_xrd_file(path)
        except Exception:
            continue
        tt = d['two_theta']
        if tt.size < 10:
            continue
        coverage = max(0.0, min(tt.max(), hi) - max(tt.min(), lo))
        key = (round(coverage, 1), tt.size)
        if best_key is None or key > best_key:
            best, best_key = path, key
    return best


def find_oqmd_cif(formula, entry_id):
    """CIF for the most stable OQMD structure at this composition."""
    cif_dir = OQMD_DIR / formula
    if not cif_dir.exists():
        return None
    if pd.notna(entry_id):
        cif = cif_dir / f'{int(entry_id)}.cif'
        if cif.exists():
            return cif
    cifs = sorted(cif_dir.glob('*.cif'))
    return cifs[0] if cifs else None


def resolve_solution(sample_number):
    """(spec, label, detail) for the structural solution, or (None, None, None).
    ``detail`` is an optional longer label for large panels.  ``spec`` is a
    Structure, an existing CIF path, or a ``Jade`` marker."""
    entry = SOLUTION.get(sample_number)
    if entry is None:
        return None, None, None
    spec, label, *rest = entry
    detail = rest[0] if rest else None
    if not isinstance(spec, (Structure, Jade)):
        spec = Path(spec)
        if not spec.exists():
            raise FileNotFoundError(
                f'SOLUTION for {sample_number:04d} points at {spec}, which does not exist')
    return spec, label, detail


def load_jade_fit(sample_number, filename=None):
    """(two_theta, calculated - background) from a Jade ``*_XRD_fit.txt`` export,
    clipped to TWO_THETA_RANGE and normalised."""
    folder = [d for d in RAW_DIR.iterdir() if d.name.startswith(f'{sample_number:04d}_')][0]
    path = folder / (filename or f'{sample_number:04d}_XRD_fit.txt')
    if not path.exists():
        raise FileNotFoundError(f'Jade export {path} not found')
    rows = []
    for line in path.read_text(encoding='latin-1').splitlines():
        parts = line.split('\t')
        if len(parts) < 5:
            continue
        try:
            rows.append([float(v) for v in parts[:5]])
        except ValueError:
            continue
    data = np.array(rows)
    tt, calc, bg = data[:, 0], data[:, 2], data[:, 4]
    lo, hi = TWO_THETA_RANGE
    m = (tt >= lo) & (tt <= hi) & (calc >= 0)
    return tt[m], normalise(np.clip(calc[m] - bg[m], 0, None)), path


# ---------------------------------------------------------------------------
# Pattern preparation
# ---------------------------------------------------------------------------


def normalise(ii):
    """Scale to a maximum of 1 with the baseline at 0."""
    lo, hi = float(np.min(ii)), float(np.max(ii))
    return (ii - lo) / (hi - lo) if hi > lo else np.zeros_like(ii)


def load_measured(sample_number):
    """Background-subtracted measured pattern, clipped to TWO_THETA_RANGE."""
    path = choose_xrd_file(sample_number)
    if path is None:
        return None, None
    data = parse_xrd_file(path)
    tt, ii = data['two_theta'], data['intensity'].astype(float)

    lo, hi = TWO_THETA_RANGE
    fit = (tt >= lo - BG_FIT_PAD) & (tt <= hi + BG_FIT_PAD)
    tt, ii = tt[fit], ii[fit]
    if tt.size < 10:
        return None, path
    ii, _ = subtract_background(tt, ii, window_deg=BG_WINDOW)

    mask = (tt >= lo) & (tt <= hi)
    if mask.sum() < 10:
        return None, path
    return (tt[mask], normalise(ii[mask])), path


def simulate(structure):
    """Simulated pattern for a CIF path or Structure; None if it cannot be done."""
    if structure is None:
        return None
    try:
        return simulate_pattern(structure, two_theta_range=TWO_THETA_RANGE,
                                fwhm=FWHM)
    except Exception as exc:
        print(f'  WARNING: simulation failed for {structure}: {exc}')
        return None


def space_group_symbol(cif_path):
    """Hermann-Mauguin symbol of the structure in a CIF, or None."""
    if cif_path is None:
        return None
    try:
        return Structure.from_file(cif_path).get_space_group_info()[0]
    except Exception:
        return None


def build_record(row):
    """Everything one panel needs, for one row of ``load_classified_samples``.

    Returns None (after printing why) if the sample has no usable measurement.
    """
    num = int(row['sample_number'])
    formula = row['formula']

    measured, path = load_measured(num)
    if measured is None:
        print(f'  SKIP {num:04d} {formula}: no usable measurement')
        return None

    oqmd_cif = find_oqmd_cif(formula, row['oqmd_entry_id'])
    predicted = simulate(oqmd_cif)
    solution_spec, solution_label, solution_detail = resolve_solution(num)
    if isinstance(solution_spec, Jade):
        jt, ji, jpath = load_jade_fit(num, solution_spec.filename)
        solution = (jt, ji)
        # The fit must belong to the scan we show: compare angular ranges.
        if abs(jt.min() - max(measured[0].min(), TWO_THETA_RANGE[0])) > 0.6 or \
                abs(jt.max() - min(measured[0].max(), TWO_THETA_RANGE[1])) > 0.6:
            print(f'  WARNING {num:04d}: Jade export {jpath.name} covers '
                  f'{jt.min():.1f}-{jt.max():.1f} deg but the scan shown covers '
                  f'{measured[0].min():.1f}-{measured[0].max():.1f} deg')
    else:
        solution = simulate(solution_spec)

    rec = {
        'sample_number': num,
        'formula': formula,
        'outcome': row['XRD Result'],
        'measured': measured,
        'measured_file': path,
        'predicted': predicted,
        'predicted_space_group': space_group_symbol(oqmd_cif),
        'solution': solution,
        'solution_label': solution_label,
        'solution_detail': solution_detail,
    }
    print(f'  {num:04d} {formula:<10} {row["XRD Result"]:<3} '
          f'{path.name:<34} predicted={"yes" if predicted else "NO":<3} '
          f'solution={"yes" if solution else "-"}')
    return rec


def build_records(df):
    """One record per row of ``df`` that has a measurement, in ``df`` order."""
    records = [build_record(row) for _, row in df.iterrows()]
    return [r for r in records if r is not None]


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------


def panel_letter(index):
    """1 -> 'a', 2 -> 'b', ...  Figures here never exceed 26 panels."""
    return string.ascii_lowercase[index - 1]


def _traces(rec):
    """(pattern, colour, linewidth, label) bottom to top."""
    traces = []
    if rec['predicted'] is not None:
        label = LEGEND_PREDICTED
        if rec['predicted_space_group']:
            label += f" ({_mathtext_space_group(rec['predicted_space_group'])})"
        traces.append((rec['predicted'], COLOR_PREDICTED, LW_SIMULATED, label))
    traces.append((rec['measured'], COLOR_MEASURED, LW_MEASURED, 'Measured'))
    if rec['solution'] is not None:
        traces.append((rec['solution'], COLOR_SOLUTION, LW_SIMULATED,
                       rec.get('solution_detail') or rec['solution_label'] or LEGEND_SOLUTION))
    return traces


def _mathtext_space_group(symbol):
    """'P6_3mc' -> '$P6_3mc$', 'Fm-3m' -> '$Fm\\bar{3}m$'."""
    out = symbol.replace('-3', r'\bar{3}').replace('-4', r'\bar{4}') \
                .replace('-6', r'\bar{6}').replace('-1', r'\bar{1}')
    return f'${out}$'


def draw_panel(ax, rec, letter, *, fontsize=6.5, label_traces=False,
               show_sample_number=True, title=None, inset=None, minor=XTICK_MINOR):
    """Draw one sample's stacked traces into ``ax``.

    ``label_traces`` writes each trace's name at its right-hand end, which
    suits large panels; small panels rely on a figure legend instead.
    ``title`` goes above the panel, after the letter.  ``inset`` = (lo, hi)
    adds a zoom of that 2-theta window with measured and predicted overlaid,
    scaled to the predicted intensity there: for showing a missing line.
    """
    traces = _traces(rec)
    offsets = [i * (1.0 + GAP) for i in range(len(traces))]
    for offset, ((tt, ii), colour, lw, _) in zip(offsets, traces):
        ax.plot(tt, ii + offset, color=colour, lw=lw, rasterized=True,
                solid_joinstyle='round')

    ax.set_xlim(*TWO_THETA_RANGE)
    # Headroom above the top trace so the corner labels sit clear of it; an
    # inset gets its own band above that so it covers no data.
    INSET_BAND = 1.0
    y_top = offsets[-1] + 1.30 + (INSET_BAND if inset is not None else 0.0)
    ax.set_ylim(-0.08, y_top)
    ax.set_yticks([])
    lo, hi = TWO_THETA_RANGE
    ax.xaxis.set_major_locator(MultipleLocator(XTICK_MAJOR))
    ax.xaxis.set_minor_locator(MultipleLocator(minor))
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda x, _: f'{x:g}' if lo < x < hi else ''))
    ax.tick_params(axis='x', which='major', direction='in', top=True, length=3.5,
                   labelbottom=False, labelsize=fontsize, pad=2)
    ax.tick_params(axis='x', which='minor', direction='in', top=True, length=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    if inset is not None and rec['predicted'] is not None:
        lo, hi = inset
        ax.axvspan(lo, hi, color='0.5', alpha=0.10, zorder=0)
        band_lo = (offsets[-1] + 1.30 + 0.08) / (y_top + 0.08)
        band_h = INSET_BAND / (y_top + 0.08)
        axins = ax.inset_axes([0.36, band_lo + 0.12 * band_h, 0.28, 0.70 * band_h])
        for (tt, ii), colour, lw in ((rec['predicted'], COLOR_PREDICTED, LW_SIMULATED),
                                     (rec['measured'], COLOR_MEASURED, LW_MEASURED)):
            m = (tt >= lo) & (tt <= hi)
            axins.plot(tt[m], ii[m], color=colour, lw=lw)
        pt, pi = rec['predicted']; mt, mi = rec['measured']
        top = max(pi[(pt >= lo) & (pt <= hi)].max(), mi[(mt >= lo) & (mt <= hi)].max())
        axins.set_xlim(lo, hi); axins.set_ylim(-0.02 * top, 1.35 * top)
        axins.set_yticks([]); axins.tick_params(axis='x', labelsize=fontsize - 2, pad=1, length=2)
        axins.set_title(f'{lo:g}–{hi:g}°, same scale', fontsize=fontsize - 2, pad=2)
        for spine in axins.spines.values():
            spine.set_linewidth(0.5)
        axins.set_facecolor('white')

    if label_traces:
        for offset, (_, colour, _, label) in zip(offsets, traces):
            if '\n' in label:
                # Multi-line labels (phase lists) hang down from the top of the
                # trace band, over its quiet high-angle end.
                ax.text(0.99, offset + 0.84, label, transform=ax.get_yaxis_transform(),
                        fontsize=fontsize - 0.5, color=colour, ha='right', va='top',
                        linespacing=1.15,
                        bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none',
                                  alpha=0.8))
            else:
                ax.text(0.99, offset + 0.55, label, transform=ax.get_yaxis_transform(),
                        fontsize=fontsize - 0.5, color=colour, ha='right', va='bottom',
                        bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none',
                                  alpha=0.8))

    formula = ax.text(0.975, 0.955, formula_to_mathtext(rec['formula']),
                      transform=ax.transAxes, fontsize=fontsize, va='top',
                      ha='right',
                      bbox=dict(boxstyle='square,pad=0.12', fc='white',
                                ec='none', alpha=0.75))
    if not label_traces and rec['solution'] is not None:
        # Name the green trace under the formula; the legend only says
        # "Solution".
        ax.annotate(rec['solution_label'] or LEGEND_SOLUTION, xycoords=formula,
                    xy=(1.0, 0.0), xytext=(0, -1.5), textcoords='offset points',
                    fontsize=fontsize - 1, color=COLOR_SOLUTION, va='top',
                    ha='right',
                    bbox=dict(boxstyle='square,pad=0.12', fc='white', ec='none',
                              alpha=0.75))

    if title is None:
        tag = ax.text(0.025, 0.955, f'({letter})', transform=ax.transAxes,
                      fontsize=fontsize + 0.5, fontweight='bold', va='top',
                      ha='left')
    else:
        ax.set_title(f'({letter}) {title}', loc='left', fontsize=fontsize + 0.5,
                     fontweight='bold', pad=3)
        tag = None

    if show_sample_number:
        text = f"{rec['sample_number']:04d}"
        if tag is not None:
            # Beside the panel letter rather than below it: the traces come
            # close to the top-left corner in the busier panels.
            ax.annotate(text, xycoords=tag, xy=(1.0, 0.0), xytext=(3, 0),
                        textcoords='offset points', fontsize=fontsize - 1.5,
                        color='#888888', va='baseline', ha='left')
        else:
            ax.text(0.025, 0.955, text, transform=ax.transAxes,
                    fontsize=fontsize - 1.5, color='#888888', va='top', ha='left')


def legend_handles(records):
    handles = [
        Line2D([0], [0], color=COLOR_PREDICTED, lw=1.2, label=LEGEND_PREDICTED),
        Line2D([0], [0], color=COLOR_MEASURED, lw=1.2, label=LEGEND_MEASURED),
    ]
    if any(rec['solution'] is not None for rec in records):
        handles.append(Line2D([0], [0], color=COLOR_SOLUTION, lw=1.2,
                              label=LEGEND_SOLUTION))
    return handles


def label_outer_axes(axes, *, fontsize):
    """x tick labels on the lowest visible panel of each column, y labels on
    the leftmost panel of every row."""
    nrows, ncols = axes.shape
    for c in range(ncols):
        for r in range(nrows - 1, -1, -1):
            if axes[r, c].get_visible():
                axes[r, c].tick_params(labelbottom=True)
                axes[r, c].set_xlabel(XLABEL, fontsize=fontsize + 1, labelpad=1)
                break
    for r in range(nrows):
        if axes[r, 0].get_visible():
            axes[r, 0].set_ylabel(YLABEL, fontsize=fontsize + 0.5, labelpad=2)
