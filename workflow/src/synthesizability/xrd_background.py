"""Background estimation and removal for measured powder XRD patterns.

The backgrounds in this study are dominated by two things that have nothing to
do with the sample's crystal structure: the amorphous hump from the PMMA slide
and the petroleum jelly / vacuum grease used to hold the powder (a broad maximum
near 2-theta ~ 12 deg, with a long tail), and a smooth instrumental / air-scatter
trend that is falling for most patterns but rising for the bulk pieces measured
without a slide.  Comparing a measured pattern to a simulated one is much easier
once that is gone.

``search_match`` did this by building one global background shape from the 10th
percentile across all files in a directory, fitting a Chebyshev polynomial to it
and then scaling that single shape per file.  That works when every pattern comes
off the same instrument in the same geometry.  It does not work here: these 33
patterns come from two diffractometers over two angular ranges, and the bulk
pieces have backgrounds that *rise* with angle, which no rescaling of a common
falling shape can reproduce.

We therefore estimate a background per pattern using SNIP (statistics-sensitive
non-linear iterative peak clipping, Ryan et al., Nucl. Instrum. Methods B 34, 396
(1988)), the standard choice for this in diffraction and gamma spectroscopy.
SNIP has one meaningful parameter -- the clipping window -- and it has a direct
physical reading: features broader than roughly twice the window are treated as
background, features narrower than it survive.  Nothing is fitted, so it cannot
diverge at the ends of the range the way a high-order polynomial does.
"""
from __future__ import annotations

import numpy as np

# Default clipping half-width in degrees 2-theta.  Bragg peaks here are at most
# ~1.5 deg FWHM even for the broadest solid solutions; the amorphous hump is
# ~10 deg wide.  4 deg sits comfortably between the two.
DEFAULT_WINDOW_DEG = 4.0

# Default smoothing applied before clipping, in degrees 2-theta.  This only
# damps point-to-point counting noise so the clipping does not chase it; it is
# far narrower than any real peak.
DEFAULT_SMOOTH_DEG = 0.06


def _lls(y):
    """Log-log-square-root transform: compresses dynamic range for SNIP."""
    return np.log(np.log(np.sqrt(y + 1.0) + 1.0) + 1.0)


def _inv_lls(z):
    return (np.exp(np.exp(z) - 1.0) - 1.0) ** 2 - 1.0


def _smooth(y, n_points):
    """Moving average over an odd number of points, with edge padding."""
    if n_points < 3:
        return y
    if n_points % 2 == 0:
        n_points += 1
    half = n_points // 2
    padded = np.pad(y, half, mode='edge')
    kernel = np.ones(n_points) / n_points
    return np.convolve(padded, kernel, mode='valid')


def snip_background(two_theta, intensity, window_deg=DEFAULT_WINDOW_DEG,
                    smooth_deg=DEFAULT_SMOOTH_DEG, decreasing=True):
    """Estimate the background of one pattern by SNIP peak clipping.

    Parameters
    ----------
    two_theta, intensity : array_like
        The measured pattern, on a monotonically increasing and (to good
        approximation) uniform 2-theta grid.
    window_deg : float
        Clipping half-width in degrees.  Features much broader than this are
        background; features much narrower are peaks.
    smooth_deg : float
        Width of the moving average applied before clipping.  Set to 0 to skip.
    decreasing : bool
        Run the iterations with a shrinking window (the usual "adaptive" SNIP).
        This removes broad structure first and is less prone to eating the
        flanks of strong peaks than the increasing-window form.

    Returns
    -------
    ndarray
        The background, on the same grid and in the same units as ``intensity``.
    """
    two_theta = np.asarray(two_theta, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    if two_theta.size != intensity.size:
        raise ValueError('two_theta and intensity must have the same length')
    if two_theta.size < 5:
        return np.full_like(intensity, intensity.min())

    step = float(np.median(np.diff(two_theta)))
    if step <= 0:
        raise ValueError('two_theta must be increasing')

    # SNIP works on a non-negative signal; shift so the minimum sits at zero.
    floor = intensity.min()
    y = intensity - floor

    if smooth_deg > 0:
        y = _smooth(y, int(round(smooth_deg / step)))

    max_window = max(int(round(window_deg / step)), 1)

    # SNIP holds the outermost `w` points fixed, which pins the background to
    # the data at both ends of the range.  On a pattern whose range starts part
    # way down the flank of the amorphous hump that is badly wrong: the estimate
    # can only descend from the pinned value, so the flank survives subtraction.
    # Reflecting the pattern about each end gives the clipping something to work
    # against, and the reflected part is discarded afterwards.
    pad = max_window
    z = _lls(np.concatenate([y[pad:0:-1], y, y[-2:-pad - 2:-1]]))

    windows = range(max_window, 0, -1) if decreasing else range(1, max_window + 1)
    for w in windows:
        # z[i] <- min(z[i], mean(z[i-w], z[i+w]))
        shifted_mean = 0.5 * (z[:-2 * w] + z[2 * w:])
        z[w:-w] = np.minimum(z[w:-w], shifted_mean)

    return _inv_lls(z[pad:pad + y.size]) + floor


def subtract_background(two_theta, intensity, clip_negative=True, **kwargs):
    """Return ``(background_subtracted, background)`` for one pattern.

    ``clip_negative`` floors the result at zero.  SNIP is a lower envelope, so
    the residual is non-negative except where counting noise dips below it; the
    clipped-away part is at the level of the noise.
    """
    background = snip_background(two_theta, intensity, **kwargs)
    subtracted = np.asarray(intensity, dtype=float) - background
    if clip_negative:
        subtracted = np.maximum(subtracted, 0.0)
    return subtracted, background
