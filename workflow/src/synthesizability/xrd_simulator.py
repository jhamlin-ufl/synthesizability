"""Powder XRD pattern simulation from a CIF.

Ported from ``search_match/xrd_simulator.py`` (github.com/jhamlin-ufl/search_match)
so that this repository is a self-contained record of the analysis.  The physics
is unchanged -- pymatgen's ``XRDCalculator`` for structure factors, a Caglioti
FWHM, Gaussian peak profiles -- but the synthetic background, synthetic noise
and the interactive test suite of the original have been dropped: here we
compare against *measured* patterns that have had their own background removed,
so a fabricated background would only have to be subtracted again.

The single knob that matters for the figures is ``fwhm``.  With U = V = 0 the
Caglioti expression reduces to FWHM = sqrt(W), so ``fwhm`` is set directly
rather than through W, which is what the original code did indirectly (and
inconsistently: two call sites used W = 0.2*(fwhm/0.5) and W = 0.2*(fwhm/0.5)^2
for the same intended width).
"""
from __future__ import annotations

import numpy as np
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser

# Cu K-alpha wavelengths (Angstrom)
WAVELENGTHS = {
    'CuKa': 1.5418,    # Ka1/Ka2 intensity-weighted average
    'CuKa1': 1.5406,
}


def _resolve_wavelength(wavelength):
    if isinstance(wavelength, str):
        try:
            return WAVELENGTHS[wavelength]
        except KeyError:
            raise ValueError(f'unknown wavelength {wavelength!r}; '
                             f'use a float or one of {sorted(WAVELENGTHS)}')
    return float(wavelength)


def _load_structure(cif):
    if isinstance(cif, Structure):
        return cif
    return CifParser(str(cif)).parse_structures(primitive=True)[0]


def caglioti_fwhm(two_theta, U=0.0, V=0.0, W=0.25):
    """Gaussian FWHM in degrees 2-theta: sqrt(U tan^2(th) + V tan(th) + W)."""
    tan_theta = np.tan(np.radians(np.asarray(two_theta, dtype=float) / 2))
    fwhm_sq = U * tan_theta ** 2 + V * tan_theta + W
    return np.sqrt(np.maximum(fwhm_sq, 1e-4))


def peak_list(cif, two_theta_range=(10, 135), wavelength='CuKa'):
    """Bragg peak positions and relative intensities (max = 100)."""
    structure = _load_structure(cif)
    calculator = XRDCalculator(wavelength=_resolve_wavelength(wavelength))
    pattern = calculator.get_pattern(structure, two_theta_range=two_theta_range)

    angles = np.asarray(pattern.x, dtype=float)
    intensities = np.asarray(pattern.y, dtype=float)
    if intensities.size and intensities.max() > 0:
        intensities = intensities / intensities.max() * 100.0
    return angles, intensities


def simulate_pattern(cif, two_theta_range=(20, 120), step_size=0.02,
                     wavelength='CuKa', fwhm=0.5, U=0.0, V=0.0):
    """Simulate a background-free powder pattern.

    Peaks are computed over a range padded by 5 x fwhm on each side so that a
    reflection just outside the plotted window still contributes its tail.

    Returns
    -------
    (two_theta, intensity) : tuple of ndarray
        ``intensity`` is normalised to a maximum of 1 (all zeros if the
        structure has no reflections in range).
    """
    lo, hi = two_theta_range
    pad = 5 * fwhm
    angles, rel_intensities = peak_list(
        cif, two_theta_range=(max(lo - pad, 0.5), min(hi + pad, 179.5)),
        wavelength=wavelength)

    two_theta = np.arange(lo, hi + step_size, step_size)
    intensity = np.zeros_like(two_theta)

    W = fwhm ** 2
    for centre, height in zip(angles, rel_intensities):
        width = float(caglioti_fwhm(centre, U=U, V=V, W=W))
        sigma = width / (2 * np.sqrt(2 * np.log(2)))
        intensity += height * np.exp(-0.5 * ((two_theta - centre) / sigma) ** 2)

    if intensity.max() > 0:
        intensity = intensity / intensity.max()
    return two_theta, intensity
