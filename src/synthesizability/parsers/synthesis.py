# src/synthesizability/parsers/synthesis.py
"""
Parser for SYNTHESIS files containing synthesis procedure and mass loss information.
"""

import re
from pymatgen.core import Composition, Element


def parse_synthesis_file(content: str, formula: str = None) -> dict:
    """
    Parse structured data from SYNTHESIS file.

    Args:
        content: String content of SYNTHESIS file
        formula: Formula string from directory name (e.g. 'MoTiTa2'), used to
                 compute expected mole fractions for composition deviation check

    Returns:
        dict with keys:
            - 'mass_loss_percent': float or None
            - 'initial_mass_g': float or None
            - 'final_mass_g': float or None
            - 'composition_max_deviation': float or None
            - 'composition_euclidean_deviation': float or None
            - 'composition_measured_fractions': dict or None (pickle only)
            - 'composition_expected_fractions': dict or None (pickle only)
    """
    base = {
        'mass_loss_percent': None,
        'initial_mass_g': None,
        'final_mass_g': None,
        'composition_max_deviation': None,
        'composition_euclidean_deviation': None,
        'composition_measured_fractions': None,
        'composition_expected_fractions': None,
        'composition_ok': None
    }

    if not content:
        return base

    mass_loss, initial_mass, final_mass = _extract_mass_data(content)
    base['mass_loss_percent'] = mass_loss
    base['initial_mass_g'] = initial_mass
    base['final_mass_g'] = final_mass

    if formula:
        deviation_data = _compute_composition_deviation(content, formula)
        base.update(deviation_data)

    return base


def _extract_mass_data(synthesis_content: str) -> tuple:
    """Extract mass loss percentage, initial and final masses."""
    mass_loss = None
    match = re.search(r'loss:\s*([\d.]+)%', synthesis_content, re.IGNORECASE)
    if match:
        mass_loss = float(match.group(1))

    initial_mass = None
    match = re.search(r'initial mass:\s*([\d.]+)\s*g', synthesis_content, re.IGNORECASE)
    if match:
        initial_mass = float(match.group(1))

    final_mass = None
    match = re.search(r'final mass:\s*([\d.]+)\s*g', synthesis_content, re.IGNORECASE)
    if match:
        final_mass = float(match.group(1))

    return mass_loss, initial_mass, final_mass


def _parse_measured_masses(content: str) -> dict | None:
    """
    Extract measured masses (in grams) per element from SYNTHESIS file content.

    Handles known format variations:
        - 'El: 0.1234 g'  (standard)
        - 'El: 0.1234g'   (no space before g)
        - 'El 0.1234 g:'  (colon after value)
        - 'GeL 0.1234g'   (typo: L instead of colon)
        - 'Total: ...'    (ignored)

    Returns:
        dict mapping element symbol -> mass in grams, or None if not parseable
    """
    # Find the measured masses block
    match = re.search(
        r'measured\s+masses?\s*:?\s*\n(.*?)(?:\n\s*\n|\ninitial|\nfinal|\nloss|\nheating|\nprocedure|\ntarget)',
        content,
        re.IGNORECASE | re.DOTALL
    )
    if not match:
        return None

    block = match.group(1)

    # Split on commas to get individual element entries
    entries = [e.strip() for e in block.split(',')]

    masses = {}
    for entry in entries:
        if not entry:
            continue

        # Skip total field
        if re.match(r'total', entry, re.IGNORECASE):
            continue

        # Standard: 'El: 0.1234 g' or 'El: 0.1234g'
        m = re.match(r'^([A-Z][a-z]?)\s*:\s*([\d.]+)\s*g', entry)
        if m:
            masses[m.group(1)] = float(m.group(2))
            continue

        # Swapped: 'El 0.1234 g:' (colon after value)
        m = re.match(r'^([A-Z][a-z]?)\s+([\d.]+)\s*g:', entry)
        if m:
            masses[m.group(1)] = float(m.group(2))
            continue

        # Typo: 'ElL 0.1234g' (L instead of colon, e.g. GeL)
        m = re.match(r'^([A-Z][a-z]?)L\s+([\d.]+)\s*g', entry)
        if m:
            masses[m.group(1)] = float(m.group(2))
            continue

    return masses if masses else None


def _compute_composition_deviation(content: str, formula: str) -> dict:
    """
    Compute deviation between measured and expected mole fractions.

    Returns:
        dict with composition_max_deviation, composition_euclidean_deviation,
        composition_measured_fractions, composition_expected_fractions
    """
    null_result = {
        'composition_max_deviation': None,
        'composition_euclidean_deviation': None,
        'composition_measured_fractions': None,
        'composition_expected_fractions': None,
        'composition_ok': None
    }

    measured_masses = _parse_measured_masses(content)
    if not measured_masses:
        return null_result

    # Compute expected mole fractions from formula
    try:
        comp = Composition(formula)
        expected_fractions = {
            str(el): amt / comp.num_atoms
            for el, amt in comp.items()
        }
    except Exception:
        return null_result

    # Check that measured elements match expected
    if set(measured_masses.keys()) != set(expected_fractions.keys()):
        return null_result

    # Convert measured masses to moles using pymatgen atomic masses
    try:
        moles = {
            el: mass / float(Element(el).atomic_mass)
            for el, mass in measured_masses.items()
        }
    except Exception:
        return null_result

    total_moles = sum(moles.values())
    if total_moles == 0:
        return null_result

    measured_fractions = {el: n / total_moles for el, n in moles.items()}

    # Compute deviations over the shared element set
    elements = list(expected_fractions.keys())
    max_dev = max(
        abs(measured_fractions[el] - expected_fractions[el])
        for el in elements
    )
    euclidean_dev = sum(
        (measured_fractions[el] - expected_fractions[el]) ** 2
        for el in elements
    ) ** 0.5

    return {
        'composition_max_deviation': round(max_dev, 6),
        'composition_euclidean_deviation': round(euclidean_dev, 6),
        'composition_measured_fractions': measured_fractions,
        'composition_expected_fractions': expected_fractions,
        'composition_ok': max_dev <= 0.02,
    }
