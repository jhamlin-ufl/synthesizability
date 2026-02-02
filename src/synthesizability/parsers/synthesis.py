# src/synthesizability/parsers/synthesis.py
"""
Parser for SYNTHESIS files containing synthesis procedure and mass loss information.
"""

import re


def parse_synthesis_file(content: str) -> dict:
    """
    Parse structured data from SYNTHESIS file.
    
    Args:
        content: String content of SYNTHESIS file
        
    Returns:
        dict with keys:
            - 'mass_loss_percent': float or None - mass loss percentage
            - 'initial_mass_g': float or None - initial mass in grams
            - 'final_mass_g': float or None - final mass in grams
    """
    if not content:
        return {
            'mass_loss_percent': None,
            'initial_mass_g': None,
            'final_mass_g': None,
        }
    
    mass_loss, initial_mass, final_mass = _extract_mass_data(content)
    
    return {
        'mass_loss_percent': mass_loss,
        'initial_mass_g': initial_mass,
        'final_mass_g': final_mass,
    }


def _extract_mass_data(synthesis_content: str) -> tuple:
    """Extract mass loss percentage, initial and final masses."""
    # Extract mass loss percentage (case insensitive)
    mass_loss = None
    match = re.search(r'loss:\s*([\d.]+)%', synthesis_content, re.IGNORECASE)
    if match:
        mass_loss = float(match.group(1))
    
    # Extract initial mass (case insensitive)
    initial_mass = None
    match = re.search(r'initial mass:\s*([\d.]+)\s*g', synthesis_content, re.IGNORECASE)
    if match:
        initial_mass = float(match.group(1))
    
    # Extract final mass (case insensitive)
    final_mass = None
    match = re.search(r'final mass:\s*([\d.]+)\s*g', synthesis_content, re.IGNORECASE)
    if match:
        final_mass = float(match.group(1))
    
    return mass_loss, initial_mass, final_mass