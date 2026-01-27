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
            - 'has_powder_premelting': bool - whether powders were pre-melted
            - 'air_sensitive_handling': bool - whether air-sensitive procedures used
    """
    if not content:
        return {
            'mass_loss_percent': None,
            'initial_mass_g': None,
            'final_mass_g': None,
            'has_powder_premelting': False,
            'air_sensitive_handling': False
        }
    
    mass_loss, initial_mass, final_mass = _extract_mass_data(content)
    
    # Check for powder premelting (case insensitive)
    has_premelting = bool(re.search(r'(arc melted|arced|pelleted).+(powder|pellet).+before', 
                                     content, re.IGNORECASE))
    
    # Check for air-sensitive handling (case insensitive)
    has_air_sensitive = bool(re.search(r'(air sensitive|glovebox|outside of glovebox)', 
                                        content, re.IGNORECASE))
    
    return {
        'mass_loss_percent': mass_loss,
        'initial_mass_g': initial_mass,
        'final_mass_g': final_mass,
        'has_powder_premelting': has_premelting,
        'air_sensitive_handling': has_air_sensitive
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