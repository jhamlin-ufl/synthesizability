# src/synthesizability/parsers/xrd.py
"""
XRD file parser for Siemens D500 (.txt) and Panalytical (.xy) formats.
"""

import re
from pathlib import Path
import numpy as np


def is_xrd_file(filepath: Path) -> bool:
    """
    Check if a file contains XRD data by inspecting its contents.
    
    Returns:
        True if file appears to be XRD data (2theta vs intensity)
    """
    # Exclude chi (magnetic susceptibility) files
    fname_lower = filepath.name.lower()
    if 'chi' in fname_lower:
        return False

    # Exclude Jade WPF fitting output files (Angle/I(o)/I(c)/I(d)/I(b) format)
    if 'xrd_fit' in fname_lower:
        return False
    
    # Only process .txt and .xy files
    if filepath.suffix not in ['.txt', '.xy']:
        return False
    
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
            
        if len(lines) < 10:
            return False
        
        # Check if it's a Siemens RAW file
        if lines[0].strip().startswith(';RAW'):
            return True
        
        # For other files, check numerical data range
        # XRD has 2theta typically in range 5-150 degrees
        sample_lines = [l for l in lines[:50] if l.strip()]
        numerical_lines = 0
        valid_xrd_range = 0
        
        for line in sample_lines:
            try:
                parts = line.strip().split()
                if len(parts) >= 2:
                    x_val = float(parts[0])
                    numerical_lines += 1
                    
                    # XRD 2theta is typically 5-150 degrees
                    if 5 <= x_val <= 150:
                        valid_xrd_range += 1
                        
            except (ValueError, IndexError):
                continue
        
        # Require most data points in valid XRD 2theta range
        if numerical_lines >= 10 and valid_xrd_range / numerical_lines > 0.7:
            return True
            
        return False
        
    except Exception:
        return False


def parse_xrd_file(filepath: Path) -> dict:
    """
    Parse XRD file (auto-detects format).
    
    Args:
        filepath: Path to XRD file
        
    Returns:
        dict with XRD pattern data
        
    Raises:
        ValueError: If file is not recognized as XRD data
    """
    filepath = Path(filepath)
    
    if not is_xrd_file(filepath):
        raise ValueError(f"File does not appear to contain XRD data: {filepath}")
    
    # Try to determine format
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        first_line = f.readline().strip()
    
    if first_line.startswith(';RAW'):
        return _parse_siemens_txt(filepath)
    else:
        # Assume Panalytical-style two-column format
        return _parse_panalytical_xy(filepath)


def _parse_siemens_txt(filepath: Path) -> dict:
    """Parse Siemens D500 .txt file in RAW4.00 format."""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    
    # Extract metadata from header
    metadata = _extract_siemens_metadata(lines)
    
    # Find where data starts (look for lines with comma-separated numbers)
    data_start = None
    for i, line in enumerate(lines):
        # Skip empty lines and headers
        if not line.strip() or line.startswith('[') or line.startswith(';') or '=' in line:
            continue
        # Check if this looks like data (contains comma and numbers)
        if ',' in line and re.search(r'\d+\.\d+', line):
            data_start = i
            break
    
    if data_start is None:
        raise ValueError(f"Could not find data section in {filepath}")
    
    # Parse data
    two_theta = []
    intensity = []
    
    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        
        try:
            parts = line.split(',')
            if len(parts) >= 2:
                two_theta.append(float(parts[0].strip()))
                intensity.append(float(parts[1].strip()))
        except (ValueError, IndexError):
            continue
    
    two_theta = np.array(two_theta)
    intensity = np.array(intensity)
    
    # Calculate statistics
    step_sizes = np.diff(two_theta)
    avg_step_size = np.mean(step_sizes) if len(step_sizes) > 0 else 0.0
    
    return {
        'two_theta': two_theta,
        'intensity': intensity,
        'instrument': 'Siemens',
        'two_theta_min': float(np.min(two_theta)) if len(two_theta) > 0 else None,
        'two_theta_max': float(np.max(two_theta)) if len(two_theta) > 0 else None,
        'n_points': len(two_theta),
        'step_size': float(avg_step_size),
        'date': metadata.get('date'),
        'anode': metadata.get('anode')
    }


def _parse_panalytical_xy(filepath: Path) -> dict:
    """Parse Panalytical .xy file (simple two-column format)."""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    
    two_theta = []
    intensity = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        try:
            parts = line.split()
            if len(parts) >= 2:
                two_theta.append(float(parts[0]))
                intensity.append(float(parts[1]))
        except (ValueError, IndexError):
            continue
    
    two_theta = np.array(two_theta)
    intensity = np.array(intensity)
    
    # Calculate statistics
    step_sizes = np.diff(two_theta)
    avg_step_size = np.mean(step_sizes) if len(step_sizes) > 0 else 0.0
    
    return {
        'two_theta': two_theta,
        'intensity': intensity,
        'instrument': 'Panalytical',
        'two_theta_min': float(np.min(two_theta)) if len(two_theta) > 0 else None,
        'two_theta_max': float(np.max(two_theta)) if len(two_theta) > 0 else None,
        'n_points': len(two_theta),
        'step_size': float(avg_step_size),
        'date': None,
        'anode': None
    }


def _extract_siemens_metadata(lines: list) -> dict:
    """Extract metadata from Siemens file header."""
    metadata = {}
    
    for line in lines:
        line = line.strip()
        
        # Look for Date=
        if line.startswith('Date='):
            metadata['date'] = line.split('=', 1)[1].strip()
        
        # Look for Anode=
        if line.startswith('Anode='):
            metadata['anode'] = line.split('=', 1)[1].strip()
        
        # Stop when we hit the data section
        if ',' in line and re.search(r'\d+\.\d+', line):
            break
    
    return metadata


def get_xrd_summary(xrd_dict: dict) -> dict:
    """
    Extract dataframe-friendly summary from XRD pattern dict.
    
    Args:
        xrd_dict: Output from parse_xrd_file()
        
    Returns:
        dict with summary info suitable for dataframe columns
    """
    return {
        'xrd_two_theta_range': (xrd_dict['two_theta_min'], xrd_dict['two_theta_max']),
        'xrd_n_points': xrd_dict['n_points'],
        'xrd_instrument': xrd_dict['instrument']
    }