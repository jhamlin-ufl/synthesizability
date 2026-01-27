# src/synthesizability/parsers/xrd.py
"""
XRD file parser for Siemens D500 (.txt) and Panalytical (.xy) formats.
"""

import re
from pathlib import Path
import numpy as np


def parse_xrd_file(filepath: Path) -> dict:
    """
    Parse XRD file (Siemens .txt or Panalytical .xy).
    
    Args:
        filepath: Path to XRD file
        
    Returns:
        dict with keys:
            - 'two_theta': np.ndarray - 2theta angles
            - 'intensity': np.ndarray - intensity values
            - 'instrument': str - 'Siemens D500' or 'Panalytical'
            - 'two_theta_min': float - minimum 2theta
            - 'two_theta_max': float - maximum 2theta
            - 'n_points': int - number of data points
            - 'step_size': float - average step size in 2theta
            - 'date': str or None - measurement date (Siemens only)
            - 'anode': str or None - X-ray anode material (Siemens only)
    """
    filepath = Path(filepath)
    
    if filepath.suffix == '.txt':
        return _parse_siemens_txt(filepath)
    elif filepath.suffix == '.xy':
        return _parse_panalytical_xy(filepath)
    else:
        raise ValueError(f"Unknown XRD file type: {filepath.suffix}")


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
        'instrument': 'Siemens D500',
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