# src/synthesizability/parsers/status.py
"""
Parser for STATUS files containing synthesis outcome information.
"""

import re


def parse_status_file(content: str) -> dict:
    """
    Parse structured data from STATUS file.
    
    Args:
        content: String content of STATUS file
        
    Returns:
        dict with keys:
            - 'superconductivity': str - raw superconductivity text
            - 'tc_kelvin': float or None - extracted Tc value
            - 'xrd_type': str or None - 'Bulk' or 'Powder'
            - 'xrd_instrument': str or None - XRD instrument used
            - 'xrd_result': str or None - phase description
            - 'prediction_list': str or None - prediction source/method
    """
    if not content:
        return {
            'superconductivity': None,
            'tc_kelvin': None,
            'xrd_type': None,
            'xrd_instrument': None,
            'xrd_result': None,
            'prediction_list': None
        }
    
    data = {}
    
    # Extract superconductivity line (case insensitive)
    sc_match = re.search(r'superconductivity:\s*(.+?)(?:\n|$)', content, re.IGNORECASE)
    data['superconductivity'] = sc_match.group(1).strip() if sc_match else None
    data['tc_kelvin'] = _extract_tc_value(data['superconductivity']) if data['superconductivity'] else None
    
    # Extract XRD info (case insensitive)
    xrd_match = re.search(r'xrd:\s*(.+?)(?:\n|$)', content, re.IGNORECASE)
    if xrd_match:
        xrd_text = xrd_match.group(1).strip()
        
        # XRD type (case insensitive)
        if re.search(r'bulk', xrd_text, re.IGNORECASE):
            data['xrd_type'] = 'Bulk'
        elif re.search(r'powder', xrd_text, re.IGNORECASE):
            data['xrd_type'] = 'Powder'
        else:
            data['xrd_type'] = None
        
        # XRD instrument (case insensitive)
        if re.search(r'nrf\s+xrd', xrd_text, re.IGNORECASE):
            data['xrd_instrument'] = 'NRF XRD'
        elif re.search(r'hamlin\s+xrd', xrd_text, re.IGNORECASE):
            data['xrd_instrument'] = 'Hamlin XRD'
        else:
            data['xrd_instrument'] = None
        
        # XRD result (everything after the instrument)
        parts = xrd_text.split(',')
        if len(parts) >= 3:
            data['xrd_result'] = parts[2].strip()
        elif len(parts) == 2:
            data['xrd_result'] = parts[1].strip()
        else:
            data['xrd_result'] = None
    else:
        data['xrd_type'] = None
        data['xrd_instrument'] = None
        data['xrd_result'] = None
    
    # Extract List category (case insensitive)
    list_match = re.search(r'list:\s*(.+?)(?:\n|$)', content, re.IGNORECASE)
    data['prediction_list'] = list_match.group(1).strip() if list_match else None
    
    return data


def _extract_tc_value(superconductivity_text: str) -> float:
    """Extract Tc value in Kelvin from superconductivity text."""
    if not superconductivity_text or 'not' in superconductivity_text.lower():
        return None
    
    # Look for pattern like "Tc of X.X K" or "Tc onset ~ X.X K" (case insensitive)
    match = re.search(r'tc\s+(?:of|onset)\s*~?\s*([\d.]+)\s*k', superconductivity_text, re.IGNORECASE)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            # Handle malformed numbers like "5.1."
            cleaned = match.group(1).rstrip('.')
            try:
                return float(cleaned)
            except ValueError:
                return None
    
    return None