# src/synthesizability/io/dataframe.py
"""
Build and analyze dataframe from synthesis experiment data.
"""

import re
from pathlib import Path
import pandas as pd
import numpy as np

from ..parsers import parse_status_file, parse_synthesis_file, parse_xrd_file
from ..parsers.xrd import is_xrd_file
from ..formula import enrich_with_formula_properties


def extract_sample_info(dir_path: Path) -> tuple:
    """Extract sample number and formula from directory name."""
    dir_name = dir_path.name
    
    # Extract sample number (leading digits)
    match = re.match(r'^(\d+)_', dir_name)
    sample_number = int(match.group(1)) if match else None
    
    # Extract formula (everything after the second underscore)
    parts = dir_name.split('_', 2)
    formula = parts[2] if len(parts) > 2 else None
    
    return sample_number, formula


def is_chi_file(filepath: Path) -> bool:
    """Check if file is a chi (magnetic susceptibility) file."""
    return 'chi' in filepath.name.lower() and filepath.suffix == '.txt'


def parse_chi_field_from_filename(filename: str) -> float:
    """
    Extract magnetic field value from chi filename.
    
    Example: '20251106_HM499GaHf2Nb_chiAC_vs_T_B_1T.001.txt' -> 1.0
    
    Returns:
        Field value in Tesla, or None if not found
    """
    match = re.search(r'B_([0-9.]+)T', filename)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None


def parse_chi_files(dir_path: Path) -> dict:
    """
    Parse chi (magnetic susceptibility) files in a directory.
    
    Returns:
        dict with chi file info
    """
    chi_files = []
    chi_fields = []
    
    for file in dir_path.iterdir():
        if not file.is_file():
            continue
        
        if is_chi_file(file):
            chi_files.append(file.name)
            
            # Extract field value
            field = parse_chi_field_from_filename(file.name)
            if field is not None and field not in chi_fields:
                chi_fields.append(field)
    
    # Sort fields
    chi_fields.sort()
    
    return {
        'chi_files': chi_files,
        'chi_n_files': len(chi_files),
        'chi_has_high_field': any(f > 0 for f in chi_fields),
        'chi_fields': chi_fields
    }


def parse_xrd_files(dir_path: Path) -> list:
    """
    Parse all XRD files in a directory by content inspection.
    
    Returns:
        List of XRD pattern dicts from parse_xrd_file()
    """
    xrd_patterns = []
    
    # Check all files in directory
    for file in dir_path.iterdir():
        if not file.is_file():
            continue
        
        # Skip obviously non-XRD files
        if file.suffix in ['.pptx', '.pdf', '.png', '.jpg']:
            continue
        
        # Try to parse as XRD
        try:
            if is_xrd_file(file):
                pattern = parse_xrd_file(file)
                pattern['filename'] = file.name  # Add filename to pattern
                xrd_patterns.append(pattern)
        except Exception as e:
            # Silently skip files that aren't XRD data
            pass
    
    return xrd_patterns


def get_xrd_summary_columns(xrd_patterns: list) -> dict:
    """
    Extract summary info from XRD patterns for dataframe columns.
    
    Returns:
        dict with summary columns suitable for CSV export
    """
    if not xrd_patterns:
        return {
            'xrd_files': [],
            'xrd_n_files': 0,
            'xrd_instrument': None,
            'xrd_two_theta_min': None,
            'xrd_two_theta_max': None,
        }
    
    # Collect filenames
    filenames = [p['filename'] for p in xrd_patterns]
    
    # Collect ranges from all patterns
    two_theta_mins = [p['two_theta_min'] for p in xrd_patterns if p['two_theta_min'] is not None]
    two_theta_maxs = [p['two_theta_max'] for p in xrd_patterns if p['two_theta_max'] is not None]
    
    # Get overall range
    overall_min = min(two_theta_mins) if two_theta_mins else None
    overall_max = max(two_theta_maxs) if two_theta_maxs else None
    
    # Determine instruments
    instruments = set(p['instrument'] for p in xrd_patterns)
    if len(instruments) > 1:
        instrument_str = "Both"
    elif len(instruments) == 1:
        instrument_str = list(instruments)[0]
    else:
        instrument_str = None
    
    return {
        'xrd_files': filenames,
        'xrd_n_files': len(xrd_patterns),
        'xrd_instrument': instrument_str,
        'xrd_two_theta_min': overall_min,
        'xrd_two_theta_max': overall_max,
    }


def build_dataframe(data_raw_dir: Path) -> pd.DataFrame:
    """Build dataframe from all experiment directories."""
    
    records = []
    
    # Iterate through all directories
    for dir_path in sorted(data_raw_dir.iterdir()):
        if not dir_path.is_dir():
            continue
        
        # Extract basic info
        sample_number, formula = extract_sample_info(dir_path)
        
        # Get list of files
        files = [f.name for f in dir_path.iterdir()]
        
        # Read STATUS file
        status_path = dir_path / "STATUS"
        status_content = status_path.read_text(encoding='utf-8').strip() if status_path.exists() else None
        
        # Read SYNTHESIS file
        synthesis_path = dir_path / "SYNTHESIS"
        synthesis_content = synthesis_path.read_text(encoding='utf-8').strip() if synthesis_path.exists() else None
        
        # Parse structured data using parsers
        status_data = parse_status_file(status_content)
        synthesis_data = parse_synthesis_file(synthesis_content, formula)
        
        # Parse XRD files
        xrd_patterns = parse_xrd_files(dir_path)
        xrd_summary = get_xrd_summary_columns(xrd_patterns)
        
        # Parse chi files
        chi_data = parse_chi_files(dir_path)
        
        # Check for summary file
        has_summary = any(f.endswith('.pptx') for f in files)
        
        # Build record
        record = {
            'sample_number': sample_number,
            'sample_id': dir_path.name,
            'formula': formula,
            'files': files,
            'status_content': status_content,
            'synthesis_content': synthesis_content,
            'xrd_patterns': xrd_patterns,  # Full patterns (pickle only)
            'has_summary': has_summary,
            **status_data,
            **synthesis_data,
            **xrd_summary,
            **chi_data
        }
        
        records.append(record)
    
    df = pd.DataFrame(records)
    
    # Enrich with formula-derived properties
    df = enrich_with_formula_properties(df)
    
    return df


def analyze_field_statistics(df: pd.DataFrame):
    """Print statistics for each field in the dataframe."""
    
    print("\n" + "="*80)
    print("DATAFRAME FIELD STATISTICS")
    print("="*80)
    print(f"\nTotal samples: {len(df)}\n")
    
    # Analyze each column
    for col in df.columns:
        if col in ['files', 'status_content', 'synthesis_content', 'xrd_patterns',
           'xrd_files', 'chi_files', 'chi_fields',
           'composition_measured_fractions', 'composition_expected_fractions']:
            # Skip the large/complex fields
            continue
        
        print(f"\n{col}")
        print("-" * 40)
        
        # Count non-null values
        non_null = df[col].notna().sum()
        null_count = df[col].isna().sum()
        pct_populated = (non_null / len(df)) * 100
        
        print(f"Populated: {non_null}/{len(df)} ({pct_populated:.1f}%)")
        print(f"Missing: {null_count}")
        
        # Type-specific statistics
        if df[col].dtype == bool:
            true_count = df[col].sum()
            true_pct = (true_count / len(df)) * 100
            print(f"True: {true_count} ({true_pct:.1f}%)")
        
        elif df[col].dtype in [np.float64, np.int64]:
            if non_null > 0:
                print(f"Min: {df[col].min():.2f}")
                print(f"Max: {df[col].max():.2f}")
                print(f"Mean: {df[col].mean():.2f}")
                print(f"Median: {df[col].median():.2f}")
        
        elif df[col].dtype == object and col not in ['formula']:
            # For string fields, show unique values
            unique_vals = df[col].dropna().unique()
            n_unique = len(unique_vals)
            print(f"Unique values: {n_unique}")
            
            if n_unique <= 10:
                # Show all values if there aren't many
                value_counts = df[col].value_counts()
                for val, count in value_counts.items():
                    pct = (count / len(df)) * 100
                    print(f"  '{val}': {count} ({pct:.1f}%)")
            else:
                # Show top 5
                value_counts = df[col].value_counts().head(5)
                print("Top 5 values:")
                for val, count in value_counts.items():
                    pct = (count / len(df)) * 100
                    print(f"  '{val}': {count} ({pct:.1f}%)")
    
    # Summary of potentially problematic fields
    print("\n" + "="*80)
    print("FIELDS WITH LOW POPULATION (< 50%)")
    print("="*80)
    
    for col in df.columns:
        if col in ['files', 'status_content', 'synthesis_content', 'xrd_patterns',
           'xrd_files', 'chi_files', 'chi_fields',
           'composition_measured_fractions', 'composition_expected_fractions']:
            continue
        
        non_null = df[col].notna().sum()
        pct_populated = (non_null / len(df)) * 100
        
        if pct_populated < 50:
            print(f"{col}: {pct_populated:.1f}% populated")


def show_missing_samples(df: pd.DataFrame):
    """Show which samples are missing data for well-populated fields."""
    
    # Fields with >90% population
    fields_to_check = [
        'superconductivity',
        'xrd_type',
        'xrd_instrument',
        'xrd_result',
        'prediction_list',
        'mass_loss_percent',
        'initial_mass_g',
        'final_mass_g'
    ]
    
    print("\n" + "="*80)
    print("SAMPLES MISSING DATA FOR WELL-POPULATED FIELDS")
    print("="*80)
    print("")
    
    for field in fields_to_check:
        if field not in df.columns:
            continue
            
        missing_mask = df[field].isna()
        missing_count = missing_mask.sum()
        
        if missing_count > 0:
            print(f"{field} (missing {missing_count}):")
            print("-" * 40)
            
            missing_samples = df[missing_mask][['sample_number', 'formula']]
            for _, row in missing_samples.iterrows():
                print(f"  {row['sample_number']:04d} - {row['formula']}")
            
            print("")
    
    print("=" * 80)


def get_xrd_pattern(df: pd.DataFrame, sample_number: int, pattern_index: int = 0) -> dict:
    """
    Extract XRD pattern for a specific sample.
    
    Args:
        df: Dataframe with xrd_patterns column
        sample_number: Sample number to retrieve
        pattern_index: Index of pattern if multiple XRD files exist (default: 0)
        
    Returns:
        XRD pattern dict or None if not found
    """
    row = df[df['sample_number'] == sample_number]
    
    if len(row) == 0:
        print(f"Sample {sample_number} not found")
        return None
    
    patterns = row.iloc[0]['xrd_patterns']
    
    if not patterns or len(patterns) == 0:
        print(f"No XRD patterns found for sample {sample_number}")
        return None
    
    if pattern_index >= len(patterns):
        print(f"Pattern index {pattern_index} out of range (found {len(patterns)} patterns)")
        return None
    
    return patterns[pattern_index]
