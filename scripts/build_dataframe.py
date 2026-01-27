#!/usr/bin/env python3
"""
Build a dataframe from synthesis experiment data and analyze field populations.
"""

import re
from pathlib import Path
import pandas as pd
import numpy as np


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


def extract_tc_value(superconductivity_text: str) -> float:
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


def extract_mass_loss(synthesis_content: str) -> tuple:
    """Extract mass loss percentage, initial and final masses."""
    if not synthesis_content:
        return None, None, None
    
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


def parse_status_file(content: str) -> dict:
    """Parse structured data from STATUS file."""
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
    data['tc_kelvin'] = extract_tc_value(data['superconductivity']) if data['superconductivity'] else None
    
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


def parse_synthesis_file(content: str) -> dict:
    """Parse structured data from SYNTHESIS file."""
    if not content:
        return {
            'mass_loss_percent': None,
            'initial_mass_g': None,
            'final_mass_g': None,
            'has_powder_premelting': False,
            'air_sensitive_handling': False
        }
    
    mass_loss, initial_mass, final_mass = extract_mass_loss(content)
    
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
        
        # Parse structured data
        status_data = parse_status_file(status_content)
        synthesis_data = parse_synthesis_file(synthesis_content)
        
        # Check for XRD data files
        has_siemens_xrd = any(f.endswith('.txt') and 'circular' not in f.lower() 
                               and 'chi' not in f.lower() for f in files)
        has_panalytical_xrd = any(f.endswith('.xy') for f in files)
        has_summary = any(f.endswith('.pptx') for f in files)
        
        # Build record
        record = {
            'sample_number': sample_number,
            'formula': formula,
            'files': files,
            'status_content': status_content,
            'synthesis_content': synthesis_content,
            'has_siemens_xrd': has_siemens_xrd,
            'has_panalytical_xrd': has_panalytical_xrd,
            'has_summary': has_summary,
            **status_data,
            **synthesis_data
        }
        
        records.append(record)
    
    return pd.DataFrame(records)


def analyze_field_statistics(df: pd.DataFrame):
    """Print statistics for each field in the dataframe."""
    
    print("\n" + "="*80)
    print("DATAFRAME FIELD STATISTICS")
    print("="*80)
    print(f"\nTotal samples: {len(df)}\n")
    
    # Analyze each column
    for col in df.columns:
        if col in ['files', 'status_content', 'synthesis_content']:
            # Skip the large text fields
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
        if col in ['files', 'status_content', 'synthesis_content']:
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


def main():
    # Get the data/raw directory
    data_raw_dir = Path.cwd() / "data" / "raw"
    
    if not data_raw_dir.exists():
        print(f"Error: Directory not found: {data_raw_dir}")
        return
    
    print("Building dataframe...")
    df = build_dataframe(data_raw_dir)
    
    print(f"Created dataframe with {len(df)} samples and {len(df.columns)} columns")
    
    # Analyze statistics
    analyze_field_statistics(df)
    
    # Show missing data
    show_missing_samples(df)
    
    # Save to CSV (excluding the large text fields for now)
    output_path = Path.cwd() / "data" / "processed" / "synthesis_data.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save a version without the large text fields
    df_to_save = df.drop(columns=['files', 'status_content', 'synthesis_content'])
    df_to_save.to_csv(output_path, index=False)
    print(f"\n✓ Saved dataframe to: {output_path}")
    
    # Also save full dataframe as pickle
    pickle_path = Path.cwd() / "data" / "processed" / "synthesis_data.pkl"
    df.to_pickle(pickle_path)
    print(f"✓ Saved full dataframe (with text fields) to: {pickle_path}")
    
    return df


if __name__ == "__main__":
    df = main()