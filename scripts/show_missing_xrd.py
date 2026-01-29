#!/usr/bin/env python3
"""
Show which samples are missing XRD data.
"""

from pathlib import Path
import pandas as pd


def show_missing_xrd(df: pd.DataFrame):
    """Show samples without XRD data."""
    
    print("="*80)
    print("SAMPLES WITHOUT XRD DATA")
    print("="*80)
    print(f"\nTotal samples: {len(df)}")
    
    # Find samples with and without XRD
    has_xrd = df['xrd_n_files'] > 0
    no_xrd = ~has_xrd
    
    n_with_xrd = has_xrd.sum()
    n_without_xrd = no_xrd.sum()
    
    print(f"Samples with XRD: {n_with_xrd} ({n_with_xrd/len(df)*100:.1f}%)")
    print(f"Samples without XRD: {n_without_xrd} ({n_without_xrd/len(df)*100:.1f}%)")
    
    # Show samples without XRD
    print("\n" + "="*80)
    print("SAMPLES MISSING XRD DATA:")
    print("="*80)
    
    missing_df = df[no_xrd][['sample_number', 'formula', 'superconductivity', 
                              'xrd_result', 'prediction_list']].copy()
    
    # Group by superconductivity status
    print("\n### Superconducting samples without XRD:")
    print("-"*80)
    sc_samples = missing_df[missing_df['tc_kelvin'].notna()] if 'tc_kelvin' in missing_df.columns else pd.DataFrame()
    
    if len(sc_samples) == 0:
        # Check for any "Tc of" in superconductivity text
        sc_samples = missing_df[missing_df['superconductivity'].str.contains('Tc of', na=False)]
    
    if len(sc_samples) > 0:
        for _, row in sc_samples.iterrows():
            print(f"  {row['sample_number']:04d} - {row['formula']}")
            print(f"       SC: {row['superconductivity']}")
            print(f"       XRD result: {row['xrd_result']}")
            print()
    else:
        print("  (None)")
    
    print("\n### Non-superconducting samples without XRD:")
    print("-"*80)
    non_sc = missing_df[~missing_df.index.isin(sc_samples.index)]
    
    # Sort by prediction list
    for pred_list in non_sc['prediction_list'].dropna().unique():
        samples = non_sc[non_sc['prediction_list'] == pred_list]
        if len(samples) > 0:
            print(f"\n  {pred_list}:")
            for _, row in samples.iterrows():
                print(f"    {row['sample_number']:04d} - {row['formula']}")
                if pd.notna(row['superconductivity']):
                    print(f"         SC: {row['superconductivity']}")
    
    # Samples with no prediction list
    no_pred = non_sc[non_sc['prediction_list'].isna()]
    if len(no_pred) > 0:
        print(f"\n  (No prediction list):")
        for _, row in no_pred.iterrows():
            print(f"    {row['sample_number']:04d} - {row['formula']}")
    
    print("\n" + "="*80)


def show_xrd_summary(df: pd.DataFrame):
    """Show summary of XRD data availability."""
    
    print("\n" + "="*80)
    print("XRD DATA SUMMARY")
    print("="*80)
    
    # Count by instrument
    siemens_count = df['has_siemens_xrd'].sum()
    panalytical_count = df['has_panalytical_xrd'].sum()
    both_count = (df['has_siemens_xrd'] & df['has_panalytical_xrd']).sum()
    
    print(f"\nSiemens D500 (.txt) files: {siemens_count}")
    print(f"Panalytical (.xy) files: {panalytical_count}")
    print(f"Samples with both: {both_count}")
    
    # XRD result status
    print("\n" + "-"*80)
    print("XRD Analysis Status:")
    print("-"*80)
    
    has_xrd = df['xrd_n_files'] > 0
    xrd_samples = df[has_xrd]
    
    if len(xrd_samples) > 0:
        result_counts = xrd_samples['xrd_result'].value_counts()
        for result, count in result_counts.items():
            pct = (count / len(xrd_samples)) * 100
            print(f"  {result}: {count} ({pct:.1f}%)")
    
    print("\n" + "="*80)


def main():
    # Load the dataframe
    pickle_path = Path.cwd() / "data" / "processed" / "synthesis_data.pkl"
    
    if not pickle_path.exists():
        print(f"Error: File not found: {pickle_path}")
        print("Run scripts/build_dataframe.py first")
        return
    
    df = pd.read_pickle(pickle_path)
    
    # Show XRD summary
    show_xrd_summary(df)
    
    # Show missing XRD
    show_missing_xrd(df)


if __name__ == "__main__":
    main()