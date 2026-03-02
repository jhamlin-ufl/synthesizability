# scripts/build_dataframe.py
#!/usr/bin/env python3
"""
Build dataframe from synthesis experiment data.
"""

from pathlib import Path
import pandas as pd
from synthesizability.io import build_dataframe, analyze_field_statistics, show_missing_samples


def main():
    # Get the data/raw directory
    data_raw_dir = Path.cwd() / "data" / "raw"
    
    if not data_raw_dir.exists():
        print(f"Error: Directory not found: {data_raw_dir}")
        return
    
    print("Building dataframe...")
    df = build_dataframe(data_raw_dir)
    
    print(f"Created dataframe with {len(df)} samples and {len(df.columns)} columns")
    
    # Merge OQMD hull data if it exists
    oqmd_hull_path = Path.cwd() / "data" / "processed" / "oqmd_hull_data.csv"
    if oqmd_hull_path.exists():
        print("\nMerging OQMD hull data...")
        oqmd_df = pd.read_csv(oqmd_hull_path)
        df = df.merge(oqmd_df, on='formula', how='left')
        print(f"✓ Merged OQMD data for {df['oqmd_stability'].notna().sum()}/{len(df)} compositions")
    else:
        print("\n⚠ OQMD hull data not found, skipping merge")
    
    # Analyze statistics
    analyze_field_statistics(df)
    
    # Show missing data
    show_missing_samples(df)
    
    # Save to CSV (excluding the large fields)
    output_path = Path.cwd() / "data" / "processed" / "synthesis_data.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save a version without the large fields
    df_to_save = df.drop(columns=['files', 'status_content', 'synthesis_content', 'xrd_patterns',
                               'composition_measured_fractions', 'composition_expected_fractions'])
    df_to_save.to_csv(output_path, index=False)
    print(f"\n✓ Saved dataframe to: {output_path}")
    
    # Also save full dataframe as pickle
    pickle_path = Path.cwd() / "data" / "processed" / "synthesis_data.pkl"
    df.to_pickle(pickle_path)
    print(f"✓ Saved full dataframe (with text fields and XRD patterns) to: {pickle_path}")
    
    return df


if __name__ == "__main__":
    df = main()
