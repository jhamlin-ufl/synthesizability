# scripts/query_oqmd_hulls.py
"""
Query OQMD database for hull distances for all compositions in synthesis data.
"""
import pandas as pd
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from synthesizability.oqmd import (
    parse_formula_to_oqmd,
    query_formation_energies,
    check_database_exists
)

def main():
    # Check database
    if not check_database_exists():
        print("ERROR: OQMD database not found!")
        print("Please run: poetry run python scripts/validate_oqmd_database.py")
        sys.exit(1)
    
    # Load synthesis data
    print("Loading synthesis data...")
    df = pd.read_csv('data/processed/synthesis_data_no_disorder.csv')
    
    results = []
    for idx, row in df.iterrows():
        formula = row['formula']
        oqmd_formula = parse_formula_to_oqmd(formula)
        
        print(f"Querying {formula} -> {oqmd_formula} ({idx+1}/{len(df)})...")
        
        entries = query_formation_energies(oqmd_formula)
        
        if entries:
            # Take minimum stability (most stable polymorph)
            min_entry = min(entries, key=lambda x: x['stability'] if x['stability'] is not None else float('inf'))
            results.append({
                'formula': formula,
                'oqmd_formula': oqmd_formula,
                'oqmd_delta_e': min_entry['delta_e'],
                'oqmd_stability': min_entry['stability'],
                'oqmd_entry_id': min_entry['entry_id'],
                'oqmd_n_polymorphs': len(entries)
            })
        else:
            results.append({
                'formula': formula,
                'oqmd_formula': oqmd_formula,
                'oqmd_delta_e': None,
                'oqmd_stability': None,
                'oqmd_entry_id': None,
                'oqmd_n_polymorphs': 0
            })
    
    # Save results
    results_df = pd.DataFrame(results)
    output_path = Path('data/processed/oqmd_hull_data.csv')
    results_df.to_csv(output_path, index=False)
    print(f"\nSaved results to {output_path}")
    print(f"Found OQMD data for {results_df['oqmd_stability'].notna().sum()}/{len(results_df)} compositions")

if __name__ == "__main__":
    main()