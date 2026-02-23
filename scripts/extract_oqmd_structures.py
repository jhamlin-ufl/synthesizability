# scripts/extract_oqmd_structures.py
"""
Extract CIF structure files from OQMD database for compositions with hull data.
"""
import pandas as pd
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from synthesizability.oqmd import (
    get_structure_from_db,
    check_database_exists
)
from pymatgen.io.cif import CifWriter

def main():
    # Check database
    if not check_database_exists():
        print("ERROR: OQMD database not found!")
        print("Please run: poetry run python scripts/validate_oqmd_database.py")
        sys.exit(1)
    
    # Load OQMD hull data
    print("Loading OQMD hull data...")
    df = pd.read_csv('data/processed/oqmd_hull_data.csv')
    
    # Filter to entries with data
    df_valid = df[df['oqmd_entry_id'].notna()].copy()
    print(f"Found {len(df_valid)} compositions with OQMD data")
    
    # Create output directory
    output_dir = Path('data/external/oqmd_structures')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract structures
    extracted = 0
    failed = []
    
    for idx, row in df_valid.iterrows():
        formula = row['formula']
        entry_id = int(row['oqmd_entry_id'])
        
        print(f"Extracting {formula} (entry {entry_id}) [{idx+1}/{len(df_valid)}]...")
        
        # Get structure from database
        structure = get_structure_from_db(entry_id)
        
        if structure is None:
            print(f"  ⚠ Failed to extract structure for {formula}")
            failed.append(formula)
            continue
        
        # Create formula-specific directory
        formula_dir = output_dir / formula
        formula_dir.mkdir(exist_ok=True)
        
        # Write CIF file
        cif_path = formula_dir / f"{entry_id}.cif"
        writer = CifWriter(structure)
        writer.write_file(str(cif_path))
        
        extracted += 1
        print(f"  ✓ Saved to {cif_path}")
    
    print(f"\n{'='*60}")
    print(f"Extraction complete:")
    print(f"  Successfully extracted: {extracted}/{len(df_valid)}")
    print(f"  Failed: {len(failed)}")
    if failed:
        print(f"  Failed compositions: {', '.join(failed)}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()