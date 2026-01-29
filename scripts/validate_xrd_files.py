#!/usr/bin/env python3
"""
Find XRD pattern with all 2theta = 0.
"""

from pathlib import Path
import pandas as pd
import numpy as np


def find_constant_twotheta():
    """Find patterns where 2theta doesn't vary."""
    
    pickle_path = Path.cwd() / "data" / "processed" / "synthesis_data.pkl"
    
    if not pickle_path.exists():
        print(f"Error: File not found: {pickle_path}")
        return
    
    df = pd.read_pickle(pickle_path)
    
    print("="*80)
    print("LOOKING FOR CONSTANT 2THETA PATTERNS")
    print("="*80)
    
    for idx, row in df.iterrows():
        patterns = row['xrd_patterns']
        
        if not patterns or len(patterns) == 0:
            continue
        
        sample_num = row['sample_number']
        formula = row['formula']
        
        for i, pattern in enumerate(patterns):
            two_theta = pattern['two_theta']
            intensity = pattern['intensity']
            
            if len(two_theta) == 0:
                continue
            
            # Check if all 2theta values are the same
            unique_2theta = np.unique(two_theta)
            
            if len(unique_2theta) == 1:
                print(f"\n⚠⚠ FOUND CONSTANT 2THETA!")
                print(f"Sample: {sample_num:04d} - {formula}")
                print(f"Pattern index: {i}")
                print(f"Instrument: {pattern['instrument']}")
                print(f"All {len(two_theta)} points have 2theta = {unique_2theta[0]}")
                print(f"Intensity range: {np.min(intensity):.2f} to {np.max(intensity):.2f}")
                print(f"First 10 intensities: {intensity[:10]}")
                print(f"Last 10 intensities: {intensity[-10:]}")
                
                # Find the source file
                sample_dir = Path.cwd() / "data" / "raw" / f"{sample_num:04d}_{formula}"
                if sample_dir.exists():
                    xrd_files = list(sample_dir.glob("*.xy")) + list(sample_dir.glob("*.txt"))
                    xrd_files = [f for f in xrd_files if 'chi' not in f.name.lower()]
                    if i < len(xrd_files):
                        print(f"Source file: {xrd_files[i].name}")


def main():
    find_constant_twotheta()


if __name__ == "__main__":
    main()