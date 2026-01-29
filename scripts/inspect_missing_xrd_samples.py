#!/usr/bin/env python3
"""
Inspect files in samples reported as missing XRD data.
"""

from pathlib import Path


def inspect_sample_files(sample_numbers: list, data_raw_dir: Path):
    """Inspect all files in given sample directories."""
    
    for sample_num in sample_numbers:
        # Find matching directory
        matching_dirs = list(data_raw_dir.glob(f"{sample_num:04d}_*"))
        
        if not matching_dirs:
            print(f"\n{sample_num:04d}: Directory not found")
            continue
        
        dir_path = matching_dirs[0]
        print(f"\n{'='*80}")
        print(f"{dir_path.name}")
        print('='*80)
        
        files = sorted(dir_path.iterdir())
        
        if not files:
            print("  (empty directory)")
            continue
        
        for file in files:
            if file.is_file():
                size = file.stat().st_size
                print(f"  {file.name:<40} {size:>10,} bytes")
                
                # For text files, show first few lines
                if file.suffix in ['.txt', '.xy'] or size < 10000:
                    try:
                        with open(file, 'r', encoding='utf-8', errors='replace') as f:
                            first_lines = [f.readline() for _ in range(5)]
                        print("    First lines:")
                        for i, line in enumerate(first_lines, 1):
                            print(f"      {i}: {line.rstrip()[:70]}")
                    except Exception as e:
                        print(f"    (Could not read: {e})")
                print()


def main():
    data_raw_dir = Path.cwd() / "data" / "raw"
    
    if not data_raw_dir.exists():
        print(f"Error: Directory not found: {data_raw_dir}")
        return
    
    # Superconducting samples without XRD
    sc_samples = [449, 451, 454, 455, 457, 462, 463, 464]
    
    # Non-superconducting samples without XRD
    non_sc_samples = [465, 512, 468]
    
    print("="*80)
    print("SUPERCONDUCTING SAMPLES REPORTED AS MISSING XRD")
    print("="*80)
    inspect_sample_files(sc_samples, data_raw_dir)
    
    print("\n\n")
    print("="*80)
    print("NON-SUPERCONDUCTING SAMPLES REPORTED AS MISSING XRD")
    print("="*80)
    inspect_sample_files(non_sc_samples, data_raw_dir)


if __name__ == "__main__":
    main()