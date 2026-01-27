#!/usr/bin/env python3
"""
Inspect XRD files to understand their format and develop parsing strategy.
"""

from pathlib import Path
import random


def inspect_xrd_files(data_raw_dir: Path, n_samples: int = 5):
    """Inspect a random sample of XRD files to understand their format."""
    
    # Find all XRD files
    txt_files = []
    xy_files = []
    
    for dir_path in data_raw_dir.iterdir():
        if not dir_path.is_dir():
            continue
        
        for file in dir_path.iterdir():
            if file.suffix == '.txt' and 'circular' not in file.name.lower() and 'chi' not in file.name.lower():
                txt_files.append(file)
            elif file.suffix == '.xy':
                xy_files.append(file)
    
    print("="*80)
    print(f"Found {len(txt_files)} Siemens D500 .txt files")
    print(f"Found {len(xy_files)} Panalytical .xy files")
    print("="*80)
    
    # Sample Siemens files
    if txt_files:
        print("\n" + "="*80)
        print("SIEMENS D500 .txt FILE SAMPLES")
        print("="*80)
        
        sample_txt = random.sample(txt_files, min(n_samples, len(txt_files)))
        
        for txt_file in sample_txt:
            print(f"\n{'='*80}")
            print(f"File: {txt_file.parent.name}/{txt_file.name}")
            print(f"Size: {txt_file.stat().st_size} bytes")
            print("-"*80)
            
            try:
                with open(txt_file, 'r', encoding='utf-8', errors='replace') as f:
                    lines = f.readlines()
                    print(f"Total lines: {len(lines)}")
                    print("\nFirst 30 lines:")
                    print("-"*80)
                    for i, line in enumerate(lines[:30], 1):
                        print(f"{i:3d}: {line.rstrip()}")
                    
                    if len(lines) > 30:
                        print("\n...")
                        print(f"\nLast 5 lines:")
                        print("-"*80)
                        for i, line in enumerate(lines[-5:], len(lines)-4):
                            print(f"{i:3d}: {line.rstrip()}")
                            
            except Exception as e:
                print(f"Error reading file: {e}")
    
    # Sample Panalytical files
    if xy_files:
        print("\n" + "="*80)
        print("PANALYTICAL .xy FILE SAMPLES")
        print("="*80)
        
        sample_xy = random.sample(xy_files, min(n_samples, len(xy_files)))
        
        for xy_file in sample_xy:
            print(f"\n{'='*80}")
            print(f"File: {xy_file.parent.name}/{xy_file.name}")
            print(f"Size: {xy_file.stat().st_size} bytes")
            print("-"*80)
            
            try:
                with open(xy_file, 'r', encoding='utf-8', errors='replace') as f:
                    lines = f.readlines()
                    print(f"Total lines: {len(lines)}")
                    print("\nFirst 30 lines:")
                    print("-"*80)
                    for i, line in enumerate(lines[:30], 1):
                        print(f"{i:3d}: {line.rstrip()}")
                    
                    if len(lines) > 30:
                        print("\n...")
                        print(f"\nLast 5 lines:")
                        print("-"*80)
                        for i, line in enumerate(lines[-5:], len(lines)-4):
                            print(f"{i:3d}: {line.rstrip()}")
                            
            except Exception as e:
                print(f"Error reading file: {e}")


def main():
    data_raw_dir = Path.cwd() / "data" / "raw"
    
    if not data_raw_dir.exists():
        print(f"Error: Directory not found: {data_raw_dir}")
        return
    
    inspect_xrd_files(data_raw_dir, n_samples=3)


if __name__ == "__main__":
    main()