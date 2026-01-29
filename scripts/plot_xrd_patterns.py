#!/usr/bin/env python3
"""
Plot all XRD patterns on a single figure.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_all_xrd_patterns(df: pd.DataFrame, output_path: Path = None, 
                          offset: float = 0, normalize: bool = False):
    """
    Plot all XRD patterns on a single figure.
    
    Args:
        df: Dataframe with xrd_patterns column
        output_path: Optional path to save figure
        offset: Vertical offset between patterns (if 0, plots overlap)
        normalize: Whether to normalize intensities to [0, 1] for each pattern
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    n_patterns = 0
    
    # Iterate through samples
    for idx, row in df.iterrows():
        patterns = row['xrd_patterns']
        
        if not patterns or len(patterns) == 0:
            continue
        
        sample_num = row['sample_number']
        formula = row['formula']
        
        # Plot each pattern for this sample
        for i, pattern in enumerate(patterns):
            two_theta = pattern['two_theta']
            intensity = pattern['intensity']
            
            if len(two_theta) == 0:
                continue
            
            # Normalize if requested
            if normalize and intensity.max() > 0:
                intensity = intensity / intensity.max()
            
            # Apply offset
            intensity_plot = intensity + (n_patterns * offset)
            
            # Create label
            label = f"{sample_num:04d} - {formula}"
            if len(patterns) > 1:
                label += f" ({i+1})"
            
            # Plot
            ax.plot(two_theta, intensity_plot, linewidth=0.5, alpha=0.7, label=label)
            
            n_patterns += 1
    
    ax.set_xlabel('2θ (degrees)', fontsize=12)
    ax.set_ylabel('Intensity (a.u.)', fontsize=12)
    ax.set_title(f'All XRD Patterns (n={n_patterns})', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Only show legend if reasonable number of patterns
    if n_patterns <= 30:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        plt.tight_layout()
    else:
        plt.tight_layout()
        print(f"Too many patterns ({n_patterns}) for legend")
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved figure to: {output_path}")
    
    plt.show()


def plot_xrd_patterns_stacked(df: pd.DataFrame, output_path: Path = None):
    """
    Plot XRD patterns stacked vertically with offsets.
    """
    # Calculate appropriate offset based on max intensity
    all_intensities = []
    for _, row in df.iterrows():
        if row['xrd_patterns']:
            for pattern in row['xrd_patterns']:
                if len(pattern['intensity']) > 0:
                    all_intensities.append(pattern['intensity'].max())
    
    if all_intensities:
        max_intensity = max(all_intensities)
        offset = max_intensity * 1.2  # 20% spacing between patterns
    else:
        offset = 1000
    
    plot_all_xrd_patterns(df, output_path=output_path, offset=offset, normalize=False)


def plot_xrd_patterns_overlaid(df: pd.DataFrame, output_path: Path = None):
    """
    Plot XRD patterns overlaid (no offset), normalized.
    """
    plot_all_xrd_patterns(df, output_path=output_path, offset=0, normalize=True)


def main():
    # Load the dataframe
    pickle_path = Path.cwd() / "data" / "processed" / "synthesis_data.pkl"
    
    if not pickle_path.exists():
        print(f"Error: File not found: {pickle_path}")
        print("Run scripts/build_dataframe.py first")
        return
    
    df = pd.read_pickle(pickle_path)
    
    # Count patterns
    n_samples_with_xrd = sum(1 for _, row in df.iterrows() if row['xrd_patterns'])
    total_patterns = sum(len(row['xrd_patterns']) for _, row in df.iterrows() if row['xrd_patterns'])
    
    print(f"Found {n_samples_with_xrd} samples with XRD data")
    print(f"Total XRD patterns: {total_patterns}")
    print()
    
    # Create output directory
    figures_dir = Path.cwd() / "figures"
    figures_dir.mkdir(exist_ok=True)
    
    # Plot stacked
    print("Creating stacked plot...")
    plot_xrd_patterns_stacked(df, output_path=figures_dir / "xrd_all_stacked.png")
    
    # Plot overlaid
    print("\nCreating overlaid plot...")
    plot_xrd_patterns_overlaid(df, output_path=figures_dir / "xrd_all_overlaid.png")


if __name__ == "__main__":
    main()