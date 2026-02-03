# scripts/analyze_susceptibility.py
"""
Analyze AC susceptibility data for all samples.

Generates:
- susceptibility_real_part.pdf
- susceptibility_imaginary_part.pdf
- hc2_with_fits.pdf
- hc2_fit_parameters.csv
"""

from pathlib import Path
import sys
import pandas as pd

sys.path.insert(0, "src")

from synthesizability.io import build_dataframe
from synthesizability.susceptibility import (
    load_all_chi_data,
    extract_tc_values,
    fit_hc2_models,
    plot_chi_real_grid,
    plot_chi_imaginary_grid,
    plot_hc2_grid
)


def main():
    print("="*80)
    print("AC Susceptibility Analysis")
    print("="*80)
    
    # Load dataframe to get samples with chi data
    data_dir = Path("data/raw")
    df = build_dataframe(data_dir)
    
    samples_with_chi = df[df['chi_n_files'] > 0].sort_values('sample_number')
    
    print(f"\nFound {len(samples_with_chi)} samples with chi data")
    
    # Load chi data for all samples
    print("\nLoading chi data...")
    samples_chi_data = {}
    samples_tc_data = {}
    
    for _, row in samples_with_chi.iterrows():
        sample_id = row['sample_id']
        sample_dir = data_dir / sample_id
        
        print(f"  Loading {sample_id}...")
        chi_data = load_all_chi_data(sample_dir)
        
        if len(chi_data) > 0:
            samples_chi_data[sample_id] = chi_data
            
            # Extract Tc values
            tc_data = extract_tc_values(chi_data)
            if len(tc_data) > 0:
                samples_tc_data[sample_id] = tc_data
    
    print(f"\nSuccessfully loaded chi data for {len(samples_chi_data)} samples")
    print(f"Extracted Tc data for {len(samples_tc_data)} samples")
    
    # Sample order for plots
    sample_order = [row['sample_id'] for _, row in samples_with_chi.iterrows()]
    
    # Create output directory
    output_dir = Path("results/susceptibility")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Plot real part
    print("\nGenerating real part plot...")
    fig_real = plot_chi_real_grid(samples_chi_data, sample_order)
    fig_real.savefig(output_dir / 'susceptibility_real_part.pdf', 
                     dpi=200, bbox_inches='tight')
    print(f"  Saved to {output_dir / 'susceptibility_real_part.pdf'}")
    
    # Plot imaginary part
    print("\nGenerating imaginary part plot...")
    fig_imag = plot_chi_imaginary_grid(samples_chi_data, sample_order)
    fig_imag.savefig(output_dir / 'susceptibility_imaginary_part.pdf',
                     dpi=200, bbox_inches='tight')
    print(f"  Saved to {output_dir / 'susceptibility_imaginary_part.pdf'}")
    
    # Plot Hc2 fits
    print("\nGenerating Hc2 fits plot...")
    fig_hc2, fit_results = plot_hc2_grid(samples_tc_data, sample_order)
    fig_hc2.savefig(output_dir / 'hc2_with_fits.pdf',
                    dpi=200, bbox_inches='tight')
    print(f"  Saved to {output_dir / 'hc2_with_fits.pdf'}")
    
    # Save fit parameters
    print("\nSaving fit parameters...")
    summary_data = []
    for composition, results in fit_results.items():
        summary_data.append({
            'Composition': composition,
            'Hc2(0) Linear (T)': results['linear']['Hc2_0'],
            'Tc Linear (K)': results['linear']['Tc'],
            'Hc2(0) Quadratic (T)': results['quadratic']['Hc2_0'],
            'Tc Quadratic (K)': results['quadratic']['Tc']
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_dir / 'hc2_fit_parameters.csv', index=False)
    print(f"  Saved to {output_dir / 'hc2_fit_parameters.csv'}")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
    
    # Print summary table
    print("\nFit parameters summary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()