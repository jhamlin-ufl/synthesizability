# scripts/compute_disorder_probabilities.py
"""
Compute disorder probabilities with caching.
"""

import sys
from pathlib import Path
import pandas as pd
import time

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from synthesizability.disorder import predict_disorder


def main():
    print("="*80)
    print("DISORDER PROBABILITY COMPUTATION WITH CACHING")
    print("="*80)
    
    overall_start = time.time()
    
    # Paths
    formulas_txt = Path("data/processed/formulas.txt")
    cache_csv = Path("data/processed/disorder_cache.csv")

    # Load formula list
    if not formulas_txt.exists():
        print(f"\nERROR: {formulas_txt} not found!")
        print("Run 'snakemake data/processed/formulas.txt' first")
        return 1

    print(f"\nLoading formulas from {formulas_txt}...")
    formulas_needed = [f.strip() for f in formulas_txt.read_text().splitlines() if f.strip()]
    print(f"Found {len(formulas_needed)} unique formulas")
    
    # Load existing cache
    if cache_csv.exists():
        print(f"\nLoading existing cache from {cache_csv}...")
        cache_df = pd.read_csv(cache_csv)
        print(f"Cache contains {len(cache_df)} formulas")
        
        # Find formulas not in cache
        cached_formulas = set(cache_df['formula'])
        new_formulas = [f for f in formulas_needed if f not in cached_formulas]
    else:
        print(f"\nNo existing cache found at {cache_csv}")
        cache_df = pd.DataFrame(columns=['formula', 'disorder_probability'])
        new_formulas = formulas_needed
    
    if not new_formulas:
        print("\n✓ All formulas already in cache! No computation needed.")
        print(f"✓ Cache saved at: {cache_csv}")
        return 0
    
    # Compute disorder for new formulas
    print(f"\n{'='*80}")
    print(f"COMPUTING DISORDER FOR {len(new_formulas)} NEW FORMULAS")
    print(f"{'='*80}")
    print(f"\nEstimated time: ~{len(new_formulas) * 2:.0f}-{len(new_formulas) * 4:.0f} seconds\n")
    
    compute_start = time.time()
    results, failed = predict_disorder(new_formulas, batch_size=50, verbose=True)
    compute_time = time.time() - compute_start
    
    # Report results
    print(f"\n{'='*80}")
    print("COMPUTATION SUMMARY")
    print(f"{'='*80}")
    print(f"Successful: {len(results)}")
    print(f"Failed: {len(failed)}")
    print(f"Computation time: {compute_time:.1f} seconds ({compute_time/60:.1f} minutes)")
    if len(new_formulas) > 0:
        print(f"Average time per formula: {compute_time/len(new_formulas):.2f} seconds")
    
    if failed:
        print(f"\nFailed formulas:")
        for f in failed:
            print(f"  - {f}")
    
    # Add new results to cache
    if results:
        new_results_df = pd.DataFrame(results)
        cache_df = pd.concat([cache_df, new_results_df], ignore_index=True)
        
        # Save updated cache
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        cache_df.to_csv(cache_csv, index=False)
        print(f"\n✓ Updated cache saved: {cache_csv}")
        print(f"  Total formulas in cache: {len(cache_df)}")
    
    total_time = time.time() - overall_start
    print(f"\n✓ Total runtime: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
