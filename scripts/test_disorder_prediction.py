# scripts/test_disorder_prediction.py
"""
Test disorder prediction on a few sample formulas from your dataset.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from synthesizability.disorder import predict_disorder, predict_disorder_single


def main():
    print("="*80)
    print("TESTING DISORDER PREDICTION")
    print("="*80)
    
    # Test formulas from your dataset
    test_formulas = [
        'MoTiTa2',
        'MoTi2Ta', 
        'HfTa4Zr',
        'HfMoTa2',
        'MoTaZr2'
    ]
    
    print(f"\nTesting {len(test_formulas)} formulas...\n")
    
    # Batch prediction
    results, failed = predict_disorder(test_formulas, batch_size=10, verbose=True)
    
    print("\n" + "="*80)
    print("RESULTS")
    print("="*80)
    
    if results:
        print("\nSuccessful predictions:")
        print(f"{'Formula':<15} {'Disorder Probability':>20}")
        print("-"*40)
        for r in results:
            print(f"{r['formula']:<15} {r['disorder_probability']:>20.4f}")
    
    if failed:
        print(f"\nFailed formulas: {failed}")
    
    # Test single prediction
    print("\n" + "="*80)
    print("TESTING SINGLE PREDICTION")
    print("="*80)
    
    formula = 'NbTa2Zr'
    prob = predict_disorder_single(formula)
    if prob is not None:
        print(f"\n{formula}: {prob:.4f}")
    else:
        print(f"\n{formula}: FAILED")


if __name__ == "__main__":
    main()