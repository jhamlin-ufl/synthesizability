# src/synthesizability/disorder.py
"""
Disorder probability prediction for chemical formulas.

Uses the Jakob et al. RNN classifier to predict substitutional disorder probability.
"""

import sys
import torch
import json
import pandas as pd
import warnings
from pathlib import Path
from pymatgen.core import Composition

from .disorder_core import RepresentationGenerator, RNNDisorderClassifier

warnings.filterwarnings('ignore')


def formula_to_composition_string(formula):
    """
    Convert formula to spaced format expected by representation generator.
    
    Args:
        formula: String like 'Mo4Tc2' or 'MoTiTa2'
        
    Returns:
        String like 'Mo4 Tc2' or None if parsing fails
    """
    try:
        comp = Composition(formula)
        
        # Build composition string with explicit stoichiometry
        elements = []
        for element, amount in comp.items():
            elements.append(f"{element}{amount}")
        
        return ' '.join(elements)
    
    except Exception:
        return None


def load_model_and_config():
    """
    Load the pre-trained disorder classifier model and its configuration.
    
    Returns:
        tuple: (config_dict, rep_kwargs, dim, model_kwargs)
    """
    model_dir = Path(__file__).parent.parent.parent / 'data' / 'external' / 'disorder_model'
    
    config_path = model_dir / 'hyperopt_config.json'
    model_path = model_dir / 'hyperopt_best_model.pt'
    
    # Load configuration
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Extract relevant parameters
    rep_kws = {k: v for k, v in config.items() if k in ['rep_type', 'embedding']}
    dim = config.get('dim', 2)
    mod_kws = {k: v for k, v in config.items() if k in ['nl', 'nh']}
    
    return config, rep_kws, dim, mod_kws


def predict_disorder(formulas, batch_size=50, verbose=True):
    """
    Predict disorder probabilities for a list of chemical formulas.
    
    IMPORTANT: The model applies sigmoid internally, so we do NOT apply it again.
    
    Args:
        formulas: List of chemical formula strings
        batch_size: Number of formulas to process at once
        verbose: Print progress messages
        
    Returns:
        tuple: (results_list, failed_formulas)
            results_list: List of dicts with 'formula' and 'disorder_probability'
            failed_formulas: List of formulas that failed to process
    """
    import time
    
    # Load model configuration
    config, rep_kws, dim, mod_kws = load_model_and_config()
    
    all_results = []
    failed_formulas = []
    model = None
    
    # Process in batches
    n_batches = (len(formulas) + batch_size - 1) // batch_size
    
    for batch_idx in range(0, len(formulas), batch_size):
        batch_formulas = formulas[batch_idx:batch_idx + batch_size]
        
        if verbose:
            batch_num = batch_idx // batch_size + 1
            print(f"\nBatch {batch_num}/{n_batches}: Processing {len(batch_formulas)} formulas...")
            batch_start = time.time()
        
        # Convert formulas to composition strings
        composition_strings = []
        valid_formulas = []
        
        for i, formula in enumerate(batch_formulas, 1):
            comp_str = formula_to_composition_string(formula)
            if comp_str is not None:
                composition_strings.append(comp_str)
                valid_formulas.append(formula)
                if verbose and i % 10 == 0:
                    print(f"  Parsed {i}/{len(batch_formulas)} formulas...", end='\r')
            else:
                if verbose:
                    print(f"  WARNING: Failed to parse formula: {formula}")
                failed_formulas.append(formula)
        
        if verbose and len(batch_formulas) > 0:
            print(f"  Parsed {len(batch_formulas)}/{len(batch_formulas)} formulas    ")
        
        if not composition_strings:
            continue
        
        try:
            if verbose:
                print(f"  Generating representations...")
            
            # Create dataframe for representation generation
            temp_df = pd.DataFrame({
                'composition': composition_strings,
                'disordered': [-1] * len(composition_strings)  # Dummy labels
            })
            
            # Generate representations
            gen = RepresentationGenerator(temp_df, dim=dim)
            X, y = gen.get_representations(**rep_kws)
            
            if verbose:
                print(f"  Running neural network prediction...")
            
            # Initialize model on first batch
            if model is None:
                model = RNNDisorderClassifier(
                    nin=X.shape[-1],
                    nout=1,
                    batched=True,
                    **mod_kws
                )
                model_path = Path(__file__).parent.parent.parent / 'data' / 'external' / 'disorder_model' / 'hyperopt_best_model.pt'
                state_dict = torch.load(model_path, map_location='cpu')
                model.load_state_dict(state_dict)
                model.eval()
                
                if verbose:
                    print(f"  Model loaded: input_dim={X.shape[-1]}, hidden={mod_kws['nh']}, layers={mod_kws['nl']}")
            
            # Predict (model already applies sigmoid internally)
            with torch.no_grad():
                probabilities = model(X).numpy().flatten()
            
            # Store results
            for formula, prob in zip(valid_formulas, probabilities):
                all_results.append({
                    'formula': formula,
                    'disorder_probability': float(prob)
                })
            
            if verbose:
                batch_time = time.time() - batch_start
                print(f"  ✓ Batch completed in {batch_time:.1f}s ({batch_time/len(valid_formulas):.2f}s per formula)")
        
        except Exception as e:
            if verbose:
                print(f"  ERROR in batch {batch_num}: {str(e)}")
            failed_formulas.extend(valid_formulas)
    
    return all_results, failed_formulas


def predict_disorder_single(formula):
    """
    Convenience function to predict disorder for a single formula.
    
    Args:
        formula: Chemical formula string
        
    Returns:
        float: Disorder probability or None if prediction fails
    """
    results, failed = predict_disorder([formula], verbose=False)
    
    if results:
        return results[0]['disorder_probability']
    else:
        return None
