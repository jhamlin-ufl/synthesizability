# test_disorder_import.py
"""
Quick test to verify disorder model files are accessible.
Run from synthesizability root: python test_disorder_import.py
"""

from pathlib import Path
import json
import torch

def test_model_files():
    print("Testing disorder model files...")
    
    model_dir = Path("data/external/disorder_model")
    
    # Check model file
    model_path = model_dir / "hyperopt_best_model.pt"
    assert model_path.exists(), f"Model file not found: {model_path}"
    print(f"✓ Found model: {model_path}")
    
    # Load model to verify it's valid
    state_dict = torch.load(model_path, map_location='cpu')
    print(f"✓ Model loads successfully, {len(state_dict)} parameters")
    
    # Check config file
    config_path = model_dir / "hyperopt_config.json"
    assert config_path.exists(), f"Config file not found: {config_path}"
    print(f"✓ Found config: {config_path}")
    
    # Load config
    with open(config_path) as f:
        config = json.load(f)
    print(f"✓ Config keys: {list(config.keys())}")
    print(f"  dim={config['dim']}, nl={config['nl']}, nh={config['nh']}")
    
    print("\n✓ All model files verified!")

if __name__ == "__main__":
    test_model_files()