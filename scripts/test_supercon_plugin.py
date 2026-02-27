# scripts/test_supercon_plugin.py
import sys, json
from pathlib import Path
sys.path.insert(0, "src")
from synthesizability.dashboard_plugins.supercon import (
    _composition_distance, _n_elements, _load_hits
)

# HfTa4Zr has measured Tc=7.3K, good test case
hits = _load_hits("HfTa4Zr")
print(f"HfTa4Zr: {len(hits)} hits")
for h in hits[:5]:
    dist = _composition_distance("HfTa4Zr", h["formula"])
    print(f"  dist={dist:.3f}  n_elem={_n_elements(h['formula'])}  {h['formula']}  Tc={h['tc']}")