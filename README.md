# Synthesizability

Analysis of the gap between DFT-predicted stable crystal structures and experimentally synthesizable materials.

## Project Structure

- `data/` - Experimental and reference data
  - `raw/` - Raw experimental data (XRD, magnetic susceptibility)
  - `processed/` - Cleaned/processed data
  - `external/` - Third-party reference databases
- `scripts/` - Analysis scripts (see `scripts/README.md`)
- `src/synthesizability/` - Reusable Python modules
- `notebooks/` - Jupyter notebooks for analysis and figures
- `figures/` - Generated figures for publication
- `results/` - Analysis results

## Setup

```bash
# Clone repository
git clone https://github.com/jhamlin-ufl/synthesizability.git
cd synthesizability

# Install dependencies
poetry install
```

## Usage

See `scripts/README.md` for information on available analysis scripts.

## Citation

[Paper citation will go here]

## License

[License info]
