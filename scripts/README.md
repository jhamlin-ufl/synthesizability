# Analysis Scripts

## find_relevant_SCs.py

Search the SuperCon database for superconducting materials containing only specified elements.

### Usage

```bash
# Find superconductors containing only Mo, Nb, and Ta
poetry run python scripts/find_relevant_SCs.py Mo Nb Ta

# Include oxygen-containing compounds
poetry run python scripts/find_relevant_SCs.py --include-oxygen La Ba Cu

# See all options
poetry run python scripts/find_relevant_SCs.py --help
```

### Output

The script returns a table with:
- **Material**: Common name or chemical formula
- **Formula**: Detailed chemical composition
- **Elements**: Elements present (excluding oxygen by default)
- **Tc (K)**: Superconducting critical temperature
- **Reference**: Journal citation with clickable link
  - Direct DOI links for Physical Review, J.Phys.Soc.Jpn
  - Google Scholar search links for other journals

### Notes

- Tc values may include measurements under pressure - consult original references
- The script automatically downloads the SuperCon database on first run (~4 MB)
- Data is cached in `data/external/supercon/`

### Data Source

MDR SuperCon Datasheet Ver.240322 from the National Institute for Materials Science (NIMS)
- DOI: https://doi.org/10.48505/nims.3837
- License: CC BY 4.0
- See `data/external/supercon/README.md` for details

## Other Scripts

[Documentation for other scripts will go here]
