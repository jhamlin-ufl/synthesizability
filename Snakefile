# Snakefile for synthesizability pipeline
from pathlib import Path

# Dynamically find all Python files in src
SRC_FILES = [str(p) for p in Path("src/synthesizability").rglob("*.py")]

# Track all files in data/raw (not just directories)
RAW_DATA_FILES = [str(p) for p in Path("data/raw").rglob("*") if p.is_file()]

# Track reference data files
REFERENCE_DATA_FILES = [
    "data/external/reference/element_prices.csv",
    "data/external/reference/element_vapor_pressures.csv"
]

print(f"Tracking {len(SRC_FILES)} source files")
print(f"Tracking {len(RAW_DATA_FILES)} raw data files")
print(f"Tracking {len(REFERENCE_DATA_FILES)} reference data files")

rule all:
    input:
        "data/processed/synthesis_data.csv",
        "data/processed/synthesis_data.pkl",
        "results/susceptibility/susceptibility_real_part.pdf",
        "results/susceptibility/susceptibility_imaginary_part.pdf",
        "results/susceptibility/hc2_with_fits.pdf",
        "results/susceptibility/hc2_fit_parameters.csv"

rule build_dataframe:
    input:
        script="scripts/build_dataframe.py",
        src=SRC_FILES,
        data=RAW_DATA_FILES,
        reference=REFERENCE_DATA_FILES
    output:
        csv="data/processed/synthesis_data.csv",
        pkl="data/processed/synthesis_data.pkl"
    shell:
        "poetry run python {input.script}"

rule analyze_susceptibility:
    input:
        script="scripts/analyze_susceptibility.py",
        src=SRC_FILES,
        data=RAW_DATA_FILES
    output:
        real="results/susceptibility/susceptibility_real_part.pdf",
        imag="results/susceptibility/susceptibility_imaginary_part.pdf",
        hc2="results/susceptibility/hc2_with_fits.pdf",
        params="results/susceptibility/hc2_fit_parameters.csv"
    shell:
        "poetry run python {input.script}"