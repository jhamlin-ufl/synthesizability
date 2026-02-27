# Snakefile for synthesizability pipeline
from pathlib import Path

# Dynamically find all Python files in src
SRC_FILES = [str(p) for p in Path("src/synthesizability").rglob("*.py")]

# Track all files in data/raw
RAW_DATA_FILES = [str(p) for p in Path("data/raw").rglob("*") if p.is_file()]

# Track specific file types for targeted dependencies
CHI_DATA_FILES = [str(p) for p in Path("data/raw").rglob("*chiAC*.txt")]
STATUS_FILES = [str(p) for p in Path("data/raw").rglob("STATUS")]
SYNTHESIS_FILES = [str(p) for p in Path("data/raw").rglob("SYNTHESIS")]
XRD_FILES = [str(p) for p in Path("data/raw").rglob("*.xy")] + \
            [str(p) for p in Path("data/raw").rglob("*.txt") if "chiAC" not in p.name]

# Core data files for dataframe building
DATAFRAME_INPUT_FILES = STATUS_FILES + SYNTHESIS_FILES + XRD_FILES

# Track reference data files
REFERENCE_DATA_FILES = [
    "data/external/reference/element_prices.csv",
    "data/external/reference/element_vapor_pressures.csv"
]

# Track disorder model files
DISORDER_MODEL_FILES = [
    "data/external/disorder_model/hyperopt_best_model.pt",
    "data/external/disorder_model/hyperopt_config.json"
]

print(f"Tracking {len(SRC_FILES)} source files")
print(f"Tracking {len(RAW_DATA_FILES)} raw data files")
print(f"  - {len(CHI_DATA_FILES)} chi data files")
print(f"  - {len(DATAFRAME_INPUT_FILES)} dataframe input files (STATUS/SYNTHESIS/XRD)")

rule all:
    input:
        "data/processed/synthesis_data.csv",
        "data/processed/synthesis_data.pkl",
        "data/processed/oqmd_hull_data.csv",
        "data/external/oqmd_structures/.extracted",
        "results/susceptibility/susceptibility_real_part.pdf",
        "results/susceptibility/susceptibility_imaginary_part.pdf",
        "results/susceptibility/hc2_with_fits.pdf",
        "results/susceptibility/hc2_fit_parameters.csv",
        "results/dashboard/index.html"

# Build initial dataframe without disorder to extract formulas
checkpoint build_dataframe_for_formulas:
    input:
        script="scripts/build_dataframe.py",
        src=SRC_FILES,
        data=STATUS_FILES + SYNTHESIS_FILES,
        reference=REFERENCE_DATA_FILES
    output:
        csv="data/processed/synthesis_data_no_disorder.csv",
        formulas="data/processed/formulas.txt"
    run:
        import shutil
        import pandas as pd

        # Temporarily hide disorder cache if it exists
        cache_path = Path("data/processed/disorder_cache.csv")
        backup_path = Path("data/processed/disorder_cache.csv.backup")

        if cache_path.exists():
            shutil.move(str(cache_path), str(backup_path))

        # Build dataframe
        shell("poetry run python {input.script}")

        # Move output and restore cache
        shutil.move("data/processed/synthesis_data.csv", str(output.csv))

        if backup_path.exists():
            shutil.move(str(backup_path), str(cache_path))

        # Write sorted unique formula list
        df = pd.read_csv(output.csv)
        formulas = sorted(df['formula'].dropna().unique().tolist())
        Path(output.formulas).write_text('\n'.join(formulas) + '\n')

rule compute_disorder_cache:
    input:
        script="scripts/compute_disorder_probabilities.py",
        src=SRC_FILES,
        model_files=DISORDER_MODEL_FILES,
        formulas="data/processed/formulas.txt"
    output:
        cache="data/processed/disorder_cache.csv"
    shell:
        "poetry run python {input.script}"

rule validate_oqmd_database:
    output:
        validation_marker=touch("data/processed/.oqmd_validated")
    log:
        "logs/validate_oqmd_database.log"
    shell:
        """
        poetry run python scripts/validate_oqmd_database.py > {log} 2>&1
        """

rule query_oqmd_hulls:
    input:
        script="scripts/query_oqmd_hulls.py",
        src=SRC_FILES,
        validation="data/processed/.oqmd_validated",
        csv="data/processed/synthesis_data_no_disorder.csv"
    output:
        csv="data/processed/oqmd_hull_data.csv"
    log:
        "logs/query_oqmd_hulls.log"
    shell:
        """
        poetry run python {input.script} > {log} 2>&1
        """

rule extract_oqmd_structures:
    input:
        script="scripts/extract_oqmd_structures.py",
        src=SRC_FILES,
        validation="data/processed/.oqmd_validated",
        hull_data="data/processed/oqmd_hull_data.csv"
    output:
        marker=touch("data/external/oqmd_structures/.extracted")
    log:
        "logs/extract_oqmd_structures.log"
    shell:
        """
        poetry run python {input.script} > {log} 2>&1
        """

rule build_dataframe:
    input:
        script="scripts/build_dataframe.py",
        src=SRC_FILES,
        data=DATAFRAME_INPUT_FILES,
        reference=REFERENCE_DATA_FILES,
        disorder_cache="data/processed/disorder_cache.csv",
        oqmd_hulls="data/processed/oqmd_hull_data.csv"
    output:
        csv="data/processed/synthesis_data.csv",
        pkl="data/processed/synthesis_data.pkl"
    shell:
        "poetry run python {input.script}"

rule analyze_susceptibility:
    input:
        script="scripts/analyze_susceptibility.py",
        src=SRC_FILES,
        data=CHI_DATA_FILES
    output:
        real="results/susceptibility/susceptibility_real_part.pdf",
        imag="results/susceptibility/susceptibility_imaginary_part.pdf",
        hc2="results/susceptibility/hc2_with_fits.pdf",
        params="results/susceptibility/hc2_fit_parameters.csv"
    shell:
        "poetry run python {input.script}"

rule generate_dashboard:
    input:
        script="scripts/generate_dashboard.py",
        data="data/processed/synthesis_data.pkl",
        params="results/susceptibility/hc2_fit_parameters.csv",
        chi_data=CHI_DATA_FILES,
        src=SRC_FILES
    output:
        index="results/dashboard/index.html"
    shell:
        "poetry run python {input.script}"