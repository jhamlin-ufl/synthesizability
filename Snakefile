# Snakefile for synthesizability pipeline
from pathlib import Path

# ---------------------------------------------------------------------------
# Data file lists
# ---------------------------------------------------------------------------

RAW_DATA_FILES = [str(p) for p in Path("data/raw").rglob("*") if p.is_file()]
CHI_DATA_FILES = [str(p) for p in Path("data/raw").rglob("*chiAC*.txt")]
STATUS_FILES = [str(p) for p in Path("data/raw").rglob("STATUS")]
SYNTHESIS_FILES = [str(p) for p in Path("data/raw").rglob("SYNTHESIS")]
XRD_FILES = (
    [str(p) for p in Path("data/raw").rglob("*.xy")] +
    [str(p) for p in Path("data/raw").rglob("*.txt") if "chiAC" not in p.name]
)
DATAFRAME_INPUT_FILES = STATUS_FILES + SYNTHESIS_FILES + XRD_FILES

REFERENCE_DATA_FILES = [
    "data/external/reference/element_prices.csv",
    "data/external/reference/element_vapor_pressures.csv",
]
DISORDER_MODEL_FILES = [
    "data/external/disorder_model/hyperopt_best_model.pt",
    "data/external/disorder_model/hyperopt_config.json",
]
SUPERCON_DATA_FILES = [
    "data/external/supercon/primary.tsv",
]

# ---------------------------------------------------------------------------
# Source file lists — targeted per rule to avoid spurious reruns
# ---------------------------------------------------------------------------

SRC_IO = (
    [str(p) for p in Path("src/synthesizability/io").rglob("*.py")] +
    [str(p) for p in Path("src/synthesizability/parsers").rglob("*.py")] +
    ["src/synthesizability/formula.py"]
)

SRC_DISORDER = (
    [str(p) for p in Path("src/synthesizability/disorder_core").rglob("*.py")] +
    ["src/synthesizability/disorder.py"]
)

SRC_SUSCEPTIBILITY = [
    "src/synthesizability/susceptibility.py",
]

SRC_OQMD = [
    "src/synthesizability/oqmd.py",
]

SRC_DASHBOARD = (
    [str(p) for p in Path("src/synthesizability/dashboard_plugins").rglob("*.py")]
)

print(f"Tracking {len(RAW_DATA_FILES)} raw data files")
print(f"  - {len(CHI_DATA_FILES)} chi data files")
print(f"  - {len(DATAFRAME_INPUT_FILES)} dataframe input files (STATUS/SYNTHESIS/XRD)")

# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

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


checkpoint build_dataframe_for_formulas:
    input:
        script="scripts/build_dataframe.py",
        src=SRC_IO,
        data=STATUS_FILES + SYNTHESIS_FILES,
        reference=REFERENCE_DATA_FILES,
    output:
        csv="data/processed/synthesis_data_no_disorder.csv",
        formulas="data/processed/formulas.txt",
    run:
        import shutil
        import pandas as pd

        cache_path = Path("data/processed/disorder_cache.csv")
        backup_path = Path("data/processed/disorder_cache.csv.backup")

        if cache_path.exists():
            shutil.move(str(cache_path), str(backup_path))

        shell("poetry run python {input.script}")

        shutil.move("data/processed/synthesis_data.csv", str(output.csv))

        if backup_path.exists():
            shutil.move(str(backup_path), str(cache_path))

        df = pd.read_csv(output.csv)
        formulas = sorted(df['formula'].dropna().unique().tolist())
        Path(output.formulas).write_text('\n'.join(formulas) + '\n')


rule compute_disorder_cache:
    input:
        script="scripts/compute_disorder_probabilities.py",
        src=SRC_DISORDER,
        model_files=DISORDER_MODEL_FILES,
        formulas="data/processed/formulas.txt",
    output:
        cache="data/processed/disorder_cache.csv",
    shell:
        "poetry run python {input.script}"


rule compute_supercon_cache:
    input:
        script="scripts/compute_supercon_cache.py",
        supercon=SUPERCON_DATA_FILES,
        formulas="data/processed/formulas.txt",
    output:
        marker=touch("data/processed/.supercon_cached"),
    shell:
        "poetry run python {input.script}"


rule validate_oqmd_database:
    output:
        validation_marker=touch("data/processed/.oqmd_validated"),
    log:
        "logs/validate_oqmd_database.log"
    shell:
        "poetry run python scripts/validate_oqmd_database.py > {log} 2>&1"


rule query_oqmd_hulls:
    input:
        script="scripts/query_oqmd_hulls.py",
        src=SRC_OQMD,
        validation="data/processed/.oqmd_validated",
        csv="data/processed/synthesis_data_no_disorder.csv",
    output:
        csv="data/processed/oqmd_hull_data.csv",
    log:
        "logs/query_oqmd_hulls.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule extract_oqmd_structures:
    input:
        script="scripts/extract_oqmd_structures.py",
        src=SRC_OQMD,
        validation="data/processed/.oqmd_validated",
        hull_data="data/processed/oqmd_hull_data.csv",
    output:
        marker=touch("data/external/oqmd_structures/.extracted"),
    log:
        "logs/extract_oqmd_structures.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule build_dataframe:
    input:
        script="scripts/build_dataframe.py",
        src=SRC_IO,
        data=DATAFRAME_INPUT_FILES,
        reference=REFERENCE_DATA_FILES,
        disorder_cache="data/processed/disorder_cache.csv",
        oqmd_hulls="data/processed/oqmd_hull_data.csv",
    output:
        csv="data/processed/synthesis_data.csv",
        pkl="data/processed/synthesis_data.pkl",
    shell:
        "poetry run python {input.script}"


rule analyze_susceptibility:
    input:
        script="scripts/analyze_susceptibility.py",
        src=SRC_IO + SRC_SUSCEPTIBILITY,
        data=CHI_DATA_FILES,
    output:
        real="results/susceptibility/susceptibility_real_part.pdf",
        imag="results/susceptibility/susceptibility_imaginary_part.pdf",
        hc2="results/susceptibility/hc2_with_fits.pdf",
        params="results/susceptibility/hc2_fit_parameters.csv",
    shell:
        "poetry run python {input.script}"


rule generate_dashboard:
    input:
        script="scripts/generate_dashboard.py",
        src=SRC_DASHBOARD + SRC_IO + SRC_SUSCEPTIBILITY + SRC_OQMD,
        data="data/processed/synthesis_data.pkl",
        params="results/susceptibility/hc2_fit_parameters.csv",
        chi_data=CHI_DATA_FILES,
        supercon_cache="data/processed/.supercon_cached",
    output:
        index="results/dashboard/index.html",
    shell:
        "poetry run python {input.script}"