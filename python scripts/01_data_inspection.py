"""Step 1: Inspect the CheMixHub drug-solubility CSV files without changing them.

Run from the repository root:
    python scripts/01_data_inspection.py

Writes column-level reports into reports/data_cleaning/.
"""

from pathlib import Path
import pandas as pd

# This script should be stored in the 'scripts' folder in your Git repository.
REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT.parent / "chemixhub" / "datasets" / "drug-solubility" / "processed_data"
REPORT_DIR = REPO_ROOT / "reports" / "data_cleaning"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "solubility": DATA_DIR / "processed_DrugSolubilityData.csv",
    "compounds": DATA_DIR / "compounds.csv",
}

for name, filepath in FILES.items():
    if not filepath.is_file():
        raise FileNotFoundError(
            f"Could not find {filepath}\n"
            "Put this script in your repository's 'scripts' folder "
            "and check that the CSV file paths match your project."
        )

    df = pd.read_csv(filepath)

    print(f"\n{'=' * 60}")
    print(f"DATASET: {name} ({filepath.name})")
    print(f"Rows: {len(df):,} | Columns: {len(df.columns)}")
    print(f"Exactly duplicated rows: {df.duplicated().sum():,}")
    print("\nFirst 5 rows:")
    print(df.head().to_string(index=False))

    # Each row of this report describes one column in the CSV.
    column_report = pd.DataFrame({
        "column": df.columns,
        "data_type": df.dtypes.astype(str).values,
        "missing_count": df.isna().sum().values,
        "missing_percent": (100 * df.isna().mean()).round(2).values,
        "unique_values": df.nunique(dropna=True).values,
    })

    print("\nColumn summary:")
    print(column_report.to_string(index=False))
    output = REPORT_DIR / f"01_{name}_column_report.csv"
    column_report.to_csv(output, index=False)
    print(f"\nSaved: {output}")

    # The processed solubility file should contain these numeric variables.
    if name == "solubility":
        for column in ["value", "Temperature, K"]:
            if column not in df.columns:
                print(f"WARNING: Expected column not found: {column}")
                continue
            numeric = pd.to_numeric(df[column], errors="coerce")
            invalid_count = (df[column].notna() & numeric.isna()).sum()
            print(
                f"{column}: minimum={numeric.min()}, maximum={numeric.max()}, "
                f"non-numeric entries={invalid_count}"
            )

print("\nFinished. No source CSV files were modified.")
print(f"Review the reports in: {REPORT_DIR}")
