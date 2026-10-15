# IQR screening of the strictly binary three-drug dataset

All 612 observations remain in input order. Flags identify statistical review candidates; they are not confirmed errors. No source-based treatment decisions, removals, transformations or Git operations were made.

## Input and verification

- Input: `data/processed/subset_drugs_binary_only_DrugSolubilityData.csv`; SHA-256 `dc51b62244bb2017216f4a6d0f3c9e1bb82373e32b3c9663813006078a063797`.
- Annotated output: `data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv`.
- Validated counts: ID 5 = 360, ID 43 = 162, ID 86 = 90; total = 612.
- Every record has one valid target drug, two distinct paired solvent IDs, and two finite strictly positive fractions no greater than 1. Solvent fractions sum to 1 with absolute tolerance 1e-6 and relative tolerance zero. Drug fractions are valid and aligned.
- Stored value and kelvin temperature are finite numeric values; every temperature is positive.
- Complete saved-file read-back checks passed: all original eight cells and column order, 612-row order, source indices 0–611, drug counts, summary/flag agreement, and unique review union.
- Input, original 684-record subset, other datasets and existing processing scripts were protected by before/after hashes and remained unchanged.
- `source_row_index` is the original zero-based position in this 612-row INPUT ONLY. It is not a row position or raw Label in the full dataset.

## Method

For each group and variable, sort a temporary numerical copy. Q2 is the median of all values. For even n, split into equal halves. For odd n, exclude the overall middle observation. Q1 and Q3 are the medians of the lower and upper halves. IQR = Q3 − Q1; fences are Q1 − 1.5 × IQR and Q3 + 1.5 × IQR. Flag only strict inequalities; fence equality is retained without a flag. No interpolated percentile defaults are used.

Decimal arithmetic uses precision derived from the input (at least 50 digits), with no rounding in flag comparisons. Tables in this report round only for display; the summary CSV retains the calculated quartiles and fences. Percentages use that particular group's n.

If IQR is zero, the interval has zero width and the same strict rule applies, without epsilon. Constant groups have no flags. Nonconstant zero-IQR groups can flag departures from the common value and require caution.

Primary comparisons combine all rows of each drug separately. Supplementary comparisons use drug plus an order-independent sorted solvent-ID pair. The n >= 10 requirement is a project screening convention, not a universal statistical minimum. Smaller groups are 'not screened—small group'; their supplementary flags are missing, not False. Their counts/ranges remain reported and their primary screens remain active.

IDs are categorical. Composition is an experimental condition: no solvent-fraction IQR or automatic composition-based treatment was applied. The lower solvent ID is the reference for composition plots, with its own attached fraction.

Calculation checks: 7 PASS (even, odd, constant, nonconstant zero-IQR, exact-fence, ID/fraction alignment and missing small-group flags).

## Primary results

| Drug ID | Variable | n | Distinct | Min | Q1 | Median | Q3 | Max | IQR | Lower fence | Upper fence | Low | High | Total | % | Zero IQR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 5 | value | 360 | 357 | -1.312 | -0.3574 | 0.6175 | 1.27405 | 1.5587 | 1.63145 | -2.804575 | 3.721225 | 0 | 0 | 0 | 0.00 | False |
| 5 | Temperature, K | 360 | 5 | 288.15 | 293.15 | 298.15 | 303.15 | 308.15 | 10 | 278.15 | 318.15 | 0 | 0 | 0 | 0.00 | False |
| 43 | value | 162 | 162 | -2.0608 | -1.2972 | -0.94725 | -0.6985 | -0.2508 | 0.5987 | -2.19525 | 0.19955 | 0 | 0 | 0 | 0.00 | False |
| 43 | Temperature, K | 162 | 9 | 273.15 | 283.15 | 293.15 | 303.15 | 313.15 | 20 | 253.15 | 333.15 | 0 | 0 | 0 | 0.00 | False |
| 86 | value | 90 | 90 | -2.12 | -0.7049 | 0.68895 | 1.1261 | 1.5116 | 1.831 | -3.4514 | 3.8726 | 0 | 0 | 0 | 0.00 | False |
| 86 | Temperature, K | 90 | 5 | 293.15 | 298.15 | 303.15 | 308.15 | 313.15 | 10 | 283.15 | 323.15 | 0 | 0 | 0 | 0.00 | False |

Temperature flags mean uncommon experimental temperatures within the comparison group; they are not automatically incorrect measurements or removal reasons.

| Drug ID | Distinct temperatures, K |
|---|---|
| 5 | 288.15, 293.15, 298.15, 303.15, 308.15 |
| 43 | 273.15, 278.15, 283.15, 288.15, 293.15, 298.15, 303.15, 308.15, 313.15 |
| 86 | 293.15, 298.15, 303.15, 308.15, 313.15 |

Unique primary flagged rows (overlaps counted once):

| Drug ID | Target only | Temperature only | Both | Either |
|---|---:|---:|---:|---:|
| 5 | 0 | 0 | 0 | 0 |
| 43 | 0 | 0 | 0 | 0 |
| 86 | 0 | 0 | 0 | 0 |
| All | 0 | 0 | 0 | 0 |

## Supplementary drug–solvent-system results

Primary and supplementary flags are separate columns and distinct comparisons.

| Drug ID | Solvent pair | Variable | n | Range | Q1 | Q3 | Fences | Low | High | Total | % | Status |
|---|---|---|---:|---|---:|---:|---|---:|---:|---:|---:|---|
| 5 | 127/134 | value | 45 | -0.9137 to 1.5587 | -0.20425 | 1.422 | -2.643625 to 3.861375 | 0 | 0 | 0 | 0 | screened |
| 5 | 127/134 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 5 | 127/135 | value | 45 | -0.9108 to 1.512 | -0.32425 | 1.3545 | -2.842375 to 3.872625 | 0 | 0 | 0 | 0 | screened |
| 5 | 127/135 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 5 | 130/134 | value | 45 | -1.312 to 1.4843 | -0.396 | 1.28915 | -2.923725 to 3.816875 | 0 | 0 | 0 | 0 | screened |
| 5 | 130/134 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 5 | 130/135 | value | 45 | -1.16 to 1.4696 | -0.53315 | 1.1926 | -3.121775 to 3.781225 | 0 | 0 | 0 | 0 | screened |
| 5 | 130/135 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 5 | 131/134 | value | 45 | -1.2633 to 1.4182 | -0.27935 | 1.2766 | -2.613275 to 3.610525 | 0 | 0 | 0 | 0 | screened |
| 5 | 131/134 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 5 | 131/135 | value | 45 | -1.2531 to 1.4438 | -0.5061 | 1.19845 | -3.062925 to 3.755275 | 0 | 0 | 0 | 0 | screened |
| 5 | 131/135 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 5 | 134/147 | value | 45 | -1.2933 to 1.4225 | -0.16845 | 1.27245 | -2.3298 to 3.4338 | 0 | 0 | 0 | 0 | screened |
| 5 | 134/147 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 5 | 135/147 | value | 45 | -1.2612 to 1.4371 | -0.4428 | 1.20345 | -2.912175 to 3.672825 | 0 | 0 | 0 | 0 | screened |
| 5 | 135/147 | Temperature, K | 45 | 288.15 to 308.15 | 293.15 | 303.15 | 278.15 to 318.15 | 0 | 0 | 0 | 0 | screened |
| 43 | 130/145 | value | 81 | -1.3493 to -0.2508 | -0.8659 | -0.54925 | -1.340875 to -0.074275 | 1 | 0 | 1 | 1.234568 | screened |
| 43 | 130/145 | Temperature, K | 81 | 273.15 to 313.15 | 283.15 | 303.15 | 253.15 to 333.15 | 0 | 0 | 0 | 0 | screened |
| 43 | 140/145 | value | 81 | -2.0608 to -0.7143 | -1.5027 | -1.05595 | -2.172825 to -0.385825 | 0 | 0 | 0 | 0 | screened |
| 43 | 140/145 | Temperature, K | 81 | 273.15 to 313.15 | 283.15 | 303.15 | 253.15 to 333.15 | 0 | 0 | 0 | 0 | screened |
| 86 | 129/141 | value | 45 | -2.12 to 1.5116 | -0.78865 | 1.31925 | -3.9505 to 4.4811 | 0 | 0 | 0 | 0 | screened |
| 86 | 129/141 | Temperature, K | 45 | 293.15 to 313.15 | 298.15 | 308.15 | 283.15 to 323.15 | 0 | 0 | 0 | 0 | screened |
| 86 | 129/144 | value | 45 | -1.742 to 1.2164 | -0.19925 | 1.06015 | -2.08835 to 2.94925 | 0 | 0 | 0 | 0 | screened |
| 86 | 129/144 | Temperature, K | 45 | 293.15 to 313.15 | 298.15 | 308.15 | 283.15 to 323.15 | 0 | 0 | 0 | 0 | screened |

Supplementary unique flagged rows: **1**. Unique union of ALL completed screens in the review table: **1**.

| source_row_index | Input data record (1-based) | Drug ID | Solvent pair | Stored value | Temperature, K | Triggers | Review status |
|---:|---:|---|---|---|---|---|---|
| 360 | 361 | 43 | 130/145 | -1.3493 | 273.15 | value_by_system:low | not reviewed |

## Metadata and unresolved measurement definitions

| Field | Distinct stored cells | Constant? |
|---|---|---|
| property | `["Log solubility"]` | True |
| unit | `["g/L"]` | True |
| cmp_mole_fractions_drug | `["[1.0]"]` | True |

Drug-fraction cells were safely parsed and checked against the singleton drug lists; the confirmed constant field was not subjected to IQR. Neither metadata nor categorical IDs receive IQR flags.

The local upstream `../chemixhub/datasets/drug-solubility/processed_data/data_processing.py` renames `LogS` to `value` and assigns `property = Log solubility` and `unit = g/L`. These operations and constant labels do not document the logarithm base or concentration denominator. Those definitions remain unresolved. Values were screened exactly as stored: no extra logarithm or back-transformation. Figures say **Stored value**.

Compound display-name source: `../chemixhub/datasets/drug-solubility/processed_data/compounds.csv` (available; explicit compound IDs). Names are plot labels, not chemical-structure validation.

Flags depend on comparison scope and can reflect experimental composition/temperature patterns. Inspect the composition plots and source evidence before any later treatment. No observation has been removed or declared erroneous; review_status remains not reviewed.

## Figures and reproduction

Boxplots use precomputed Q1, median and Q3 passed to Matplotlib bxp, not percentile defaults. Whiskers end at the most extreme observed values inside the exact fences; fence lines are drawn separately. Composition plots identify the lower-ID reference solvent, color all observations by temperature, mark primary value flags with red crosses and supplementary value flags with orange diamonds. Lines connect observed points at each temperature for context; they are not fitted models or treatment rules.

Dependencies used: Python 3.10.21, Matplotlib 3.10.9, NumPy 2.2.6.
Numerical analysis, validation and CSV handling otherwise use the standard library.

Run from the repository root (defaults resolve relative to the script's repository):

```powershell
uv run --no-project --with matplotlib==3.10.9 --with numpy==2.2.6 python -B "python scripts/07_identify_iqr_outliers.py"
```

Identical rerun artifacts are reused unchanged; any conflicting existing artifact is refused before publishing. Temporary staging and plotting caches are removed.

Generated outputs:

- `data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv`
- `reports/iqr/iqr_summary.csv`
- `reports/iqr/iqr_review_rows.csv`
- `reports/iqr/iqr_report.md`
- `reports/iqr/figures/drug_5_value_boxplot.png`
- `reports/iqr/figures/drug_5_temperature_boxplot.png`
- `reports/iqr/figures/drug_43_value_boxplot.png`
- `reports/iqr/figures/drug_43_temperature_boxplot.png`
- `reports/iqr/figures/drug_86_value_boxplot.png`
- `reports/iqr/figures/drug_86_temperature_boxplot.png`
- `reports/iqr/figures/drug_5_solvents_127_134_composition.png`
- `reports/iqr/figures/drug_5_solvents_127_135_composition.png`
- `reports/iqr/figures/drug_5_solvents_130_134_composition.png`
- `reports/iqr/figures/drug_5_solvents_130_135_composition.png`
- `reports/iqr/figures/drug_5_solvents_131_134_composition.png`
- `reports/iqr/figures/drug_5_solvents_131_135_composition.png`
- `reports/iqr/figures/drug_5_solvents_134_147_composition.png`
- `reports/iqr/figures/drug_5_solvents_135_147_composition.png`
- `reports/iqr/figures/drug_43_solvents_130_145_composition.png`
- `reports/iqr/figures/drug_43_solvents_140_145_composition.png`
- `reports/iqr/figures/drug_86_solvents_129_141_composition.png`
- `reports/iqr/figures/drug_86_solvents_129_144_composition.png`

Saved CSV SHA-256 values:

- `data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv`: `0a603e3b05e1243beb268440d781d0eb1ec2f1545e3a1f317c8431586f5896e5`
- `reports/iqr/iqr_summary.csv`: `88671e1ca4e82a64edea18dddd12bc5917804137fc07006e31df19737ac4c632`
- `reports/iqr/iqr_review_rows.csv`: `f5b58f3ca360d1c5a5541dbcea7d9645a14deed0feffda1e632bcb5e35f2d667`
