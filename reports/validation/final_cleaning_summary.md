# Final validation and data-cleaning summary

**STRUCTURAL / PRESERVATION CHECKS PASS; TARGET DEFINITION ESTABLISHED; INHERITED UNIT METADATA CORRECTION PENDING**

This is final checking and documentation of the existing strictly binary analysis population. No filtering, correction, outlier removal, screening rerun or new cleaned CSV was performed.

Documentation amendment: sources reviewed on **2026-10-14 (America/Toronto)** established the target definition and updated the four documentation questions below. That documentation-only amendment did not rerun structural validation. The original **37 PASS, 0 FAIL, 4 UNRESOLVED** record remains unchanged in the validation-check CSV. Numerical values, inherited CSV labels, screening results and saved review decisions are unchanged.

## Analysis population

- Analysis dataset: `data/processed/subset_drugs_binary_only_DrugSolubilityData.csv`.
- Input SHA-256: `dc51b62244bb2017216f4a6d0f3c9e1bb82373e32b3c9663813006078a063797`.
- Companion annotation: `data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv`.
- Current review outcomes: `reports/iqr/iqr_review_decisions.csv` and `.md`.

| Drug ID | Drug | Observed records |
|---|---|---:|
| 5 | Itraconazole | 360 |
| 43 | Posaconazole | 162 |
| 86 | Ketoconazole | 90 |
| Total | | 612 |

## Integrity and duplicates

Empty required cells by original column: {"Temperature, K": 0, "cmp_ids_drug": 0, "cmp_ids_solvent": 0, "cmp_mole_fractions_drug": 0, "cmp_mole_fractions_solvent": 0, "property": 0, "unit": 0, "value": 0}.

Structural checks cover finite stored values and positive finite kelvin temperatures, exactly one target drug, exactly two distinct solvent IDs with aligned strictly positive finite fractions <=1, absolute fraction-sum tolerance 1e-6 (relative tolerance zero), and aligned valid drug fractions. No fraction normalization, imputation or unit change occurs. Each failure is recorded in the check CSV.

Exact duplicate groups across the original eight decoded cells: **0**.
Equivalent groups after in-memory paired solvent canonicalization: **0**; groups differing in solvent order: **0**.
Audit identifiers are excluded from all duplicate keys. All other cells, including stored value, remain literal; ID/fraction numeric comparisons use exact equality, no tolerance or rounding. Repeated conditions alone do not establish erroneous copies. No candidate is automatically removed.

## Verified research-scope selection

Original subset: `data/processed/subset_drugs_DrugSolubilityData.csv`, identified from INPUT_NAME in the filtering script.
Actual reconciliation: **684 original - 72 pure-solvent scope exclusions = 612 binary observations**.

| Drug ID | Original | Pure excluded | Binary retained | Excluded original-subset data records |
|---|---:|---:|---:|---|
| 5 | 390 | 30 | 360 | 622–651 |
| 43 | 189 | 27 | 162 | 442–450, 652–669 |
| 86 | 105 | 15 | 90 | 670–684 |

Those locations are **1-based data records in the original 684-row subset**, excluding its header; they are calculated from actual compositions, not used as selection instructions. Every exclusion has exactly one positive solvent fraction, including 9 pure ethyl-acetate Posaconazole records. All remaining records have two active solvents, and their complete original cells match the binary input in order. These are scope exclusions, not outlier removals.

Excluded counts by drug ID:positive solvent ID: {"43:130": 9, "43:140": 9, "43:145": 9, "5:127": 5, "5:130": 5, "5:131": 5, "5:134": 5, "5:135": 5, "5:147": 5, "86:129": 5, "86:141": 5, "86:144": 5}.

Filtering evidence: `python scripts/06_remove_pure_solvents.py:45`; `python scripts/06_remove_pure_solvents.py:137`; `python scripts/06_remove_pure_solvents.py:158`; `python scripts/06_remove_pure_solvents.py:164`.

## Existing IQR results and completed review

The existing procedure uses median-of-halves quartiles, excluding the overall middle observation for odd n, and strict 1.5 x IQR fences. Primary comparisons are per drug; supplementary comparisons are per drug and sorted solvent pair for n >= 10, a project convention, not a universal statistical minimum. Zero-IQR groups use the unchanged strict rule. See `reports/iqr/iqr_report.md` (Method section) and `python scripts/07_identify_iqr_outliers.py` (`quartiles` and `make_screen` sections). The screening script was not executed or changed.

Existing unique flag counts: {"temperature_by_drug": 0, "temperature_by_system": 0, "value_by_drug": 0, "value_by_system": 1}.
Review table unique observations: **1**; source_row_index values: [360].
Summary group-variable rows reconciled: **30**. Stored group counts, percentages, flag columns, ranges and fence algebra agree. No new quartiles or flags were calculated.

The separate decision is authoritative: **reviewed_retained / retain**, made by the project analyst (user) after visual inspection of the composition–value plot and temperature series.

Observation: **source_row_index 360** (zero-based position in the 612-row binary input, not in the original full dataset), Posaconazole ID 43, ethanol ID 130 fraction 0.945091, ethyl acetate ID 145 fraction 0.054909, 273.15 K, stored value -1.3493.

The recorded input SHA-256 matches the current binary file. The supplementary low-value flag remains **True**, and the primary flag remains **False**. The unchanged automatically generated review snapshot's initial 'not reviewed' entry is superseded by the separate decision record. Retention does not cancel a statistical flag.

> The observation was flagged as a low stored value by the supplementary IQR screen within the Posaconazole–ethanol–ethyl acetate system. Visual inspection showed that it appeared consistent with the surrounding composition trend at the same temperature, rather than forming an isolated deviation. The observation was retained because the statistical flag and visual review did not establish a measurement error. Its accuracy has not been independently confirmed against the original experimental source.

Its accuracy was **not independently verified against the original experimental source**. It is neither proven correct nor a confirmed error or false-positive flag.

Frost's “Guidelines for Dealing with Outliers” supports retaining observations unless exclusion has a specific justification; statistical extremeness alone is insufficient. [4] This is external methodological guidance, not evidence of this measurement's accuracy. The project-specific assessment comes from the user's saved visual review of the [composition–value plot](../iqr/figures/drug_43_solvents_130_145_composition.png) and the [separate retention decision](../iqr/iqr_review_decisions.md).

**Zero observations were removed during IQR screening/review**, confirmed by complete original-cell/order preservation.

## Subset-specific earlier history

Content-verified mapping covers all **684** selected observations, using `reports/data_cleaning/04_duplicate_removal/04_duplicate_row_mapping.csv`, the cumulative dataset selected by script 05, saved ID-repaired and upstream processed tables, raw Labels, raw drug SMILES, temperatures and target values. Labels were read, not inferred from positions.

- Logged solvent-ID-only corrections in the original 684-row subset: **9**.
- Such corrections in the retained 612-row binary population: **0**.
- Earlier exact-duplicate removals involving these three drug IDs: **0**.
- Correction membership and scope are established by the record details below; values, temperatures and fractions were not repaired.
- These counts are derived from the selected records, repair change log and duplicate-removal mapping/log. Full-dataset totals are not assigned to the three-drug subset. Evidence: `reports/data_cleaning/03_repair/03_id_mismatch_change_log.csv` and `reports/data_cleaning/04_duplicate_removal/04_duplicate_removal_log.csv`.
- Raw nonzero solvent ID/fraction correspondence uses the existing 1e-12 serialization tolerance only for provenance matching, not duplicate detection or measurement similarity.

Historical correction details: [{"active_solvents": 1, "drug_id": 43, "original_data_row": 9418, "positive_solvent_ids": [145], "raw_label": "9417", "subset_data_row": 442}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9419, "positive_solvent_ids": [145], "raw_label": "9418", "subset_data_row": 443}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9420, "positive_solvent_ids": [145], "raw_label": "9419", "subset_data_row": 444}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9421, "positive_solvent_ids": [145], "raw_label": "9420", "subset_data_row": 445}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9422, "positive_solvent_ids": [145], "raw_label": "9421", "subset_data_row": 446}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9423, "positive_solvent_ids": [145], "raw_label": "9422", "subset_data_row": 447}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9424, "positive_solvent_ids": [145], "raw_label": "9423", "subset_data_row": 448}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9425, "positive_solvent_ids": [145], "raw_label": "9424", "subset_data_row": 449}, {"active_solvents": 1, "drug_id": 43, "original_data_row": 9426, "positive_solvent_ids": [145], "raw_label": "9425", "subset_data_row": 450}].

## Target provenance versus measurement definitions

### Four documentation questions: current status after source review

| Question | Current status |
|---|---|
| 1. What logarithm base was used? | Established from the original authors' published code. |
| 2. What concentration basis was used? | Original target definition established; inconsistent inherited g/L metadata documented; metadata correction pending team discussion. |
| 3. What transformations occurred before value was created? | Main processing documented; earlier record-by-record literature extraction and preparation not independently audited. |
| 4. Why was the flagged observation retained? | Retention decision complete; independent verification against the original experimental publication not performed. |

### 1. Logarithm base

The original authors' external notebook, “Unit conversion,” `calculate_logs()`, explicitly applies `np.log10(df_updated['Solubility (g/100g)'])`. This establishes base 10 from the published code, rather than from the paper alone. [2]

### 2. Concentration basis and inherited unit discrepancy

Bao et al., **Results → Dataset overview** and **Methods → Data collection**, report LogS with S in g/100 g and describe standardizing literature measurements. [1] In `calculate_logs()`, drug mass is divided by the combined mass of drug and both solvents, then multiplied by 100. The denominator is total solution mass. [2]

The documented target is **LogS: the base-10 logarithm of the numerical solubility expressed in g/100 g of solution.** Thus `value` contains a logarithmic target, not an ordinary concentration measurement. [1, 2]

The provided Chemixhub code's `name_processing()` renames `LogS` to `value`, assigns `property = "Log solubility"`, and assigns `unit = "g/L"`. No corresponding numerical target conversion to g/L is implemented. Its division by 100 concerns composition percentages, not the target. [3]

**Metadata correction is pending discussion at the next group meeting.** This documentation update preserves all numerical values and inherited CSV labels, including `g/L`; it performs no conversion or metadata correction. The assigned label does not establish a g/L basis for the stored logarithmic target. [2, 3] The supplied `../chemixhub/datasets/drug-solubility/processed_data/croissant.json:51` also gives the generic example `g/100g solvent`; that wording differs from the total-solution denominator in the authors' code and is not used to override it. The dataset READMEs describe fields generically and do not resolve this discrepancy.

### 3. Documented transformation chain and provenance limit

The paper describes collection and standardization of literature solubility measurements. [1] The notebook's “Unit conversion” section implements solvent mole/mass-fraction interconversion in `calculate_fractions()`, followed by this target calculation in `calculate_logs()`: **standardized solubility input in mol/mol → component masses using composition and molecular weights → drug mass fraction → g/100 g of solution → base-10 logarithm, LogS**. [2] Chemixhub subsequently renames `LogS` to `value`. [3]

The raw header `../chemixhub/datasets/drug-solubility/raw_data/DrugSolubilityData.csv:1` names `LogS`. Existing content-verified project mappings establish numerical equality of raw `LogS` and stored `value` for all 684 selected records. This is internal preservation evidence, not independent verification of every source experiment or historical conversion.

The main calculation is documented; earlier record-by-record literature extraction and preparation of the standardized input have **not** been independently audited. Original experimental tables/SI and record-level preparation records would be needed for that audit. This is a provenance limitation, not evidence that the measurements are wrong.

### 4. Retention of source_row_index 360

The saved decision is **retain** for Posaconazole ID 43 in ethanol–ethyl acetate, 273.15 K, stored value -1.3493. `source_row_index = 360` is the zero-based position in the **612-row binary-only input**, not the original full dataset. The supplementary per-system low-value flag remains **True**; the primary per-drug flag is **False**.

The user reported that the point appeared consistent with the surrounding composition trend at the same temperature. The [decision CSV](../iqr/iqr_review_decisions.csv), [decision explanation](../iqr/iqr_review_decisions.md), and [supporting plot](../iqr/figures/drug_43_solvents_130_145_composition.png) are the project evidence. Frost supports the treatment principle that extremeness alone is not sufficient justification for removal. [4] Neither that article nor this visual decision independently verifies the original measurement. Its accuracy is unconfirmed, and the flag has not been cancelled or established as a false positive.

Actual distinct stored metadata: {"cmp_mole_fractions_drug": ["[1.0]"], "property": ["Log solubility"], "unit": ["g/L"]}. Constant metadata and valid drug fractions establish consistency, not physical definitions.

## Historical validation outcomes and safeguards

Original saved result: **37 PASS, 0 FAIL, 4 UNRESOLVED**. The unchanged [validation-check CSV](final_validation_checks.csv) preserves the original check definitions, statuses and evidence; the later source review does not relabel them as a new all-PASS result.

Protected inputs/companions at the original validation start: 61 files, including the validator itself. The original before/after hashes excluded only that validator's two output reports; protected contents were checked before and after publication.

Current validation status totals: {"PASS": 37, "UNRESOLVED": 4}. The same validation checks are rerun; their historical documentation questions are retained separately from the approved source-review conclusions. A rerun does not undo those later findings.

This run protects 61 pre-existing inputs/companions, including the validator itself. Their hashes are checked before and after report publication; existing datasets, screening artifacts, decisions, scripts and figures are not modified by validation.

Failed checks in this run:

- None.

The four original documentation questions were `retained_observation_original_accuracy`, `target_logarithm_base`, `target_concentration_basis`, and `target_precompilation_conversion`. Their original `UNRESOLVED` rows remain intact. The current status table above records the later source-review findings: the base and original basis are established, the inherited unit label remains inconsistent pending team discussion, the main transformation is documented but earlier record-level preparation is unaudited, and retention is complete without original-source accuracy verification.

Earlier IQR reports, plots and machine-readable validation records remain unchanged historical artifacts. Their original statements that definitions were unresolved are superseded **for this documentation review** by the sourced target definition here; their counts, flags and review history are preserved. Chemical-structure validation is outside this reporting task.

## References and supporting sources

Sources were inspected on **2026-10-14 (America/Toronto)** during the recorded source review. Normal validation reruns do not access or reverify external sources. The published article was used; no local `1147_paper.pdf` was located. All four required sources were accessible. No source files were downloaded to the project.

**[1] Original research publication.** Bao, Z.; Tom, G.; Cheng, A.; Watchorn, J.; Aspuru-Guzik, A.; Allen, C. (2024). [Towards the prediction of drug solubility in binary solvent mixtures at various temperatures using machine learning](https://doi.org/10.1186/s13321-024-00911-3). *Journal of Cheminformatics*, 16, Article 117. [Published full text](https://link.springer.com/article/10.1186/s13321-024-00911-3). Relevant locations: **Results → Dataset overview**; **Methods → Data collection**; and the standalone **Data availability** heading, which links the authors' [Solubility_ML repository](https://github.com/Christine-Allen-Lab/Solubility_ML). Data availability is not a subsection of Supplementary Information. The explicit logarithm base is confirmed by code [2].

**[2] External supporting source code linked by the original publication.** Christine-Allen-Lab (n.d.). [Dataset_processing_overview.ipynb, Solubility_ML](https://github.com/Christine-Allen-Lab/Solubility_ML/blob/main/Codes_and_datasets/1_Codes_and_datasets_original/Exp_1_Raw_Dataset/Dataset_processing_overview.ipynb). Inspected [version e1b5c16712d14a398bf7229d37fb7775f31ea166](https://github.com/Christine-Allen-Lab/Solubility_ML/blob/e1b5c16712d14a398bf7229d37fb7775f31ea166/Codes_and_datasets/1_Codes_and_datasets_original/Exp_1_Raw_Dataset/Dataset_processing_overview.ipynb), read as notebook JSON in memory. Relevant location: **Unit conversion**, `calculate_fractions()` and `calculate_logs()` (notebook cells 38 and 39, zero-based; section heading in cell 37). This is externally accessed author code, not a notebook originally supplied in the Chemixhub dataset.

**[3] Provided dataset processing code / provenance evidence.** Chemixhub (n.d.). `data_processing.py`, supplied drug-solubility dataset. Dataset-repository-relative path: `datasets/drug-solubility/processed_data/data_processing.py`. Actual local location, relative to our project root: `../chemixhub/datasets/drug-solubility/processed_data/data_processing.py`, established from `08_final_validation.py`'s `DEFAULT_SOURCE` and recorded provenance. [Original repository URL](https://github.com/chemcognition-lab/chemixhub/blob/main/datasets/drug-solubility/processed_data/data_processing.py); [supplied checkout version 356c16b99e38b56a91b5219d33db410fc372d538](https://github.com/chemcognition-lab/chemixhub/blob/356c16b99e38b56a91b5219d33db410fc372d538/datasets/drug-solubility/processed_data/data_processing.py). The inspected local copy matches that checkout's tracked source after CRLF/LF normalization; its actual byte SHA-256 is `938f1a604477dcde28ddcc720bbe1b0cfdd03f619d28b182a30bf4d231aa7b17`. Relevant function: `name_processing()`, including composition-fraction handling, target renaming and label assignments. The local provided copy is the evidence for our processing history. The current online `main` version was not inspected or assumed identical. This is neither independent literature evidence nor code written by our project team.

**[4] External methodological reference for outlier treatment.** Frost, J. (n.d.). [Guidelines for Removing and Handling Outliers in Data](https://statisticsbyjim.com/basics/remove-outliers/). *Statistics By Jim*. Relevant location: **Guidelines for Dealing with Outliers**. Used only for the exclusion/retention principle; the project plot and saved review record supply the observation-specific assessment. No article publication date was established.

## Internal project evidence

These project-generated records support our counts, preservation checks and visual-review decision; external publications are not cited as evidence for those project results.

- Cleaning and selection: [03_repair_solvent_pairs.py](../../python%20scripts/03_repair_solvent_pairs.py), [04_remove_duplicates.py](../../python%20scripts/04_remove_duplicates.py), [05_subset_drugs_selection.py](../../python%20scripts/05_subset_drugs_selection.py), and [06_remove_pure_solvents.py](../../python%20scripts/06_remove_pure_solvents.py).
- Earlier preservation evidence: [ID-repair change log](../data_cleaning/03_repair/03_id_mismatch_change_log.csv), [ID-repair validation summary](../data_cleaning/03_repair/03_id_mismatch_repair_validation_summary.txt), [complete duplicate-removal row mapping](../data_cleaning/04_duplicate_removal/04_duplicate_row_mapping.csv), and [duplicate-removal log](../data_cleaning/04_duplicate_removal/04_duplicate_removal_log.csv).
- Screening implementation and historical results: [07_identify_iqr_outliers.py](../../python%20scripts/07_identify_iqr_outliers.py), [IQR summary](../iqr/iqr_summary.csv), [IQR report](../iqr/iqr_report.md), [unchanged screening review snapshot](../iqr/iqr_review_rows.csv), and [annotated dataset](../../data/processed/subset_drugs_binary_only_DrugSolubilityData_iqr_flagged.csv).
- Observation-specific review: [composition–value plot](../iqr/figures/drug_43_solvents_130_145_composition.png), [current decision CSV](../iqr/iqr_review_decisions.csv), and [current decision Markdown](../iqr/iqr_review_decisions.md).
- Validation: [08_final_validation.py](../../python%20scripts/08_final_validation.py) and [historical check results](final_validation_checks.csv).

## Reproduction and outputs

From the repository root:

```powershell
py -B "python scripts/08_final_validation.py"
```

**Reporting synchronization (2026-10-15, America/Toronto):** `render_summary()` now independently generates the approved source findings, references, limitations and current review outcome. Validation was rerun for this synchronization task using the unchanged checks. This fixed date records the synchronization event, not the date of every later execution.

Normal reruns generate reports from the validation inputs plus source-review text maintained in this script; the existing output Markdown is never used as a report template. No live web access is required. `publish_reports()` reuses byte-identical reports and refuses any unexpected existing report difference. The approved Markdown was synchronized once after snapshot verification and candidate-diff review; no force-overwrite option was added.

The provided local Chemixhub checkout is required at `../chemixhub/datasets/drug-solubility` by default, together with its dataset/root READMEs and the existing project preservation evidence. An existing matching checkout can be selected with `--source-dir`; no source data is downloaded.

- `python scripts/08_final_validation.py`.
- `reports/validation/final_validation_checks.csv` — unchanged machine-readable statuses, expectations, observations and evidence.
- `reports/validation/final_cleaning_summary.md` — this generated document.

The existing binary-only CSV remains the analysis dataset; the annotated CSV and separate decision files are companion audit records. No newly renamed cleaned CSV was created. No files were staged, committed or pushed; the existing branch and unrelated work were preserved.
