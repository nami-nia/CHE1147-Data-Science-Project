# IQR review decision: retain the flagged Posaconazole observation

These decision files contain the current review outcome for this observation.

Recorded at: **2026-10-14T17:44:55-04:00**, timezone **America/Toronto**.
This is the recording timestamp, not an asserted time of the earlier visual inspection.

## Observation and traceability

- Input dataset: `data/processed/subset_drugs_binary_only_DrugSolubilityData.csv`.
- Input SHA-256: `dc51b62244bb2017216f4a6d0f3c9e1bb82373e32b3c9663813006078a063797`.
- `source_row_index`: **360**, the zero-based position in the 612-row binary-only input.
- Corresponding 1-based input data-record number: **361**, excluding the header.
- This index is not a row position or raw Label in the original full dataset.
- Drug: **Posaconazole (ID 43)**.
- Solvents and their attached fractions: **ethanol (ID 130): 0.945091**; **ethyl acetate (ID 145): 0.054909**.
- Stored solvent lists are `[145.0, 130.0]` and `[0.054909, 0.945091]`; alignment was verified by pairing IDs with fractions.
- Decision CSV solvent fields use ascending solvent-ID order; each fraction remains attached to its own ID.
- Temperature: **273.15 K**; stored value: **-1.3493**.
- Screening level: `supplementary_by_system`; variable: `value`; direction: `low`.
- `value_iqr_low_by_system = True`; `value_iqr_flag_by_system = True`; `value_iqr_flag_by_drug = False`.
- Supporting figure: `reports/iqr/figures/drug_43_solvents_130_145_composition.png`.

## Current review outcome

- review_status: `reviewed_retained`.
- decision: `retain`.
- review_method: visual inspection of the composition–value plot and temperature series.
- decision_maker: project analyst (user).
- original_source_verification: not independently verified during this review.

The project analyst (user) reported visually reviewing the supporting composition–value plot and temperature series.

Decision explanation, recorded verbatim:

> The observation was flagged as a low stored value by the supplementary IQR screen within the Posaconazole–ethanol–ethyl acetate system. Visual inspection showed that it appeared consistent with the surrounding composition trend at the same temperature, rather than forming an isolated deviation. The observation was retained because the statistical flag and visual review did not establish a measurement error. Its accuracy has not been independently confirmed against the original experimental source.

## Screening snapshot and remaining limitation

`reports/iqr/iqr_review_rows.csv` remains an unchanged, automatically generated screening snapshot.
Its initial `review_status = not reviewed` entry for source_row_index 360 is superseded by this separate decision record.
The IQR report and other screening artifacts also remain unchanged historical screening results.

Retention does not cancel the statistical flag: the supplementary flag remains True.
All 612 observations remain in the binary input and annotated analysis dataset; no values, flags, quartiles, fences or counts were changed.

The visual review did not establish a measurement error. It does not prove that the observation is correct or that the flag is a confirmed false positive.
Accuracy has not been independently confirmed against the original experimental source during this review.
