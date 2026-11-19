# CHE1147: Drug Solubility Prediction

This project cleans and analyzes drug-solubility measurements from CheMixHub and compares regression models for predicting reported log solubility in pure solvents and binary solvent mixtures. Inputs combine temperature, solvent mole fractions, and molecular descriptors calculated from drug and solvent SMILES using RDKit.

The main workflow is [drug_solubility_workflow.ipynb](data_cleaning_CL/drug_solubility_workflow.ipynb). It includes preprocessing, exploratory analysis, feature generation, model training and tuning, evaluation, and figure/table generation. It is the single notebook needed to run the project.

## Repository layout

```text
CHE1147-chenyue/
├── README.md
├── requirements.txt
├── .gitattributes
└── data_cleaning_CL/
    ├── drug_solubility_workflow.ipynb   # Main end-to-end workflow
    └── results/                         # Cleaned data, reports, figures, models
```

The external input CSVs are obtained separately from CheMixHub. A local virtual environment is not part of the repository.

## Get the project and input data

You need Git, Python, and a notebook frontend such as VS Code with its Python and Jupyter extensions. Install Git LFS if you want to download the saved model caches.

Run the following from a directory where both repositories can be placed side by side. The project clone uses the `chenyue-work` branch. With Git LFS installed, cloning can download approximately 1.59 GB of model caches.

```powershell
git lfs install
git clone --branch chenyue-work https://github.com/nami-nia/CHE1147-Data-Science-Project.git CHE1147-chenyue

git clone --filter=blob:none --no-checkout https://github.com/chemcognition-lab/chemixhub.git chemixhub
git -C chemixhub sparse-checkout init --cone
git -C chemixhub sparse-checkout set datasets/drug-solubility
git -C chemixhub checkout 356c16b99e38b56a91b5219d33db410fc372d538
```

The sparse checkout retrieves the drug-solubility dataset without checking out the other dataset directories. The pinned CheMixHub revision is the revision inspected for this project. If you already have a CheMixHub checkout with other work, use a separate copy to obtain this revision.

The notebook expects this arrangement; the project directory itself can have a different name, but the sibling data checkout must be named `chemixhub`:

```text
parent-directory/
├── CHE1147-chenyue/
│   └── data_cleaning_CL/
│       └── drug_solubility_workflow.ipynb
└── chemixhub/
    └── datasets/
        └── drug-solubility/
            ├── raw_data/
            │   └── DrugSolubilityData.csv
            └── processed_data/
                ├── processed_DrugSolubilityData.csv
                └── compounds.csv
```

| Required input | Purpose |
| --- | --- |
| `raw_data/DrugSolubilityData.csv` | Original measurements, temperatures, solvent SMILES, and composition percentages used to verify alignment and repair solvent representations. |
| `processed_data/processed_DrugSolubilityData.csv` | CheMixHub-formatted observations, including the regression target `value`, temperature, compound IDs, and mole fractions. |
| `processed_data/compounds.csv` | Mapping from compound IDs to SMILES for repair and molecular descriptor calculation. |

All three files are included in the pinned upstream revision. No upstream preprocessing script or CheMixHub Python-package installation is needed: the final notebook reads these CSVs directly and performs this project's additional cleaning. Keep the raw and processed rows in their original order because repair uses their correspondence. Cell 1 lists any missing inputs.

Dataset provenance is recorded in CheMixHub's [dataset documentation](https://github.com/chemcognition-lab/chemixhub/blob/356c16b99e38b56a91b5219d33db410fc372d538/datasets/README.md) and [drug-solubility metadata](https://github.com/chemcognition-lab/chemixhub/blob/356c16b99e38b56a91b5219d33db410fc372d538/datasets/drug-solubility/processed_data/croissant.json). The metadata attributes the data to Bao et al. (2024), *Towards the prediction of drug solubility in binary solvent mixtures at various temperatures using machine learning*, and specifies CC BY 4.0 for the dataset. See the upstream metadata for attribution and reuse details.

## Recreate the Python environment

The versions in [requirements.txt](requirements.txt) were captured from the project-level environment using **Python 3.14.7 on Windows**. They include the scientific libraries, RDKit, scikit-learn, XGBoost, and an IPython kernel. A fresh installation and execution on other platforms have not been verified.

From the project repository root, create a dedicated environment. Confirm that `python --version` reports the intended Python version first.

```powershell
cd CHE1147-chenyue
python --version
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip check
```

On macOS/Linux, the environment interpreter is `.venv/bin/python` instead of `.venv\Scripts\python.exe`. The requirements provide the notebook kernel; VS Code or another Jupyter frontend supplies the notebook interface.

In VS Code, open the project folder and the final notebook, then select the newly created `.venv` Python interpreter as its kernel. If you use an existing Jupyter frontend, register the kernel with:

```powershell
.\.venv\Scripts\python.exe -m ipykernel install --user --name che1147 --display-name "CHE1147"
```

Then select **CHE1147** in that frontend. Installing packages in a different environment does not make them available to the selected notebook kernel.

## Run preprocessing, training, and analysis

Open `data_cleaning_CL/drug_solubility_workflow.ipynb`, restart its kernel, and run all cells from top to bottom. Use the repository root or `data_cleaning_CL/` as the notebook working directory. The notebook discovers its folder, checks the external inputs, creates `results/` when needed, and overwrites generated CSVs and PNGs there.

Cell numbers below refer to positions in the notebook, not execution counters.

| Cells | Workflow stage |
| --- | --- |
| 1–9 | Configure paths, read the processed input, and inspect the original data and distributions. |
| 10–17 | Compare raw/processed records, repair solvent slots and fractions, and remove duplicate observations using active compositions. |
| 18–20 | Flag outliers within drug–solvent groups, apply the retention policy, and export cleaned data and the policy report. |
| 21–24 | Generate cleaned-data plots and Pearson/Spearman correlation heatmaps for binary mixtures. |
| 25–27 | Split cleaned observations and construct 27 numeric features. |
| 28–34 | Train linear regression, tune Ridge and LASSO using five-fold cross-validation, and compare their test metrics. |
| 35–36 | Load saved XGBoost/Random Forest results or run their hyperparameter searches when caches are absent. |
| 37–39 | Evaluate the tree models on the test set and regenerate the five-model comparison table and two-panel MAE/RMSE bar chart. |

For the pinned inputs, the saved cleaning run starts with **27,166 observations**, repairs **453 list-length mismatches**, expands **2,490 single-solvent representations**, and removes **9 duplicates**, leaving **27,157 observations**. The run flags **311 observations** using IQR rules and retains all of them. An IQR flag alone does not establish an erroneous measurement; rows are removed only when explicitly entered in `confirmed_outlier_errors` with supporting evidence. That dictionary is currently empty.

The target is the existing reported log-solubility `value`; the notebook does not apply an additional logarithm. Features include temperature, two solvent fractions, and eight RDKit descriptors for each of the drug and two solvent slots. An absent `"N/A"` solvent slot receives zero descriptors. The common preprocessor scales temperature and mole fractions; Ridge and LASSO additionally scale the resulting feature matrix in their own pipelines.

The split uses 80% training and 20% test observations with seed `1`, stratified using an indicator for target values above the dataset's 90th percentile. With the pinned inputs, this yields **21,725 training** and **5,432 test** observations. Hyperparameter selection uses the training set; MAE and RMSE are reported in the target's log-solubility scale. This is an observation-level split, so the test results do not establish generalization to entirely unseen drugs or solvent systems. The notebook creates its own split rather than using upstream split files or raw-data training labels.

### Use the saved model caches

The two `.joblib` files are managed by Git LFS. From the project root, obtain their actual contents if your checkout contains only LFS pointers:

```powershell
git lfs install
git lfs pull --include="data_cleaning_CL/results/*_cleaned.joblib"
```

`random_forest_tuning_results_cleaned.joblib` is approximately 1.57 GB; `xgboost_tuning_results_cleaned.joblib` is approximately 19 MB. Allow memory for loading the models in addition to their file sizes. Run the whole notebook even when using caches: preprocessing, feature generation, linear-model fitting, evaluation, and plotting still run. Each cache contains the fitted pipeline, selected hyperparameters, and training/CV metrics; test metrics are calculated in cell 37.

The current code checks cache existence and does not validate whether it matches changed data, features, or split settings. Use the captured package versions for loading saved models, and retrain after changing those settings.

### Retrain XGBoost and Random Forest

Before running the notebook, move either or both cache files outside `results/`, or rename them with a `.bak` extension:

- `data_cleaning_CL/results/xgboost_tuning_results_cleaned.joblib`
- `data_cleaning_CL/results/random_forest_tuning_results_cleaned.joblib`

Restart the kernel and run all cells. A missing cache triggers its search and saves a replacement under the original filename. If your checkout has LFS pointer files instead of downloaded models, move those pointers aside as well; `joblib.load()` cannot read them as models.

XGBoost evaluates 200 randomized parameter combinations with five-fold CV; Random Forest evaluates 100 combinations with ten-fold CV. Each search therefore performs 1,000 CV fits plus the final refit. Both use seed `42`. Full tuning can take substantial time; no fixed runtime is assumed. Ridge and LASSO are fitted and tuned on every full run.

## Regenerate and inspect results

Running the final notebook from top to bottom regenerates these outputs in `data_cleaning_CL/results/`:

| Output | Generator | Contents |
| --- | --- | --- |
| `log_solubility_distribution.png` | Cell 8 | Distribution of reported log solubility before cleaning. |
| `temperature_distribution.png` | Cell 9 | Per-drug temperature boxplots before cleaning, arranged in three panels. |
| `cleaned_drug_solubility.csv` | Cell 20 | Repaired, deduplicated observations after the outlier policy. |
| `outlier_policy_report.csv` | Cell 20 | Flagged observations, original row indices, flags, and retention/removal decisions. Default retained rows currently have blank reason/evidence fields. |
| `temperature_vs_log_solubility.png` | Cell 21 | Cleaned-data temperature versus reported log solubility. |
| `top10_drugs_log_solubility.png` | Cell 22 | Solubility observations for the ten most common drugs. |
| `solvent_fraction_vs_log_solubility.png` | Cell 23 | Binary-mixture solubility versus the first listed solvent's mole fraction. |
| `binary_mixture_correlations.png` | Cell 24 | Pearson/Spearman heatmaps across all binary mixtures. |
| `xgboost_tuning_results_cleaned.joblib` | Cell 35, when cache is absent | Selected fitted XGBoost pipeline and search/training summary. |
| `random_forest_tuning_results_cleaned.joblib` | Cell 36, when cache is absent | Selected fitted Random Forest pipeline and search/training summary. |
| `test_model_comparison.png` | Cell 38 | Five-model test MAE/RMSE comparison table. |
| `test_mae_rmse_comparison.png` | Cell 39 | Two-panel test MAE/RMSE bar chart saved at 300 dpi, with percentage reductions relative to Linear Regression and XGBoost highlighted as the best-performing model. |

The bar chart uses Linear Regression as the baseline and displays each model's error value and percentage reduction: `(Baseline Error − Model Error) / Baseline Error × 100`. Negative percentages indicate higher error than the baseline. Both panels use the same vertical scale.

Other diagnostic tables, search summaries, and metrics are displayed in notebook outputs. Save the executed notebook to retain them.

Use regenerated outputs to report results after any code or input change. Exact strip-plot pixel positions may vary because its random jitter is generated before the notebook sets the global NumPy seed. An earlier verification of the 38-cell notebook completed in a fresh Python 3.14.7 environment on Windows installed from `requirements.txt`: `pip check` passed and all 38 cells executed successfully, saving regenerated results. That run loaded the existing XGBoost and Random Forest caches; their full hyperparameter searches were not repeated. The current notebook contains 39 cells; this earlier verification does not establish a successful fresh run of the current version.

## Troubleshooting

- **Missing CSVs:** check the sibling `chemixhub/` directory name, the pinned revision, and the three paths listed in cell 1. Do not substitute the final cleaned output for the original processed input.
- **Import errors:** confirm that the notebook kernel is the environment where `requirements.txt` was installed.
- **Model-loading errors:** check that Git LFS downloaded the actual model files and that the installed package versions match the requirements. To regenerate a cache, move it aside and rerun the notebook as described above.
- **Undefined variables after running a later cell:** restart the kernel and execute all preceding cells in order.
