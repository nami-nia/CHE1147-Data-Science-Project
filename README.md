# CHE1147: Drug Solubility Prediction

This project cleans and analyzes drug-solubility measurements from CheMixHub and compares regression models for predicting reported log solubility in pure solvents and binary solvent mixtures. Inputs combine temperature, solvent mole fractions, and molecular descriptors calculated from drug and solvent SMILES using RDKit.

The main workflow is [project_notebook.ipynb](notebook_and_results/project_notebook.ipynb). It includes preprocessing, exploratory analysis, feature generation, model training and tuning, evaluation, and figure/table generation. It is the single notebook needed to run the project and contains 38 code cells.

The project repository is [nami-nia/CHE1147-Data-Science-Project](https://github.com/nami-nia/CHE1147-Data-Science-Project), on the `chenyue-work` branch. Use `project_code_results/` as the working directory for environment setup and notebook execution.

## Repository layout

```text
project_code_results/
├── README.md
├── requirements.txt
├── .gitattributes
├── .gitignore
├── .python-version
├── pyproject.toml
├── src/                        # Package scaffold; analysis runs in the notebook
├── Data/
│   └── chemixhub/               # External data checkout (Git-ignored)
└── notebook_and_results/
    ├── project_notebook.ipynb   # Main end-to-end workflow
    └── results/                # Cleaned data, reports, figures, models
```

The external input CSVs are obtained separately from CheMixHub. The `Data/chemixhub/` checkout, virtual environments such as `.venv/` and `.venv-check/`, notebook checkpoints, and model-cache backups ending in `.joblib.bak` are ignored. The main notebook, generated CSVs/PNGs, and active `.joblib` caches remain versioned; model caches use Git LFS.

## Get the project and input data

You need Git, Python, and a notebook frontend such as VS Code with its Python and Jupyter extensions. Install Git LFS if you want to download the saved model caches.

Run the following from the directory where you want to place the project. The CheMixHub checkout goes inside the project's `Data/` folder. The project clone uses the `chenyue-work` branch. With Git LFS installed, cloning can download approximately 1.59 GB of model caches.

```powershell
git lfs install
git clone --branch chenyue-work https://github.com/nami-nia/CHE1147-Data-Science-Project.git project_code_results

New-Item -ItemType Directory -Force -Path project_code_results/Data | Out-Null
git clone --filter=blob:none --no-checkout https://github.com/chemcognition-lab/chemixhub.git project_code_results/Data/chemixhub
git -C project_code_results/Data/chemixhub sparse-checkout init --cone
git -C project_code_results/Data/chemixhub sparse-checkout set datasets/drug-solubility
git -C project_code_results/Data/chemixhub checkout 356c16b99e38b56a91b5219d33db410fc372d538
```

The sparse checkout retrieves the drug-solubility dataset without checking out the other dataset directories. The pinned CheMixHub revision is the revision inspected for this project. If you already have a CheMixHub checkout with other work, use a separate copy to obtain this revision.

The notebook expects this arrangement; the project directory itself can have a different name, but the data checkout must be at `Data/chemixhub/` inside it. This external checkout is excluded by `.gitignore`:

```text
project_code_results/
├── notebook_and_results/
│   └── project_notebook.ipynb
└── Data/
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

Use **Python 3.14.7**, the version selected by [.python-version](.python-version) and recorded in the saved notebook. [pyproject.toml](pyproject.toml) permits Python `>=3.14.7,<3.15` and declares the same ten pinned dependencies as [requirements.txt](requirements.txt): the scientific libraries, RDKit, scikit-learn, XGBoost, joblib, and an IPython kernel. These versions match the existing Windows environment. A fresh installation and execution on other platforms have not been verified.

With **uv** installed, run these commands from the project repository root (`project_code_results/`):

```powershell
uv sync
uv pip check
```

`uv sync` uses `.python-version`, installs the project and its dependencies into a project-local `.venv/`, and creates or updates `uv.lock`. No `uv.lock` is currently included; keep it versioned when generated to record the resolved dependency versions. The exact pins in `pyproject.toml` and `requirements.txt` cover the listed packages, not all transitive dependencies. Keep the two dependency lists aligned when changing the environment.

Alternatively, use **pip** from the project repository root. Confirm that `python --version` reports Python 3.14.7 first; pip does not select an interpreter using `.python-version`:

```powershell
python --version
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip check
```

On macOS/Linux, the environment interpreter is `.venv/bin/python` instead of `.venv\Scripts\python.exe`. The requirements provide the notebook kernel; VS Code or another Jupyter frontend supplies the notebook interface.

`src/data_cleaning_chenyue/` currently contains a package scaffold. There is no configured command-line entry point for running the analysis; use the notebook.

In VS Code, open the project folder and the final notebook, then select the newly created `.venv` Python interpreter as its kernel. If you use an existing Jupyter frontend, register the kernel with:

```powershell
.\.venv\Scripts\python.exe -m ipykernel install --user --name che1147 --display-name "CHE1147"
```

Then select **CHE1147** in that frontend. Installing packages in a different environment does not make them available to the selected notebook kernel.

## Run preprocessing, training, and analysis

Open `notebook_and_results/project_notebook.ipynb`, restart its kernel, and run all cells from top to bottom. Use the project root (`project_code_results/`) or `notebook_and_results/` as the notebook working directory. The notebook discovers its folder, checks the external inputs, creates `notebook_and_results/results/` when needed, and overwrites generated CSVs and PNGs there. Cell 1 resolves `Data/chemixhub/` relative to the project root.

Cell numbers below refer to positions in the notebook, not execution counters.

| Cells | Workflow stage |
| --- | --- |
| 1–9 | Configure paths, read the processed input, and inspect the original data and distributions. |
| 10–17 | Compare raw/processed records, repair solvent slots and fractions, and remove duplicate observations using active compositions. |
| 18–20 | Flag outliers within drug–solvent-pair groups, apply the retention policy, and export cleaned data and the policy report. Grouping uses recorded solvent IDs, including zero-fraction solvents and excluding `"N/A"` slots. |
| 21–24 | Generate cleaned-data plots and Pearson/Spearman correlation heatmaps for binary mixtures. |
| 25–27 | Split cleaned observations and construct 27 numeric features. |
| 28–32 | Train linear regression, define and tune Ridge and LASSO using five-fold cross-validation, and compare their test metrics. |
| 33–34 | Load saved XGBoost/Random Forest results or run their hyperparameter searches when caches are absent. |
| 35–37 | Evaluate the tree models on the test set and regenerate the five-model comparison table and two-panel MAE/RMSE bar chart. |
| 38 | Generate observed-versus-predicted parity plots for all five models on the held-out test set. |

Cell 29 defines `Ridge_Model` before running its grid search; cell 30 evaluates the selected Ridge model. Cell 31 defines `Lasso_Model`, runs its grid search, and evaluates the selected LASSO model. Both searches depend on the imports in cell 1, the split and target name in cell 25, and the features and preprocessor in cell 27. Their evaluation also uses `model_performance_metrics` from cell 28. Run these cells in order in the active kernel; saved outputs do not restore Python variables after a kernel restart.

For the pinned inputs, the saved cleaning run starts with **27,166 observations**, repairs **453 list-length mismatches**, expands **2,490 single-solvent representations**, and removes **9 duplicates**, leaving **27,157 observations**. The run flags **311 observations** using IQR rules and retains all of them. An IQR flag alone does not establish an erroneous measurement; rows are removed only when explicitly entered in `confirmed_outlier_errors` with supporting evidence. That dictionary is currently empty.

The target is the existing reported log-solubility `value`; the notebook does not apply an additional logarithm. Features include temperature, two solvent fractions, and eight RDKit descriptors for each of the drug and two solvent slots. An absent `"N/A"` solvent slot receives zero descriptors. The common preprocessor scales temperature and mole fractions; Ridge and LASSO additionally scale the resulting feature matrix in their own pipelines.

The split uses 80% training and 20% test observations with seed `1`, stratified using an indicator for target values above the dataset's 90th percentile. With the pinned inputs, this yields **21,725 training** and **5,432 test** observations. Hyperparameter selection uses the training set; MAE and RMSE are reported in the target's log-solubility scale. This is an observation-level split, so the test results do not establish generalization to entirely unseen drugs or solvent systems. The notebook creates its own split rather than using upstream split files or raw-data training labels.

### Use the saved model caches

The notebook loads the two `.joblib` caches from `notebook_and_results/results/`. The `.gitattributes` rules cover `.joblib` files in both the current and legacy folders. If your checkout contains Git LFS pointers, download the available LFS objects from the project root:

```powershell
git lfs install
git lfs pull
```

If the downloaded caches are under the legacy `data_cleaning_CL/results/` folder, copy the actual model files into `notebook_and_results/results/` before running the notebook. The notebook does not load caches from the legacy folder.

`random_forest_tuning_results_cleaned.joblib` is approximately 1.57 GB; `xgboost_tuning_results_cleaned.joblib` is approximately 19 MB. Allow memory for loading the models in addition to their file sizes. Run the whole notebook even when using caches: preprocessing, feature generation, linear-model fitting, evaluation, and plotting still run. Each cache contains the fitted pipeline, selected hyperparameters, and training/CV metrics; test metrics are calculated in cell 35.

The current code checks cache existence and does not validate whether it matches changed data, features, or split settings. Use the captured package versions for loading saved models, and retrain after changing those settings.

### Retrain XGBoost and Random Forest

Before running the notebook, move either or both cache files outside `results/`, or append `.bak` to the complete filename (for example, `random_forest_tuning_results_cleaned.joblib.bak`). These backup files are ignored by Git:

- `notebook_and_results/results/xgboost_tuning_results_cleaned.joblib`
- `notebook_and_results/results/random_forest_tuning_results_cleaned.joblib`

Restart the kernel and run all cells. A missing cache triggers its search and saves a replacement under the original filename. If your checkout has LFS pointer files instead of downloaded models, move those pointers aside as well; `joblib.load()` cannot read them as models.

XGBoost evaluates 200 randomized parameter combinations with five-fold CV; Random Forest evaluates 100 combinations with ten-fold CV. Each search therefore performs 1,000 CV fits plus the final refit. Both use seed `42`. Full tuning can take substantial time; no fixed runtime is assumed. Ridge and LASSO are fitted and tuned on every full run.

## Regenerate and inspect results

Running the final notebook from top to bottom regenerates these outputs in `notebook_and_results/results/` (model caches are written only when absent):

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
| `xgboost_tuning_results_cleaned.joblib` | Cell 33, when cache is absent | Selected fitted XGBoost pipeline and search/training summary. |
| `random_forest_tuning_results_cleaned.joblib` | Cell 34, when cache is absent | Selected fitted Random Forest pipeline and search/training summary. |
| `test_model_comparison.png` | Cell 36 | Five-model test MAE/RMSE comparison table. |
| `test_mae_rmse_comparison.png` | Cell 37 | Two-panel test MAE/RMSE bar chart saved at 300 dpi, with percentage reductions relative to Linear Regression and XGBoost highlighted as the best-performing model. |
| `Parity_Plots.png` | Cell 38 | Five observed-versus-predicted test-set scatterplots with a perfect-prediction reference line, shared axis limits, equal aspect ratios, and each model's MAE/RMSE; saved at 300 dpi. |

The bar chart uses Linear Regression as the baseline and displays each model's error value and percentage reduction: `(Baseline Error − Model Error) / Baseline Error × 100`. Negative percentages indicate higher error than the baseline. Both panels use the same vertical scale. The XGBoost highlight and `(best)` label are hard-coded in cell 37; they match the saved test results but are not selected dynamically after rerunning.

Other diagnostic tables, search summaries, and metrics are displayed in notebook outputs. Save the executed notebook to retain them.

Use regenerated outputs to report results after any code or input change. Exact strip-plot pixel positions may vary because its random jitter is generated before the notebook sets the global NumPy seed. The current notebook contains 38 code cells with saved execution counts 1–38, saved outputs through the parity plots, and no saved error outputs; its kernel metadata reports Python 3.14.7. Saved outputs do not establish a fresh environment verification or show whether the tree-model searches ran or their caches were loaded.

## Troubleshooting

- **Missing CSVs:** check the project's `Data/chemixhub/` folder, the pinned revision, and the three paths listed in cell 1. Do not substitute the final cleaned output for the original processed input.
- **Notebook not found:** use `project_code_results/` or `notebook_and_results/` as the working directory and keep the filename `project_notebook.ipynb`.
- **Import errors:** confirm that the notebook kernel is the environment where `requirements.txt` was installed.
- **Python version mismatch:** use Python 3.14.7 for the saved model caches. The uv workflow reads `.python-version`; the pip workflow uses the interpreter that created `.venv/`.
- **Model-loading errors:** check that `notebook_and_results/results/` contains actual model files rather than Git LFS pointers and that the installed package versions match the requirements. To regenerate a cache, move it aside and rerun the notebook as described above.
- **Undefined variables after running a later cell:** restart the kernel and execute all preceding cells in order.
- **`Ridge_Model` or `Lasso_Model` is not defined:** run the complete Ridge tuning cell (29) or LASSO tuning cell (31), including the pipeline definition above `GridSearchCV`, after the prerequisite cells. The current notebook defines each estimator in its own search cell.
