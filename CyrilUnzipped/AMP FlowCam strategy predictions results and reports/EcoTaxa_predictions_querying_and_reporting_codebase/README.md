# Flowcam Paper Code

This repository contains Python scripts and Jupyter Notebooks for analyzing EcoTaxa zooplankton classification predictions on DFO Flowcam data. It supports exploring different learning strategies across four aquaculture regions.

---

## Repository Structure

- **building_prediction_scenario_settings.py**  
  Generates and manages prediction scenario settings, including training/discarded class lists and strategy configurations. Produces scenario dictionaries used for batch predictions.

- **ecotaxa_api_requests.py**  
  Handles all interactions with the EcoTaxa API, including authentication, submitting prediction/export jobs, and downloading results. Uses scenario settings from `building_prediction_scenario_settings.py`.

- **reporting_prediction_results_w_md.py**  
  Main analysis and reporting script (also available as a Jupyter Notebook: `reporting_prediction_results_w_md.ipynb`). Loads, processes, and visualizes prediction results, generates summary metrics, and exports reports/figures.

- **plotting_helpers.py**  
  Contains plotting utilities for visualizing classification results, metrics, and strategy comparisons.

- **reporting_helpers.py**  
  Provides functions for processing prediction results, redistributing general categories, and preparing data for reporting.

- **scenarios_details_JSON/scenarios_implementation_details.json**  
  Stores all scenario configurations, including class lists and strategy metadata.

- **relAbundanceLabels.csv**  
  Maps full sample IDs to short labels for reporting.

- **taxon class names mapping master.xlsx**  
  Excel file mapping taxonomic class names between EcoTaxa and project-specific nomenclature.

- **exported_classification_results/**  
  (Ignored by git) Directory for downloaded EcoTaxa prediction/export results.

- **pickled_processed_exports/**  
  (Ignored by git) Stores pickled intermediate data for fast reloads.

- **strategy_predictions_results_and_reports/**  
  (Ignored by git) Output directory for generated figures, Excel reports, and summary files.

---

## Usage Workflow

1. **Configure Prediction Scenarios**
   - Edit or generate scenario settings using `building_prediction_scenario_settings.py`.
   - Scenario details are stored in `scenarios_details_JSON/scenarios_implementation_details.json`.

2. **Run EcoTaxa Predictions**
   - Use `ecotaxa_api_requests.py` to authenticate, submit prediction jobs, and download results.
   - Results are saved in `exported_classification_results/`.

3. **Process and Analyze Results**
   - Use `reporting_prediction_results_w_md.py` (or the notebook) to:
     - Load and process prediction results.
     - Map taxonomic categories and sample IDs.
     - Generate and save figures, Excel summaries, and reports in `strategy_predictions_results_and_reports/`.
     - Optionally, use pickled files in `pickled_processed_exports/` for faster reruns.

4. **Visualization and Reporting**
   - Figures and reports are automatically generated and saved.
   - Customize visualizations via `plotting_helpers.py` and `reporting_helpers.py` if needed.

---

## Requirements

- Python 3.10–3.12
- Install dependencies with [Poetry](https://python-poetry.org/):
  ```sh
  poetry install
  ```

---

## Notes

- Sensitive or large output folders are git-ignored.
- For EcoTaxa API access, ensure credentials are set up as required by `ecotaxa_py_client`.
- See code comments and docstrings for detailed parameter and function usage.

---

## Quick Start

1. **Install dependencies:**  
   `poetry install`

2. **Configure scenarios:**  
   Run or edit `building_prediction_scenario_settings.py`.

3. **Run predictions:**  
   Execute `ecotaxa_api_requests.py` to submit and download jobs.

4. **Analyze and report:**  
   Run `reporting_prediction_results_w_md.py` or open `reporting_prediction_results_w_md.ipynb` in Jupyter.

---

## Data Resource Files

This project relies on several key resource files for configuration, mapping, and reporting. Below are detailed descriptions of the main `.xlsx` and `.csv` files used:

### `taxon class names mapping master.xlsx`
- **Purpose:**  
  Provides standardized taxonomic category mappings between EcoTaxa output and project-specific nomenclature.
- **Usage:**  
  - Loaded in both `reporting_prediction_results_w_md.py` and the Jupyter notebook.
  - The sheet `'newname ecotaxa mapping'` is specifically used.
  - Columns utilized: `"newName"`, `"Ecotaxa dfo category specific"`, and `"Copepoda_mapping"`.
  - Duplicate rows are removed to ensure unique mapping entries.
- **Role:**  
  Ensures consistent taxonomic naming across all analyses and reports, allowing for accurate aggregation and comparison of results.

### `relAbundanceLabels.csv`
- **Purpose:**  
  Maps full sample identifiers to short, human-friendly labels.
- **Usage:**  
  - Imported as a DataFrame (`selected_samples_short_ID_mapping`) in `reporting_prediction_results_w_md.py` and the notebook.
  - Used for labeling in reports and visualizations, simplifying sample identification.
- **Role:**  
  Facilitates clear and concise reporting by providing readable sample names in figures and summary tables.

### `regional_prediction_strategies.xlsx`
- **Purpose:**  
  Stores configuration details for different prediction strategies and scenarios.
- **Usage:**  
  - Referenced in `building_prediction_scenario_settings.py` as `strategies_config_excel_filepath`.
  - Defines which classes are included/excluded, region-specific settings, and other scenario metadata.
- **Role:**  
  Centralizes scenario management, ensuring reproducibility and transparency in how prediction strategies are constructed and applied.

---

**Note:**  
All these files must be present in the project root or referenced paths for the scripts to function correctly. Their structure and content are critical for the reproducibility and accuracy of the workflow.