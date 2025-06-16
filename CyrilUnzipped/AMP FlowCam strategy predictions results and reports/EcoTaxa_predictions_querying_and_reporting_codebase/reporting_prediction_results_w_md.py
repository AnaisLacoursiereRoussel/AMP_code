# %% [markdown]
# #### Imports and Setup
# Importation of necessary libraries and custom modules.
# - `pickle` for object serialization/deserialization.
# - `pathlib.Path` for OS-independent path manipulations.
# - `pandas` for data manipulation and analysis.
# - `matplotlib.pyplot` for plotting.
# - Custom helper functions from local modules (`collate_models_metrics`, `collecting_results`, `plotting_helpers`,
# `reporting_helpers`) for specific tasks related to processing, visualizing, and reporting classification results.

# %%
import pickle
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

from collate_models_metrics import export_metrics_to_excel
from collecting_results import save_prediction_dataframes_to_excel
from plotting_helpers import (
    create_f1_scores_strategy_comp_by_scn_maxObj_figs,
    create_region_scn_maxObj_comp_strategy_metrics_figs,
)
from reporting_helpers import (
    batch_process_ecotaxa_prediction_jobs,
    create_classifier_metrics_dataframe,
    get_prediction_job_files_by_exported_ID_dict,
    process_all_export_jobs,
    process_all_visualizations,
)

# %% [markdown]
# #### Matplotlib Configuration
# `matplotlib` plotting parameters are established for consistent styling across all generated figures.
# The default font size, font family (iwona), LaTeX for text rendering, and a LaTeX preamble with necessary
# packages for the chosen font and encoding are specified.

# %%
params = {
    "font.size": 12,
    "font.family": "iwona",
    "text.usetex": True,
    "text.latex.preamble": r"""
    \usepackage[condensed,math]{iwona}
    \usepackage[T1]{fontenc}
    """,
}
plt.rcParams.update(params)

# %% [markdown]
# #### Define General Classes for Redistribution
# A set named `general_classes_to_be_redistributed` is defined. It contains strings representing broad
# or unidentified taxonomic classifications (e.g., "Calanoida (unid)", "Zooplankton (unid)").
# These classes may require special handling, such as redistribution to more specific categories,
# during later data processing steps.

# %%
general_classes_to_be_redistributed = {
    "Calanoida (unid)",
    "Calanoida (civ-vi)",
    "Cyclopoida (unid)",
    "Zooplankton (unid)",
}

# %% [markdown]
# #### Define Root Folder for Exported Results
# The variable `exported_results_root_folderpath` stores the string path to the root directory
# where all exported classification results from EcoTaxa are located. This path facilitates access
# to the raw prediction data.

# %%
exported_results_root_folderpath = r"exported_classification_results"

# %% [markdown]
# #### Retrieve Prediction Job Files
# The `get_prediction_job_files_by_exported_ID_dict` function from the `reporting_helpers` module is invoked.
# This function scans the `exported_results_root_folderpath` to identify and catalog all prediction job files,
# returning a dictionary, `exported_jobs_dict`, where keys are unique export job IDs and values are the file paths
# to the corresponding job data.

# %%
exported_jobs_dict = get_prediction_job_files_by_exported_ID_dict(exported_results_root_folderpath)

# %% [markdown]
# #### Define Columns to Retain
# A list called `column_to_keep` specifies the names of the columns to be selected and retained from the
# raw EcoTaxa export files during data processing. These columns typically include object identifiers,
# annotation details (original, new, EcoTaxa), sample information, and classifier output details
# (ID, name, score, timestamp).

# %%
column_to_keep = [
    "object_id",
    "object_annotation_status",
    "object_annotation_category",
    "object_annotation_hierarchy",
    "object_original_annotation_category",
    "object_newname",
    "object_stage",
    "object_ecotaxa_annotation_category",
    "sample_id",
    "sample_folder_name",
    "sample_region-year",
    "classif_id",
    "classif_auto_id",
    "classif_auto_name",
    "classif_auto_score",
    "classif_auto_when",
]

# %% [markdown]
# #### Load Taxonomic Category Mapping
# Taxonomic category mapping data is read from an Excel file named `taxon class names mapping master.xlsx`.
# Specifically, the 'newname ecotaxa mapping' sheet is loaded, and the "newName", "Ecotaxa dfo category specific",
# and "Copepoda_mapping" columns are selected. Duplicate rows are then removed from the resulting pandas DataFrame,
# `taxo_category_mapping_df`, ensuring unique mapping entries for standardizing taxonomic names.

# %%
taxo_category_mapping_df = pd.read_excel(
    r"taxon class names mapping master.xlsx",
    sheet_name="newname ecotaxa mapping",
    usecols=["newName", "Ecotaxa dfo category specific", "Copepoda_mapping"],
)
taxo_category_mapping_df.drop_duplicates(inplace=True)

# %% [markdown]
# #### Load Sample Short ID Mapping
# A CSV file named `relAbundanceLabels.csv` is imported into a pandas DataFrame called
# `selected_samples_short_ID_mapping`. This file is expected to contain mappings between full sample
# identifiers and shorter, more convenient labels, potentially for report labels or simplified sample
# identification in visualizations.

# %%
# Import selected sample's short ID mapping csv file
selected_samples_short_ID_mapping = pd.read_csv(r"relAbundanceLabels.csv")

# %% [markdown]
# #### Batch Process EcoTaxa Prediction Jobs
# Orchestration of the batch processing for the previously identified EcoTaxa prediction jobs.
# 1. A list `export_job_ids_list` is created, containing all job IDs from `exported_jobs_dict`.
# 2. The `batch_process_ecotaxa_prediction_jobs` function is then called. This function utilizes the list of job IDs,
# the dictionary of job files, the list of columns to keep, the taxonomic mapping DataFrame, and the root folder path
# for exported results.
# 3. Two dictionaries are returned:
#     - `exported_jobs_df_dict`: Containing processed data for each job as pandas DataFrames.
#     - `exported_jobs_settings_dict`: Containing settings or metadata associated with each processed job.

# %%
export_job_ids_list = list(exported_jobs_dict.keys())

exported_jobs_df_dict, exported_jobs_settings_dict = batch_process_ecotaxa_prediction_jobs(
    export_job_ids_list,
    exported_jobs_dict,
    column_to_keep,
    taxo_category_mapping_df,
    exported_results_root_folderpath,
)

# %% [markdown]
# #### Define Pickle File Paths for Processed Data
# File paths for saving (pickling) the results from the batch processing step are defined.
# - `exported_jobs_df_dict_pickle_filepath`: Path for storing the dictionary of processed job DataFrames.
# - `exported_jobs_settings_dict_pickle_filepath`: Path for storing the dictionary of job settings.
# These files will reside in a subdirectory named `pickled_processed_exports`. Using pickle files enables faster
# reloading of processed data in subsequent script runs, circumventing reprocessing.

# %%
exported_jobs_df_dict_pickle_filepath = (
    Path("pickled_processed_exports") / "exported_jobs_df_dict.pkl"
)
exported_jobs_settings_dict_pickle_filepath = (
    Path("pickled_processed_exports") / "exported_jobs_settings_dict.pkl"
)

# %% [markdown]
# #### Save Processed Data to Pickle Files
# The `exported_jobs_df_dict` (containing processed DataFrames for each job) and
# `exported_jobs_settings_dict` (containing settings for each job) are saved to disk
# using the `pickle` library. The data is serialized and written in binary mode (`'wb'`)
# to the file paths defined previously, allowing for efficient storage and retrieval of these complex Python objects.

# %%
with open(exported_jobs_df_dict_pickle_filepath, "wb") as f:
    pickle.dump(exported_jobs_df_dict, f)

with open(exported_jobs_settings_dict_pickle_filepath, "wb") as f:
    pickle.dump(exported_jobs_settings_dict, f)


# %% [markdown]
# #### Define Pickle File Path for Combined Data and Metrics
# The file path `exported_jobs_data_and_metrics_dict_pickle_filepath` is defined.
# This path points to a pickle file within the `pickled_processed_exports` directory,
# intended for storing a dictionary that contains both the processed data and calculated
# metrics for all exported jobs. This consolidated dictionary aims to facilitate easier
# access to all relevant information for subsequent analysis and reporting.

# %%
exported_jobs_data_and_metrics_dict_pickle_filepath = (
    Path("pickled_processed_exports") / "exported_jobs_data_and_metrics_dict.pkl"
)

# %% [markdown]
# #### Process All Export Jobs for Data and Metrics
# Invocation of the `process_all_export_jobs` function. This function receives the dictionaries of job DataFrames
# (`exported_jobs_df_dict`) and job settings (`exported_jobs_settings_dict`), the sample ID mappings
# (`selected_samples_short_ID_mapping`), and the file path for the output pickle file
# (`exported_jobs_data_and_metrics_dict_pickle_filepath`).
# These inputs are processed to compute various metrics and further refine the data for each job. The result is a
# comprehensive dictionary, `exported_jobs_data_and_metrics_dict`. The function may also handle caching by loading
# from or saving to the specified pickle file.

# %%
exported_jobs_data_and_metrics_dict = process_all_export_jobs(
    exported_jobs_df_dict,
    exported_jobs_settings_dict,
    selected_samples_short_ID_mapping,
    exported_jobs_data_and_metrics_dict_pickle_filepath,
)

# %% [markdown]
# #### Load Pickled Job DataFrames (Potentially Redundant)
# The `exported_jobs_df_dict` is loaded from the pickle file specified by `exported_jobs_df_dict_pickle_filepath`.
# This action reads the previously saved dictionary of processed job DataFrames from disk. This may be redundant
# if preceding operations ensure the dictionary is already in memory and current, or it could serve scripts run
# in segments or ensure data consistency if prior steps are omitted.

# %%
with open(exported_jobs_df_dict_pickle_filepath, "rb") as f:
    exported_jobs_df_dict = pickle.load(f)

# %% [markdown]
# #### Load Pickled Job Data and Metrics Dictionary
# The `exported_jobs_data_and_metrics_dict` is loaded from the pickle file specified by
# `exported_jobs_data_and_metrics_dict_pickle_filepath`. This retrieves the comprehensive dictionary
# containing detailed data and calculated metrics for each processed export job, making it available
# for subsequent visualization and reporting tasks.

# %%
with open(exported_jobs_data_and_metrics_dict_pickle_filepath, "rb") as f:
    exported_jobs_data_and_metrics_dict = pickle.load(f)

# %% [markdown]
# #### Load Pickled Job Settings Dictionary
# The `exported_jobs_settings_dict` is loaded from the pickle file specified by
# `exported_jobs_settings_dict_pickle_filepath`. This restores the dictionary containing
# the settings associated with each export job, potentially needed for context in visualizations or reports.

# %%
with open(exported_jobs_settings_dict_pickle_filepath, "rb") as f:
    exported_jobs_settings_dict = pickle.load(f)

# %% [markdown]
# #### Define Output Folder and Figure Display Flag
# Parameters for generating outputs are established:
# - `results_and_reports_folderpath`: A `Path` object pointing to the directory named
# `strategy_predictions_results_and_reports`. This folder is designated for storing all
# generated figures, reports, and Excel files.
# - `show_figure`: A boolean variable initialized to `False`. This flag likely governs whether
# plots generated by subsequent functions are displayed interactively during script execution
# or solely saved to files.

# %%
results_and_reports_folderpath = Path("strategy_predictions_results_and_reports")
show_figure = False

# %% [markdown]
# #### Process and Generate All Visualizations
# The `process_all_visualizations` function is called to generate and save a suite of plots and
# reports based on the processed data.
# It receives:
# - `exported_jobs_df_dict`: DataFrames for each job.
# - `exported_jobs_data_and_metrics_dict`: Combined data and metrics.
# - `exported_jobs_settings_dict`: Job settings.
# - `selected_samples_short_ID_mapping`: Sample ID mappings.
# - `taxo_category_mapping_df`: Taxonomic mappings.
# - `results_and_reports_folderpath`: The output directory for saving visualizations.
# - `general_classes_to_be_redistributed`: Set of general classes for special handling.
# - `show_figure`: Flag to control interactive plot display.
# The count of successfully processed jobs is returned and then printed.

# %%
successful_jobs = process_all_visualizations(
    exported_jobs_df_dict,
    exported_jobs_data_and_metrics_dict,
    exported_jobs_settings_dict,
    selected_samples_short_ID_mapping,
    taxo_category_mapping_df,
    results_and_reports_folderpath,
    general_classes_to_be_redistributed,
    show_figure,
)
print(f"Successfully processed {successful_jobs} export jobs")


# %% [markdown]
# #### Create Classifier Metrics DataFrame
# Invocation of the `create_classifier_metrics_dataframe` function. This function uses
# `exported_jobs_data_and_metrics_dict` (which contains detailed metrics for each job)
# to aggregate this information into a single, summary pandas DataFrame named `classifier_strategies_metrics_df`.
# This DataFrame is expected to organize key performance metrics (e.g., F1-score, precision, recall) across different
# classifier strategies or configurations, thereby facilitating comparisons.

# %%
classifier_strategies_metrics_df = create_classifier_metrics_dataframe(
    exported_jobs_data_and_metrics_dict
)


# %% [markdown]
# #### Generate Strategy Comparison Figures (Region, SCN, MaxObj)
# The `create_region_scn_maxObj_comp_strategy_metrics_figs` function is called.
# This function generates and saves figures comparing classifier strategy metrics, likely faceted or grouped by region,
# SCN (Scenario or Sample Class Name), and MaxObj (Maximum Objects).
# It utilizes:
# - `classifier_strategies_metrics_df`: The summary DataFrame of metrics.
# - `results_and_reports_folderpath`: The directory for saving the generated figures.
# - `show_figure=False`: To ensure figures are saved without interactive display.

# %%
create_region_scn_maxObj_comp_strategy_metrics_figs(
    classifier_strategies_metrics_df,
    results_and_reports_folderpath,
    show_figure=False,  # Change to True if you want to see the figures
)


# %% [markdown]
# #### Generate F1 Score Comparison Figures (SCN, MaxObj)
# The `create_f1_scores_strategy_comp_by_scn_maxObj_figs` function is invoked.
# This function is tasked with generating and saving figures that specifically compare F1 scores across different
# classifier strategies, detailed by SCN (Scenario or Sample Class Name) and MaxObj (Maximum Objects).
# It employs:
# - `classifier_strategies_metrics_df`: The DataFrame containing aggregated metrics, including F1 scores.
# - `results_and_reports_folderpath`: The output directory for the figures.
# - `show_figure=False`: To ensure figures are saved directly without interactive display.

# %%
create_f1_scores_strategy_comp_by_scn_maxObj_figs(
    classifier_strategies_metrics_df, results_and_reports_folderpath, show_figure=False
)


# %% [markdown]
# #### Save Prediction DataFrames to Excel
# The `save_prediction_dataframes_to_excel` function is called.
# Its purpose is to export the processed prediction data into Excel files, facilitating easier inspection and sharing.
# It requires:
# - `exported_jobs_df_dict`: Dictionary of DataFrames for each job.
# - `exported_jobs_data_and_metrics_dict`: Dictionary containing combined data and metrics.
# - `selected_samples_short_ID_mapping`: Mapping for sample IDs.
# - `results_and_reports_folderpath`: The directory where the Excel files will be saved.
# This provides a structured output of the detailed prediction results.

# %%
save_prediction_dataframes_to_excel(
    exported_jobs_df_dict,
    exported_jobs_data_and_metrics_dict,
    selected_samples_short_ID_mapping,
    results_and_reports_folderpath,
)

# %% [markdown]
# #### Export Aggregated Metrics to Excel
# The `export_metrics_to_excel` function is called.
# This function processes the `exported_jobs_data_and_metrics_dict`, which holds comprehensive data including
# calculated metrics for each job and strategy. It extracts and organizes these metrics into a tabular format
# suitable for an Excel spreadsheet. The resulting Excel file, summarizing performance metrics, is saved in the
# `results_and_reports_folderpath`.

# %%
export_metrics_to_excel(exported_jobs_data_and_metrics_dict, results_and_reports_folderpath)
