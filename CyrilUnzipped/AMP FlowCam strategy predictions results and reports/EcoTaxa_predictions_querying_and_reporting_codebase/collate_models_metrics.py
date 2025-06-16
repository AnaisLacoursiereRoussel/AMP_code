import os
from pathlib import Path
from typing import Any

import pandas as pd

from collecting_results import save_dataframes_to_excel_threaded


def create_metrics_excel_definition(
    job_id: int, metrics_data: dict[str, Any], index: int, results_base_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """
    Create metadata and DataFrame for exporting classification metrics to Excel.

    Parameters
    ----------
    job_id : int
        Unique identifier for the job/experiment.
    metrics_data : dict of str to Any
        Dictionary containing metrics and metadata for the job.
        Expected keys: 'region', 'use_scn_features', 'learning_objs_limit', 'strategy_numb', 'classif_report'.
    index : int
        Row index for the metadata DataFrame.
    results_base_path : Path
        Base directory where results should be saved.

    Returns
    -------
    excel_def : pandas.DataFrame
        DataFrame containing metadata for the Excel export.
    classif_report_df : pandas.DataFrame
        DataFrame representation of the classification report.
    sheet_name : str
        Name of the Excel sheet for this job's metrics.
    """
    region = metrics_data["region"]
    use_scn_features = metrics_data["use_scn_features"]
    learning_objs_limit = metrics_data["learning_objs_limit"]
    strategy_numb = metrics_data["strategy_numb"]
    classif_report = metrics_data["classif_report"]

    # Convert classification report to DataFrame
    classif_report_df = pd.DataFrame.from_dict(classif_report, orient="index")

    region_results_folderpath = results_base_path / region / "metrics"
    region_results_folderpath.mkdir(parents=True, exist_ok=True)

    # Format filenames with consistent naming convention
    scn_label = "With SCN" if use_scn_features else "Without SCN"
    obj_limit_label = f"{learning_objs_limit if learning_objs_limit < 999999 else 'Max'} Objects"
    excel_filename = f"{region}-{scn_label}-{obj_limit_label}-metrics data.xlsx"
    sheet_name = f"metrics_strategy#{strategy_numb}_expID{job_id}"

    # Create row in metadata dataframe
    excel_def = pd.DataFrame(
        {
            "region": [region],
            "scn_use": [scn_label],
            "max_objs": [obj_limit_label],
            "results_folderpath": [region_results_folderpath],
            "excel_filename": [excel_filename],
            "metrics_strategy_sheet_name": [sheet_name],
            "excel_filepath": [region_results_folderpath / excel_filename],
        },
        index=[index],
    )

    return excel_def, classif_report_df, sheet_name


def export_metrics_to_excel(
    exported_jobs_data_and_metrics_dict: dict[int, dict[str, Any]],
    results_and_reports_folderpath: Path,
) -> None:
    """
    Export classification metrics to Excel files, using threading for performance.

    Parameters
    ----------
    exported_jobs_data_and_metrics_dict : dict of int to dict of str to Any
        Dictionary mapping job IDs to their metrics and metadata.
    results_and_reports_folderpath : Path
        Base directory where Excel files and reports will be saved.

    Returns
    -------
    None
    """
    metrics_excel_files_definitions_df = pd.DataFrame()
    metrics_df_list = []

    # Create metadata and extract classification report DataFrames
    for i, (job_id, metrics_data) in enumerate(exported_jobs_data_and_metrics_dict.items()):
        excel_def, classif_report_df, sheet_name = create_metrics_excel_definition(
            job_id, metrics_data, i, results_and_reports_folderpath
        )

        metrics_excel_files_definitions_df = pd.concat(
            [metrics_excel_files_definitions_df, excel_def]
        )
        metrics_df_list.append(classif_report_df)

    # Add intermediate filepath column
    metrics_excel_files_definitions_df["intermediary_metrics_sheet_filepath"] = (
        metrics_excel_files_definitions_df.apply(
            lambda row: os.path.join(
                row["results_folderpath"], f"{row['metrics_strategy_sheet_name']}.xlsx"
            ),
            axis=1,
        )
    )

    # Create df_filepath_pairs for using save_dataframes_to_excel_threaded
    df_filepath_pairs = list(
        zip(
            metrics_df_list,
            metrics_excel_files_definitions_df["intermediary_metrics_sheet_filepath"].tolist(),
        )
    )

    # Use the threaded function from collecting_results.py to save DataFrames in parallel
    save_dataframes_to_excel_threaded(df_filepath_pairs)

    # Combine sheets into final Excel files
    unique_filepaths = metrics_excel_files_definitions_df["excel_filepath"].unique()
    for final_filepath in unique_filepaths:
        try:
            with pd.ExcelWriter(final_filepath) as writer:
                # Get all sheets for this filepath
                sheet_rows = metrics_excel_files_definitions_df[
                    metrics_excel_files_definitions_df["excel_filepath"] == final_filepath
                ][["intermediary_metrics_sheet_filepath", "metrics_strategy_sheet_name"]]

                # Add each sheet to the combined file
                for _, row in sheet_rows.iterrows():
                    temp_filepath = row["intermediary_metrics_sheet_filepath"]
                    sheet_name = row["metrics_strategy_sheet_name"]

                    try:
                        df = pd.read_excel(temp_filepath, index_col=0)
                        df.to_excel(writer, sheet_name=sheet_name)
                    except Exception as e:
                        print(f"Error processing {temp_filepath}: {str(e)}")
                    finally:
                        # Clean up temporary file
                        if os.path.exists(temp_filepath):
                            os.remove(temp_filepath)
        except Exception as e:
            print(f"Error creating combined file {final_filepath}: {str(e)}")
