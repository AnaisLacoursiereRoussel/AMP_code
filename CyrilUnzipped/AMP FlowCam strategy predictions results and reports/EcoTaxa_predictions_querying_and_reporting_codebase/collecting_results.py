import os
import pathlib
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Optional

import pandas as pd

from reporting_helpers import merge_sample_short_ID_to_exported_df


def save_df_to_excel(df: pd.DataFrame, filepath: str):
    """Save DataFrame to an Excel file"""
    df.to_excel(filepath)


def save_dataframes_to_excel_threaded(
    df_filepath_pairs: list[tuple[pd.DataFrame, str]], max_workers: Optional[int] = None
) -> None:
    """
    Save multiple DataFrames to Excel files using a thread pool.

    Args:
        df_filepath_pairs: List of tuples containing (DataFrame, filepath)
        max_workers: Maximum number of worker threads (None=auto-determine based on system)
    """
    # Use ThreadPoolExecutor to manage threads efficiently
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks to the executor
        futures = [
            executor.submit(save_df_to_excel, df, filepath) for df, filepath in df_filepath_pairs
        ]

        # Wait for all tasks to complete (optional - the context manager already does this)
        for future in futures:
            try:
                future.result()  # This will raise any exceptions that occurred during execution
            except Exception as e:
                print(f"Error saving DataFrame to Excel: {str(e)}")


def create_excel_export_metadata(
    export_job_id: Any, job_data: dict[str, Any], results_base_path: pathlib.Path
) -> dict[str, Any]:
    """
    Create metadata for Excel export including paths, filenames and sheet names.

    Args:
        export_job_id: The ID of the export job
        job_data: Dictionary containing job metadata
        results_base_path: Base path for storing results

    Returns:
        dict: Dictionary containing export metadata
    """
    region = job_data["region"]
    use_scn_features = job_data["use_scn_features"]
    learning_objs_limit = job_data["learning_objs_limit"]
    strategy_numb = job_data["strategy_numb"]

    # Create region directory
    region_results_folderpath = results_base_path / region / "results"
    region_results_folderpath.mkdir(parents=True, exist_ok=True)

    # Create file and sheet names
    max_objs_text = f"{learning_objs_limit if learning_objs_limit < 999999 else 'Max'} Objects"
    scn_text = "With SCN" if use_scn_features else "Without SCN"

    excel_filename = f"{region}-{scn_text}-{max_objs_text}-predictions data.xlsx"
    excel_filepath = region_results_folderpath / excel_filename

    retained_sheet_name = f"retain_cls_strat{strategy_numb}_expID{export_job_id}"
    discarded_sheet_name = f"discard_cls_strat{strategy_numb}_expID{export_job_id}"

    # Intermediate file paths (for temporary files)
    retained_temp_filepath = region_results_folderpath / f"{retained_sheet_name}.xlsx"
    discarded_temp_filepath = region_results_folderpath / f"{discarded_sheet_name}.xlsx"

    return {
        "region": region,
        "scn_use": scn_text,
        "max_objs": max_objs_text,
        "results_folderpath": region_results_folderpath,
        "excel_filename": excel_filename,
        "excel_filepath": excel_filepath,
        "retained_cls_sheet_name": retained_sheet_name,
        "discarded_cls_sheet_name": discarded_sheet_name,
        "intermediary_retain_sheet_filepath": retained_temp_filepath,
        "intermediary_discard_sheet_filepath": discarded_temp_filepath,
    }


def save_prediction_dataframes_to_excel(
    exported_jobs_df_dict: dict[Any, pd.DataFrame],
    exported_jobs_data_and_metrics_dict: dict[Any, dict[str, Any]],
    selected_samples_short_ID_mapping: pd.DataFrame,
    results_and_reports_folderpath: pathlib.Path,
) -> None:
    """
    Process prediction DataFrames and save them to Excel files organized by region and configuration.

    Args:
        exported_jobs_df_dict: Dictionary of exported job DataFrames
        exported_jobs_data_and_metrics_dict: Dictionary of job metadata and metrics
        selected_samples_short_ID_mapping: DataFrame with sample ID mapping
        results_and_reports_folderpath: Base path for results and reports
    """
    # Create DataFrame to track Excel file definitions and paths
    excel_files_definitions_df = pd.DataFrame()
    retained_classes_prediction_df_list = []
    discarded_classes_prediction_df_list = []

    # Step 1: Process each job and collect metadata
    for i, (export_job_id, prediction_df) in enumerate(list(exported_jobs_df_dict.items())):
        if prediction_df is None:
            continue

        # Process the prediction DataFrame
        prediction_df = merge_sample_short_ID_to_exported_df(
            selected_samples_short_ID_mapping, prediction_df
        )

        # Get job data and extract prediction DataFrames
        job_data = exported_jobs_data_and_metrics_dict[export_job_id]
        retained_classes_prediction_df = job_data["retained_classes_prediction_df"]
        discarded_classes_prediction_df = job_data["discarded_classes_prediction_df"]

        # Create export metadata
        export_meta = create_excel_export_metadata(
            export_job_id, job_data, results_and_reports_folderpath
        )

        # Store metadata in Excel definitions DataFrame
        for key, value in export_meta.items():
            excel_files_definitions_df.loc[i, key] = value

        # Store prediction DataFrames for later processing
        retained_classes_prediction_df_list.append(retained_classes_prediction_df)
        discarded_classes_prediction_df_list.append(discarded_classes_prediction_df)

    # Step 2: Save DataFrames to temporary files using threads (both retained and discarded classes)
    df_filepath_pairs: list[tuple[pd.DataFrame, str]] = []

    # Create pairs of (DataFrame, filepath) for retained classes
    df_filepath_pairs.extend(
        zip(
            retained_classes_prediction_df_list,
            excel_files_definitions_df["intermediary_retain_sheet_filepath"],
        )
    )

    # Create pairs of (DataFrame, filepath) for discarded classes
    df_filepath_pairs.extend(
        zip(
            discarded_classes_prediction_df_list,
            excel_files_definitions_df["intermediary_discard_sheet_filepath"],
        )
    )

    # Save all DataFrames in parallel
    save_dataframes_to_excel_threaded(df_filepath_pairs)

    # Step 3: Combine temporary files into final Excel files by region/configuration
    for excel_filepath in excel_files_definitions_df["excel_filepath"].unique():
        # Get all rows for this Excel file
        file_rows = excel_files_definitions_df[
            excel_files_definitions_df["excel_filepath"] == excel_filepath
        ][
            [
                "intermediary_retain_sheet_filepath",
                "retained_cls_sheet_name",
                "intermediary_discard_sheet_filepath",
                "discarded_cls_sheet_name",
            ]
        ]

        with pd.ExcelWriter(excel_filepath) as writer:
            sheets_written = 0

            for _, row in file_rows.iterrows():
                # Process retained and discarded files
                for file_col, name_col in [
                    ("intermediary_retain_sheet_filepath", "retained_cls_sheet_name"),
                    ("intermediary_discard_sheet_filepath", "discarded_cls_sheet_name"),
                ]:
                    try:
                        filepath = row[file_col]
                        sheet_name = row[name_col]

                        if os.path.exists(filepath):
                            df = pd.read_excel(filepath)
                            df.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written += 1
                            os.remove(filepath)
                        else:
                            print(f"Warning: Temporary file {filepath} does not exist")

                    except Exception as e:
                        print(f"Error processing {filepath}: {str(e)}")

            # Ensure workbook has at least one sheet
            if sheets_written == 0:
                pd.DataFrame().to_excel(writer, sheet_name="Empty", index=False)
