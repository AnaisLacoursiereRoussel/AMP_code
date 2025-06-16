import copy
import json
import os
import pickle
import zipfile
from functools import partial
from pathlib import Path
from random import choices
from typing import Any, Optional, Union

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from sklearn.metrics import classification_report, confusion_matrix

from plotting_helpers import (
    create_classification_figure,
    create_discarded_categories_heatmap,
    create_relative_abundance_figure,
)

broad_categ_strat_mapping = {
    ("No Calanoida, Cyclopoida, Zooplankton " "classes in learning set"): "NO Cal, Cycl, Zoo",
    (
        "With Calanoida, Cyclopoida, NO Zooplankton " "classes in learning set"
    ): "WITH Cal, Cycl, NO Zoo",
    (
        "With Calanoida, Cyclopoida and Zooplankton " "classes in learning set"
    ): "WITH Cal, Cycl, Zoo",
    (
        "No Anthoathecata, Calanoida, Copepoda, Zooplankton " "classes in learning set"
    ): "NO Antho, Cal, Cycl, Zoo",
    (
        "No Calanoida (civ-vi), Cyclopoida, Zooplankton " "classes in learning set"
    ): "NO Cal(4-6), Cycl, Zoo",
}


def read_first_tsv_from_zip(zip_path: str) -> Optional[pd.DataFrame]:
    """
    Reads the first TSV file from a ZIP archive and returns it as a pandas DataFrame.

    Parameters
    ----------
    zip_path : str
        Path to the ZIP file containing one or more TSV files.

    Returns
    -------
    df : Optional[pandas.DataFrame]
        DataFrame containing the data from the first TSV file found in the ZIP archive.
        Returns None if no TSV file is found.
    """
    with zipfile.ZipFile(zip_path, "r") as z:
        for filename in z.namelist():
            if filename.endswith(".tsv"):
                with z.open(filename) as f:
                    df = pd.read_csv(f, sep="\t")
                return df
    return None


def get_prediction_job_files_by_exported_ID_dict(
    export_folderpath: str,
) -> dict[int, tuple[str, Optional[str]]]:
    """
    Retrieves prediction job ZIP files and their corresponding JSON settings files from a folder.

    Parameters
    ----------
    export_folderpath : str
        Path to the folder containing exported prediction job files.

    Returns
    -------
    prediction_job_files_by_exported_ID_dict : dict[int, tuple[str, Optional[str]]]
        Dictionary mapping EcoTaxa exported job IDs to a tuple:
        (ZIP filename, JSON settings filename or None if not found).
    """
    exported_prediction_job_file_list: list[str] = os.listdir(export_folderpath)
    prediction_job_files_by_exported_ID_dict: dict[int, tuple[str, Optional[str]]] = {}
    for filename in exported_prediction_job_file_list:
        if filename.endswith(".zip"):
            export_id: int = int(filename[9:14])
            try:
                json_settings_filename: Optional[str] = [
                    filename
                    for filename in exported_prediction_job_file_list
                    if (
                        filename.startswith(f"exportjob{export_id}")
                        and filename.endswith("_settings.json")
                    )
                ][0]
            except IndexError:
                print(f"No JSON settings file for export job ID {export_id}")
                json_settings_filename = None
            prediction_job_files_by_exported_ID_dict[export_id] = (
                filename,
                json_settings_filename,
            )
    return prediction_job_files_by_exported_ID_dict


def sort_and_move_to_end_df_rows(
    df: pd.DataFrame, sort_column: str, move_column: str, values_to_end: list[str]
) -> pd.DataFrame:
    """
    Sorts a DataFrame by a specified column and moves rows with certain values in another column to the end.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame to sort and rearrange.
    sort_column : str
        Column name to sort by (descending).
    move_column : str
        Column whose values determine which rows are moved to the end.
    values_to_end : list of str
        List of values in move_column for which rows should be moved to the end.

    Returns
    -------
    df_final : pandas.DataFrame
        DataFrame sorted and rearranged as specified.
    """
    # Sort the dataframe by 'sort_column'
    df_sorted = df.sort_values(by=sort_column, ascending=False)

    # Split the rows you want to move to the end
    rows_to_end = df_sorted[df_sorted[move_column].isin(values_to_end)]
    rows_to_keep = df_sorted[~df_sorted[move_column].isin(values_to_end)]

    # Concatenate the two dataframes
    df_final = pd.concat([rows_to_keep, rows_to_end])

    return df_final


def process_ecotaxa_exported_job_data(
    prediction_job_export_filename: str,
    column_to_keep: list[str],
    taxo_category_mapping_df: pd.DataFrame,
    exported_results_root_folderpath: str,
) -> pd.DataFrame:
    """
    Loads and processes EcoTaxa exported job data, merging with taxonomic category mapping.

    Parameters
    ----------
    prediction_job_export_filename : str
        Filename of the EcoTaxa exported prediction job ZIP file.
    column_to_keep : list of str
        List of column names to retain from the exported data.
    taxo_category_mapping_df : pandas.DataFrame
        DataFrame containing taxonomic category mappings.
    exported_results_root_folderpath : str
        Root folder path where exported results are stored.

    Returns
    -------
    prediction_job_export_df : pandas.DataFrame
        Processed DataFrame with merged taxonomic mapping.

    Raises
    ------
    FileNotFoundError
        If no TSV file is found in the ZIP archive.
    """
    prediction_job_export_filepath = (
        exported_results_root_folderpath + r"\\" + prediction_job_export_filename
    )
    prediction_job_export_df = read_first_tsv_from_zip(prediction_job_export_filepath)
    if prediction_job_export_df is not None:
        prediction_job_export_df = prediction_job_export_df.loc[:, column_to_keep]
        prediction_job_export_df = prediction_job_export_df.merge(
            taxo_category_mapping_df,
            how="left",
            left_on="object_annotation_category",
            right_on="Ecotaxa dfo category specific",
        )
        prediction_job_export_df.drop(columns="Ecotaxa dfo category specific", inplace=True)
        prediction_job_export_df.rename(
            columns={
                "newName": "predicted_newName",
                "stage": "predicted_stage",
            },
            inplace=True,
        )

        return prediction_job_export_df
    else:
        raise FileNotFoundError("No tsv file was produced for the EcoTaxa exported job")


def batch_process_ecotaxa_prediction_jobs(
    export_job_ids_list: list[int],
    exported_jobs_dict: dict[int, tuple[str, Optional[str]]],
    column_to_keep: list[str],
    taxo_category_mapping_df: pd.DataFrame,
    exported_results_root_folderpath: str,
) -> tuple[dict[int, Optional[pd.DataFrame]], dict[int, Any]]:
    """
    Batch processes EcoTaxa prediction jobs, returning exported DataFrames and settings for each job ID.

    Parameters
    ----------
    export_job_ids_list : list of int
        List of EcoTaxa export job IDs to process.
    exported_jobs_dict : dict[int, tuple[str, Optional[str]]]
        Dictionary mapping export job IDs to (ZIP filename, JSON settings filename).
    column_to_keep : list of str
        List of column names to retain from each exported DataFrame.
    taxo_category_mapping_df : pandas.DataFrame
        DataFrame containing taxonomic category mappings.
    exported_results_root_folderpath : str
        Root folder path where exported results are stored.

    Returns
    -------
    exported_jobs_df_dict : dict[int, pandas.DataFrame or None]
        Dictionary mapping job IDs to their exported DataFrame (or None if not found).
    exported_jobs_settings_dict : dict[int, Any]
        Dictionary mapping job IDs to their settings loaded from JSON (or None if not found).
    """
    # Dictionary that should contain the exported prediction
    # dataframe for each exported job ID
    exported_jobs_df_dict: dict[int, Optional[pd.DataFrame]] = {}

    # Dictionary that should contain the prediction scenario settings
    # dataframe for each exported job ID
    exported_jobs_settings_dict: dict[int, Any] = {}

    for export_job_id in export_job_ids_list:
        (
            prediction_job_export_filename,
            prediction_job_settings_filename,
        ) = exported_jobs_dict.get(export_job_id, (None, None))

        if prediction_job_export_filename and prediction_job_settings_filename:
            try:
                prediction_job_export_df = process_ecotaxa_exported_job_data(
                    prediction_job_export_filename,
                    column_to_keep,
                    taxo_category_mapping_df,
                    exported_results_root_folderpath,
                )
                exported_jobs_df_dict[export_job_id] = prediction_job_export_df
            except FileNotFoundError:
                print(f"No data was exported for the job ID: {export_job_id}")
                exported_jobs_df_dict[export_job_id] = None

            prediction_job_settings_filepath = (
                exported_results_root_folderpath + r"\\" + prediction_job_settings_filename
            )
            with open(prediction_job_settings_filepath, "r") as json_settings_file:
                exported_jobs_settings = json.load(json_settings_file)
            exported_jobs_settings_dict[export_job_id] = exported_jobs_settings
        else:
            exported_jobs_df_dict[export_job_id] = None
            exported_jobs_settings_dict[export_job_id] = None

    return exported_jobs_df_dict, exported_jobs_settings_dict


def format_n_and_train_obj_string(
    row: pd.Series, n_col_label: str, train_col_label: str, objs_limit: int
) -> str:
    """
    Formats a string summarizing the number of objects and training objects for a category.

    Parameters
    ----------
    row : pandas.Series
        Row from a DataFrame containing the relevant columns.
    n_col_label : str
        Column label for the number of objects.
    train_col_label : str
        Column label for the number of training objects.
    objs_limit : int
        Maximum number of training objects to display.

    Returns
    -------
    summary_str : str
        Formatted string with category name, object count, and training object count.
    """
    n_value = 0 if pd.isna(row[n_col_label]) else int(row[n_col_label])
    train_value = (
        int(row[train_col_label])
        if not pd.isna(row[train_col_label]) and int(row[train_col_label]) < objs_limit
        else objs_limit
    )

    return f"{row['newName']}\n(n={n_value}-train={train_value})"


def get_samples_and_train_categ_counts_df(
    ecotaxa_training_classes: list[str],
    prediction_implementation_details_df: pd.DataFrame,
    training_set_type: str,
    learning_objs_limit: int,
):
    """
    Generates a DataFrame with sample and training category counts, including formatted labels.

    Parameters
    ----------
    ecotaxa_training_classes : list of str
        List of EcoTaxa training class names.
    prediction_implementation_details_df : pandas.DataFrame
        DataFrame with implementation details for predictions.
    training_set_type : str
        Type of training set ('Regional only' or other).
    learning_objs_limit : int
        Maximum number of training objects to display.

    Returns
    -------
    samples_and_train_categ_pred_data : pandas.DataFrame
        DataFrame with category counts, learned status, and formatted frequency labels.
    """
    if training_set_type == "Regional only":
        training_classes_count_col_name = "region_training_classes_count"
    else:
        training_classes_count_col_name = "all_regions_training_classes_count"

    # Retrieve the classes present in the selected samples
    prediction_implementation_details_df = prediction_implementation_details_df.loc[
        prediction_implementation_details_df["Category type"] != "Not Present", :
    ]

    samples_and_train_categ_pred_data = prediction_implementation_details_df.loc[
        :,
        [
            "newName",
            "stage",
            "ecotaxa_category",
            "Broad category",
            "Category type",
            "selected_samples_classes_count",
            training_classes_count_col_name,
        ],
    ]
    samples_and_train_categ_pred_data["learned_category"] = False

    samples_and_train_categ_pred_data.loc[
        samples_and_train_categ_pred_data["ecotaxa_category"].isin(ecotaxa_training_classes),
        "learned_category",
    ] = True

    # Filter out any category being "Training only" and learned_category == False
    samples_and_train_categ_pred_data = samples_and_train_categ_pred_data.loc[
        ~(
            (~samples_and_train_categ_pred_data["learned_category"])
            & (samples_and_train_categ_pred_data["Category type"] == "Training only")
        )
    ]

    samples_and_train_categ_pred_data["category_frequency_labels"] = (
        samples_and_train_categ_pred_data.apply(
            lambda row: format_n_and_train_obj_string(
                row,
                "selected_samples_classes_count",
                training_classes_count_col_name,
                learning_objs_limit,
            ),
            axis=1,
        )
    )

    samples_and_train_categ_pred_data.loc[
        ~samples_and_train_categ_pred_data["learned_category"], "category_frequency_labels"
    ] = samples_and_train_categ_pred_data.loc[
        ~samples_and_train_categ_pred_data["learned_category"]
    ].apply(
        lambda row: "$\\fontfamily{iwona}\\selectfont\\textbf{"
        + row["newName"]
        + "}$\n$\\fontfamily{iwona}\\selectfont\\textrm{n="
        + str(int(row["selected_samples_classes_count"]))
        + "}$",
        axis=1,
    )

    return samples_and_train_categ_pred_data


def replace_categ_by_values_from_weighted_distrib(
    value: str,
    replaced_category: str,
    replacing_categories,
    replacement_weights: Optional[list[float]],
):
    """
    Replaces a category value with a randomly chosen replacement based on provided weights.

    Parameters
    ----------
    value : str
        The category value to check and potentially replace.
    replaced_category : str
        The category to be replaced.
    replacing_categories : list
        List of categories to choose from as replacements.
    replacement_weights : list of float or None
        Weights for random selection among replacing_categories.

    Returns
    -------
    new_value : str
        The original value or a randomly chosen replacement if it matches replaced_category.
    """
    if value == replaced_category:
        return choices(replacing_categories, weights=replacement_weights, k=1)[0]
    return value


def redistribute_general_categories(
    prediction_df: pd.DataFrame, present_general_class_set: set[str], redistributed_col_label: str
) -> None:
    """
    Redistributes general categories in a DataFrame column based on hierarchical mappings and relative abundances.

    Parameters
    ----------
    prediction_df : pandas.DataFrame
        DataFrame containing predictions and category mappings.
    present_general_class_set : set of str
        Set of general class names to redistribute.
    redistributed_col_label : str
        Column label in prediction_df to redistribute.

    Returns
    -------
    None
        Modifies prediction_df in place by adding a new column with redistributed categories.
    """
    prediction_df[f"{redistributed_col_label}_redistr"] = prediction_df[redistributed_col_label]

    predicted_classes_count_df = prediction_df.groupby(
        [redistributed_col_label, "Copepoda_mapping"], dropna=False
    ).size()
    df = predicted_classes_count_df.to_frame().reset_index()[
        [redistributed_col_label, "Copepoda_mapping"]
    ]
    hierarchical_mappings = {
        row[redistributed_col_label]: row["Copepoda_mapping"].split(", ")
        for _, row in df.iterrows()
        if pd.notna(row["Copepoda_mapping"])
    }

    if "Zooplankton (unid)" in present_general_class_set:
        zoo_redistrib = True
        present_general_class_set.discard("Zooplankton (unid)")
    else:
        zoo_redistrib = False

    for present_general_class in present_general_class_set:
        replacing_categories = [
            k for k, v in hierarchical_mappings.items() if present_general_class in v
        ]

        relative_abundance = prediction_df.loc[
            prediction_df["Copepoda_mapping"].str.contains(
                present_general_class, na=False, regex=False
            ),
            redistributed_col_label,
        ].value_counts(normalize=True)
        replacement_weights = [
            relative_abundance.get(cat, 0) for cat in replacing_categories  # type: ignore
        ]

        replacement_function = partial(
            replace_categ_by_values_from_weighted_distrib,
            replaced_category=present_general_class,
            replacing_categories=replacing_categories,
            replacement_weights=replacement_weights,
        )

        prediction_df[f"{redistributed_col_label}_redistr"] = prediction_df[
            f"{redistributed_col_label}_redistr"
        ].apply(replacement_function)

    if zoo_redistrib:
        replacing_categories = list(
            prediction_df.loc[
                prediction_df[f"{redistributed_col_label}_redistr"] != "Zooplankton (unid)",
                f"{redistributed_col_label}_redistr",
            ].unique()
        )
        replacement_function = partial(
            replace_categ_by_values_from_weighted_distrib,
            replaced_category="Zooplankton (unid)",
            replacing_categories=replacing_categories,
            replacement_weights=None,
        )

        prediction_df[f"{redistributed_col_label}_redistr"] = prediction_df[
            f"{redistributed_col_label}_redistr"
        ].apply(replacement_function)


def merge_sample_short_ID_to_exported_df(
    selected_samples_short_ID_mapping: pd.DataFrame, exported_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merges short sample IDs into an exported DataFrame based on a mapping DataFrame.

    Parameters
    ----------
    selected_samples_short_ID_mapping : pandas.DataFrame
        DataFrame mapping sample IDs to short labels.
    exported_df : pandas.DataFrame
        DataFrame containing exported data with sample IDs.

    Returns
    -------
    merged_df : pandas.DataFrame
        DataFrame with an added 'short_sample_id' column.
    """
    exported_df = exported_df.merge(
        selected_samples_short_ID_mapping[["updatedLabel", "FlowcamCode"]],
        left_on="sample_id",
        right_on="FlowcamCode",
        how="left",
    )
    exported_df.drop(columns="FlowcamCode", inplace=True)
    exported_df.rename(columns={"updatedLabel": "short_sample_id"}, inplace=True)
    return exported_df


def correct_classif_rept_macro_avg_for_xtra_class_prev(
    report_dict: dict[str, dict[str, float]],
) -> dict[str, dict[str, float]]:
    """
    Recomputes the macro average of a classification report, excluding classes with zero support.

    Parameters
    ----------
    report_dict : dict of str to dict of str to float
        Dictionary containing classification report metrics for each class.

    Returns
    -------
    updated_report_dict : dict of str to dict of str to float
        Updated report dictionary with recalculated macro averages and removed zero-support classes.
    """

    macro_precision_sum, macro_recall_sum, macro_f1_sum = 0.0, 0.0, 0.0

    none_support_classes: dict[str, dict[str, float]] = {}
    for class_name, metrics in report_dict.items():
        if metrics["support"] == 0.0:
            none_support_classes[class_name] = metrics
            continue
        if "avg" in class_name:
            continue

        # Accumulate for macro averages
        macro_precision_sum += metrics["precision"]
        macro_recall_sum += metrics["recall"]
        macro_f1_sum += metrics["f1-score"]

    num_classes = (
        len(report_dict) - len(none_support_classes) - 3
    )  # Exclude 'avg' entries and excluded classes
    if num_classes == 0:
        return report_dict  # No valid classes to compute averages

    # Calculate new macro averages
    macro_avg = {
        "precision": macro_precision_sum / num_classes,
        "recall": macro_recall_sum / num_classes,
        "f1-score": macro_f1_sum / num_classes,
        "support": report_dict["macro avg"]["support"],
    }

    updating_dict = {"macro avg (corr)": macro_avg, "weighted avg": report_dict.pop("weighted avg")}

    # Update the report dictionary with new averages
    report_dict.pop("macro avg")
    report_dict.update(updating_dict)

    return report_dict


def correct_classif_rept_macro_avg_for_xtra_class(
    report_dict: dict[str, dict[str, float]],
) -> dict[str, dict[str, float]]:
    """
    Recomputes the Macro average of a classification report,
    excluding classes with support = 0.

    :param report_dict: Dictionary containing classification report metrics.
    :return: Updated report dictionary with recalculated averages and
                modified metrics for classes having support = 0.
    """

    # Initialize variables for macro and weighted averages
    macro_precision_sum, macro_recall_sum, macro_f1_sum = 0.0, 0.0, 0.0

    num_classes = 0
    for class_name, metrics in report_dict.items():
        if (metrics["support"] != 0.0) and (class_name not in ["macro avg", "weighted avg"]):
            num_classes += 1
            # Accumulate for macro averages
            macro_precision_sum += metrics["precision"]
            macro_recall_sum += metrics["recall"]
            macro_f1_sum += metrics["f1-score"]

    # Avoid division by zero
    if num_classes == 0:
        return report_dict  # No valid classes to compute averages

    # Calculate new macro averages
    corrected_macro_avgs = {
        "precision": macro_precision_sum / num_classes,
        "recall": macro_recall_sum / num_classes,
        "f1-score": macro_f1_sum / num_classes,
        "support": report_dict["macro avg"]["support"],
    }

    corrected_report_dict = copy.deepcopy(report_dict)
    updating_dict = {
        "macro avg (corr)": corrected_macro_avgs,
        "weighted avg": corrected_report_dict.pop("weighted avg"),
    }
    # Update the report dictionary with new averages
    corrected_report_dict.pop("macro avg")
    corrected_report_dict.update(updating_dict)

    # Modify the none support categories to display '-' instead of 0.0
    # for their associated scores

    return corrected_report_dict


def extract_job_settings(export_job_id: int, settings_dict: dict[int, Any]) -> dict[str, Any]:
    """
    Extract and organize all settings for a given export job.

    Parameters
    ----------
    export_job_id : int
        The ID of the export job whose settings are to be extracted.
    settings_dict : dict
        Dictionary containing settings for all export jobs, typically loaded from JSON.

    Returns
    -------
    settings : dict
        Dictionary containing organized settings for the specified export job, including
        region, training set type, strategy number, learning object limit, SCN feature usage,
        sample class presence, extra region training class flag, thresholds, strategy label,
        and other relevant configuration details.
    """
    job_settings = settings_dict[export_job_id][1]
    strategy_scenario = job_settings["strategy_scenario_record"]

    settings = {
        "region": job_settings["region"],
        "training_set_type": strategy_scenario["Training class set"],
        "strategy_numb": strategy_scenario["strategy_number"],
        "learning_objs_limit": job_settings["learning_objs_limit"],
        "use_scn_features": job_settings["use_scn_flags"],
        "sel_sples_present_classes": strategy_scenario["Select sples classes present"],
        "xtr_region_training_only_classes": strategy_scenario[
            "Extra regional training only classes"
        ]
        == "Yes",
        "low_training_threshold": strategy_scenario["Low training threshold"],
        "broad_categ_strategy": strategy_scenario["Broad classes strategy"],
        "scenario_title": job_settings["scenario_title"],
        "training_classes_list": job_settings["training_classes_list"],
        "discarded_training_classes_list": job_settings["discarded_training_classes_list"],
        "strategy_df": job_settings["strategy_df"],
    }

    # Compute strategy label
    broad_categ = broad_categ_strat_mapping.get(
        settings["broad_categ_strategy"], "No broad categ strat?"
    )
    settings["strategy_label"] = (
        f"Strategy {settings['strategy_numb']}\nsples: {settings['sel_sples_present_classes']}\n"
        f"Training={settings['training_set_type']}\n{broad_categ}"
    )
    if settings["xtr_region_training_only_classes"]:
        settings["strategy_label"] += "\nxtra train classes"

    return settings


def process_prediction_data(
    prediction_df: pd.DataFrame,
    settings: dict[str, Any],
    prediction_implementation_details_df: pd.DataFrame,
) -> dict[str, Any]:
    """
    Process prediction data to separate retained and discarded classes and organize related metadata.

    Parameters
    ----------
    prediction_df : pandas.DataFrame
        DataFrame containing prediction results for all objects in the export job.
    settings : dict
        Dictionary of settings for the export job, as returned by `extract_job_settings`.
    prediction_implementation_details_df : pandas.DataFrame
        DataFrame containing implementation details for the prediction scenario.

    Returns
    -------
    result : dict
        Dictionary with the following keys:
            - 'samples_and_train_categ_pred_data': DataFrame with sample and training category data.
            - 'retained_sples_and_train_categ_pred_data': DataFrame of retained (learned) categories.
            - 'retained_classes_prediction_df': DataFrame of predictions for retained classes.
            - 'discarded_classes_prediction_df': DataFrame of predictions for discarded classes.
            - 'retained_train_only_newName_classes_list': List of retained training-only class names.
    """
    # Get samples and training category data
    samples_and_train_categ_pred_data = get_samples_and_train_categ_counts_df(
        settings["training_classes_list"],
        prediction_implementation_details_df,
        settings["training_set_type"],
        settings["learning_objs_limit"],
    )

    # Get discarded classes
    discarded_samples_ecotaxa_class_list = samples_and_train_categ_pred_data.loc[
        ~samples_and_train_categ_pred_data["learned_category"], "ecotaxa_category"
    ].to_list()

    # Split prediction dataframe into retained and discarded classes
    retained_classes_prediction_df = prediction_df.loc[
        ~prediction_df["object_ecotaxa_annotation_category"].isin(
            discarded_samples_ecotaxa_class_list
        ),
        :,
    ]

    discarded_classes_prediction_df = prediction_df.loc[
        prediction_df["object_ecotaxa_annotation_category"].isin(
            discarded_samples_ecotaxa_class_list
        ),
        :,
    ]

    # Get retained training-only classes
    retained_train_only_ecotaxa_classes_list = samples_and_train_categ_pred_data.loc[
        samples_and_train_categ_pred_data["Category type"] == "Training only",
        "ecotaxa_category",
    ].to_list()

    # Sort and organize retained samples and training categories
    retained_sples_and_train_categ_pred_data = sort_and_move_to_end_df_rows(
        samples_and_train_categ_pred_data.loc[
            samples_and_train_categ_pred_data["learned_category"]
        ],
        "selected_samples_classes_count",
        "ecotaxa_category",
        retained_train_only_ecotaxa_classes_list,
    )

    # Get retained training-only newName classes
    retained_train_only_newName_classes_list = samples_and_train_categ_pred_data.loc[
        samples_and_train_categ_pred_data["ecotaxa_category"].isin(
            retained_train_only_ecotaxa_classes_list
        ),
        "newName",
    ].to_list()

    return {
        "samples_and_train_categ_pred_data": samples_and_train_categ_pred_data,
        "retained_sples_and_train_categ_pred_data": retained_sples_and_train_categ_pred_data,
        "retained_classes_prediction_df": retained_classes_prediction_df,
        "discarded_classes_prediction_df": discarded_classes_prediction_df,
        "retained_train_only_newName_classes_list": retained_train_only_newName_classes_list,
    }


def compute_classification_metrics(
    retained_classes_df: pd.DataFrame, retained_sples_data: pd.DataFrame, train_only_classes_list: list[str]
) -> tuple[np.ndarray, dict[str, dict[str, float]]]:
    """
    Compute the confusion matrix and classification report for retained classes.

    Parameters
    ----------
    retained_classes_df : pandas.DataFrame
        DataFrame containing predictions for retained classes (those present in both training and test sets).
    retained_sples_data : pandas.DataFrame
        DataFrame with metadata about retained sample categories, including class labels.
    train_only_classes_list : list
        List of class names that are present only in the training set.

    Returns
    -------
    confusion_matrix_array : numpy.ndarray
        Confusion matrix array (normalized by true labels).
    classif_report : dict
        Dictionary containing the classification report (precision, recall, f1-score, support)
        for each class and for macro/weighted averages. May be corrected for extra training-only classes.
    """
    # Create confusion matrix
    confusion_matrix_array = confusion_matrix(
        retained_classes_df["object_newname"],
        retained_classes_df["predicted_newName"],
        labels=retained_sples_data["newName"],
        normalize="true",
    )

    # Generate classification report
    classif_report = classification_report(
        retained_classes_df["object_newname"],
        retained_classes_df["predicted_newName"],
        labels=retained_sples_data["newName"],
        output_dict=True,
        zero_division=0,
    )

    # Remove unnecessary metrics
    classif_report.pop("accuracy", None)
    classif_report.pop("micro avg", None)

    # Correct classification report if necessary
    if train_only_classes_list:
        classif_report = correct_classif_rept_macro_avg_for_xtra_class(classif_report)

    return confusion_matrix_array, classif_report


def process_export_job(
    export_job_id: int,
    prediction_df: Optional[pd.DataFrame],
    exported_jobs_settings_dict: dict[int, Any],
    selected_samples_short_ID_mapping: pd.DataFrame,
) -> Optional[dict[str, Any]]:
    """
    Process a single export job, merging settings and prediction data, and compute all relevant metrics.

    Parameters
    ----------
    export_job_id : int
        The ID of the export job to process.
    prediction_df : pandas.DataFrame or None
        DataFrame containing prediction results for the export job, or None if not available.
    exported_jobs_settings_dict : dict
        Dictionary of settings for all export jobs.
    selected_samples_short_ID_mapping : pandas.DataFrame
        DataFrame mapping sample IDs to short sample labels.

    Returns
    -------
    result : dict or None
        Dictionary containing processed metrics and DataFrames for the export job, or None if
        prediction_df is None. Keys include region, use_scn_features, learning_objs_limit,
        strategy_numb, strategy_label, samples_and_train_categ_pred_data,
        retained_sples_and_train_categ_pred_data, prediction_confusion_matrix_array,
        classif_report, retained_classes_prediction_df, discarded_classes_prediction_df,
        and retained_train_only_newName_classes_list.
    """
    if prediction_df is None:
        return None

    # Merge with sample short ID
    prediction_df = merge_sample_short_ID_to_exported_df(
        selected_samples_short_ID_mapping, prediction_df
    )

    # Extract all job settings
    settings = extract_job_settings(export_job_id, exported_jobs_settings_dict)

    # Parse implementation details
    prediction_implementation_details_df = pd.DataFrame.from_dict(
        json.loads(settings["strategy_df"]),
        orient="columns",
    )

    # Process prediction data
    prediction_data = process_prediction_data(
        prediction_df, settings, prediction_implementation_details_df
    )

    # Compute classification metrics
    confusion_matrix_array, classif_report = compute_classification_metrics(
        prediction_data["retained_classes_prediction_df"],
        prediction_data["retained_sples_and_train_categ_pred_data"],
        prediction_data["retained_train_only_newName_classes_list"],
    )

    # Combine all results
    return {
        "region": settings["region"],
        "use_scn_features": settings["use_scn_features"],
        "learning_objs_limit": settings["learning_objs_limit"],
        "strategy_numb": settings["strategy_numb"],
        "strategy_label": settings["strategy_label"],
        "samples_and_train_categ_pred_data": prediction_data["samples_and_train_categ_pred_data"],
        "retained_sples_and_train_categ_pred_data": prediction_data[
            "retained_sples_and_train_categ_pred_data"
        ],
        "prediction_confusion_matrix_array": confusion_matrix_array,
        "classif_report": classif_report,
        "retained_classes_prediction_df": prediction_data["retained_classes_prediction_df"],
        "discarded_classes_prediction_df": prediction_data["discarded_classes_prediction_df"],
        "retained_train_only_newName_classes_list": prediction_data[
            "retained_train_only_newName_classes_list"
        ],
    }


def process_all_export_jobs(
    exported_jobs_df_dict: dict[int, Optional[pd.DataFrame]],
    exported_jobs_settings_dict: dict[int, Any],
    selected_samples_short_ID_mapping: pd.DataFrame,
    output_path: Union[str, Path],
) -> dict[int, dict[str, Any]]:
    """
    Process all export jobs and save the results to a pickle file.

    Parameters
    ----------
    exported_jobs_df_dict : dict
        Dictionary mapping export job IDs to their corresponding prediction DataFrames.
    exported_jobs_settings_dict : dict
        Dictionary mapping export job IDs to their settings.
    selected_samples_short_ID_mapping : pandas.DataFrame
        DataFrame mapping sample IDs to short sample labels.
    output_path : str or Path
        Path to the output pickle file where the processed results will be saved.

    Returns
    -------
    exported_jobs_data_and_metrics_dict : dict
        Dictionary mapping export job IDs to their processed metrics and DataFrames.
        The same dictionary is also saved to the specified pickle file.
    """
    exported_jobs_data_and_metrics_dict: dict[int, dict[str, Any]] = {}

    for export_job_id, prediction_df in list(exported_jobs_df_dict.items()):
        job_results = process_export_job(
            export_job_id,
            prediction_df,
            exported_jobs_settings_dict,
            selected_samples_short_ID_mapping,
        )

        if job_results:
            exported_jobs_data_and_metrics_dict[export_job_id] = job_results

    with open(output_path, "wb") as f:
        pickle.dump(exported_jobs_data_and_metrics_dict, f)

    return exported_jobs_data_and_metrics_dict


def extract_job_data(export_job_id: int, metrics_dict: dict, settings_dict: dict) -> dict:
    """
    Extract all job-related data from the metrics dictionary for a specific export job.

    Parameters
    ----------
    export_job_id : int
        ID of the export job to extract data for
    metrics_dict : dict
        Dictionary containing metrics data for all jobs
    settings_dict : dict
        Dictionary containing settings data for all jobs

    Returns
    -------
    dict
        Dictionary containing all extracted job data
    """
    job_data = metrics_dict[export_job_id]

    # Create configuration dict for classification figure
    conf_mtrx_fig_params = {
        "figure_title": settings_dict[export_job_id][1]["scenario_title"],
        "learning_objs_limit": job_data["learning_objs_limit"],
        "use_scn_features": job_data["use_scn_features"],
        "strategy_numb": job_data["strategy_numb"],
        "extra_training_only_classes": job_data["retained_train_only_newName_classes_list"],
    }

    return {
        "region": job_data["region"],
        "use_scn_features": job_data["use_scn_features"],
        "learning_objs_limit": job_data["learning_objs_limit"],
        "strategy_numb": job_data["strategy_numb"],
        "strategy_label": job_data["strategy_label"],
        "samples_data": job_data["samples_and_train_categ_pred_data"],
        "retained_samples_data": job_data["retained_sples_and_train_categ_pred_data"],
        "confusion_matrix": job_data["prediction_confusion_matrix_array"],
        "classif_report": job_data["classif_report"],
        "retained_classes_df": job_data["retained_classes_prediction_df"],
        "discarded_classes_df": job_data["discarded_classes_prediction_df"],
        "retained_train_only_classes": job_data["retained_train_only_newName_classes_list"],
        "conf_mtrx_fig_params": conf_mtrx_fig_params,
    }


def process_general_category_redistribution(
    prediction_df: pd.DataFrame, general_classes_set: set
) -> bool:
    """
    Process and redistribute general categories in the prediction dataframe.

    Parameters
    ----------
    prediction_df : pandas.DataFrame
        DataFrame containing prediction data
    general_classes_set : set
        Set of general class names to be redistributed

    Returns
    -------
    bool
        True if both original and predicted classes were redistributed, False otherwise
    """
    original_redistributed = False
    predicted_redistributed = False

    # Check original classes
    original_classes = set(prediction_df["object_newname"].unique())
    present_orig_general_classes = original_classes.intersection(general_classes_set)

    if present_orig_general_classes:
        redistribute_general_categories(
            prediction_df, present_orig_general_classes, "object_newname"
        )
        original_redistributed = True

    # Check predicted classes
    predicted_classes = set(prediction_df["predicted_newName"].unique())
    present_pred_general_classes = predicted_classes.intersection(general_classes_set)

    if present_pred_general_classes:
        redistribute_general_categories(
            prediction_df, present_pred_general_classes, "predicted_newName"
        )
        predicted_redistributed = True

    # Only return True if BOTH were redistributed
    return original_redistributed and predicted_redistributed


def save_figures_to_pdf(figures: list[Optional[Figure]], filepath: Union[str, Path]) -> None:
    """
    Save a list of figures to a PDF file.

    Parameters
    ----------
    figures : list[matplotlib.figure.Figure | None]
        List of figures to save (can contain None values)
    filepath : str or pathlib.Path
        Path to the output PDF file
    """
    with PdfPages(filepath) as pdf_file:
        for fig in figures:
            if fig is not None:  # Only save non-None figures
                pdf_file.savefig(fig)


def create_pdf_filename(job_data: dict, export_job_id: int) -> str:
    """
    Create standardized PDF filename based on job data.

    Parameters
    ----------
    job_data : dict
        Dictionary containing job data
    export_job_id : int
        ID of the export job

    Returns
    -------
    str
        Generated filename
    """
    region = job_data["region"]
    strategy_numb = job_data["strategy_numb"]
    learning_objs_limit = job_data["learning_objs_limit"]
    use_scn_features = job_data["use_scn_features"]

    return (
        f"{region}-strat{strategy_numb}-"
        f"{learning_objs_limit if learning_objs_limit < 999999 else 'Max'}Obj"
        f"{'-SCN' if use_scn_features else ''}"
        f"-expID{export_job_id}-classif report.pdf"
    )


def process_export_job_visualizations(
    export_job_id: int,
    prediction_df: Optional[pd.DataFrame],
    jobs_data_dict: dict,
    jobs_settings_dict: dict,
    sample_id_mapping: pd.DataFrame,
    taxo_category_mapping_df: pd.DataFrame,
    results_root_path: Path,
    general_classes_set: set,
    show_figure: bool = False,
) -> bool:
    """
    Process a single export job to generate and save visualization figures.

    Parameters
    ----------
    export_job_id : int
        ID of the export job to process
    prediction_df : pandas.DataFrame or None
        DataFrame containing prediction data, or None if not available
    jobs_data_dict : dict
        Dictionary containing metrics data for all jobs
    jobs_settings_dict : dict
        Dictionary containing settings data for all jobs
    sample_id_mapping : pandas.DataFrame
        DataFrame mapping sample IDs to short sample IDs
    taxo_category_mapping_df : pandas.DataFrame
        DataFrame mapping taxonomy categories
    results_root_path : pathlib.Path
        Root path for saving results
    general_classes_set : set
        Set of general class names to be redistributed
    show_figure : bool, default=False
        Whether to display figures interactively

    Returns
    -------
    bool
        True if processing was successful, False otherwise
    """
    if prediction_df is None:
        return False

    try:
        # Merge prediction data with sample ID mapping
        prediction_df = merge_sample_short_ID_to_exported_df(sample_id_mapping, prediction_df)

        # Extract job data
        job_data = extract_job_data(export_job_id, jobs_data_dict, jobs_settings_dict)

        # Create figures list to collect all generated figures
        figures: list[Optional[Figure]] = []

        # Create classification figure
        classification_fig = create_classification_figure(job_data)
        figures.append(classification_fig)
        if show_figure:
            plt.show()

        # Create discarded categories heatmap
        heatmap_fig = create_discarded_categories_heatmap(job_data, taxo_category_mapping_df)
        if heatmap_fig is not None:
            figures.append(heatmap_fig)
            if show_figure:
                plt.show()

        # Create relative abundance figure
        rel_abund_fig = create_relative_abundance_figure(prediction_df)
        if rel_abund_fig is not None:
            figures.append(rel_abund_fig)
            if show_figure:
                plt.show()

        # Process general category redistribution
        redistribution_performed = process_general_category_redistribution(
            prediction_df, general_classes_set
        )

        # Create redistributed relative abundance figure if needed
        if redistribution_performed:
            redist_rel_abund_fig = create_relative_abundance_figure(
                prediction_df, redistributed=True
            )
            if redist_rel_abund_fig is not None:  # Only add to figures list if not None
                figures.append(redist_rel_abund_fig)
                if show_figure:
                    plt.show()

        # Create output directory
        region = job_data["region"]
        region_figs_path = (
            results_root_path / region / "figures" / "Classification_strategy_figures"
        )
        region_figs_path.mkdir(parents=True, exist_ok=True)

        # Create output filename and save figures
        pdf_filename = create_pdf_filename(job_data, export_job_id)
        pdf_filepath = region_figs_path / pdf_filename
        save_figures_to_pdf(figures, pdf_filepath)

        return True

    except Exception as e:
        print(f"Error processing export job {export_job_id}: {str(e)}")
        return False


def process_all_visualizations(
    exported_jobs_df_dict: dict[int, Optional[pd.DataFrame]],
    exported_jobs_data_and_metrics_dict: dict,
    exported_jobs_settings_dict: dict,
    selected_samples_short_ID_mapping: pd.DataFrame,
    taxo_category_mapping_df: pd.DataFrame,
    results_root_path: Path,
    general_classes_set: set,
    show_figure: bool = False,
) -> int:
    """
    Process all export jobs to generate and save visualization figures.

    Parameters
    ----------
    exported_jobs_df_dict : dict
        Dictionary mapping export job IDs to prediction DataFrames
    exported_jobs_data_and_metrics_dict : dict
        Dictionary containing metrics data for all jobs
    exported_jobs_settings_dict : dict
        Dictionary containing settings data for all jobs
    selected_samples_short_ID_mapping : pandas.DataFrame
        DataFrame mapping sample IDs to short sample IDs
    taxo_category_mapping_df : pandas.DataFrame
        DataFrame mapping taxonomy categories
    results_root_path : pathlib.Path
        Root path for saving results
    general_classes_set : set
        Set of general class names to be redistributed
    show_figure : bool, default=False
        Whether to display figures interactively

    Returns
    -------
    int
        Number of successfully processed jobs
    """
    successful_jobs = 0

    for export_job_id, prediction_df in list(exported_jobs_df_dict.items()):
        success = process_export_job_visualizations(
            export_job_id,
            prediction_df,
            exported_jobs_data_and_metrics_dict,
            exported_jobs_settings_dict,
            selected_samples_short_ID_mapping,
            taxo_category_mapping_df,
            results_root_path,
            general_classes_set,
            show_figure,
        )
        if success:
            successful_jobs += 1

    # Clean up by closing all figures to free memory
    plt.close("all")

    return successful_jobs


def create_classifier_metrics_dataframe(exported_jobs_data_and_metrics_dict: dict) -> pd.DataFrame:
    """
    Create a DataFrame summarizing classifier metrics and strategy parameters.

    Parameters
    ----------
    exported_jobs_data_and_metrics_dict : dict
        Dictionary mapping export job IDs to job data. Each value should contain:
            - "region": str
            - "use_scn_features": bool
            - "learning_objs_limit": int
            - "strategy_numb": int
            - "strategy_label": str
            - "classif_report": dict with "macro avg" and "weighted avg" metrics

    Returns
    -------
    pd.DataFrame
        DataFrame with multi-level columns for strategy parameters and classification metrics.
    """
    strategy_multi_idx = pd.MultiIndex.from_product(
        [
            ["Strategy"],
            [
                "region",
                "use_scn_features",
                "learning_objs_limit",
                "strategy_numb",
                "strategy_label",
            ],
        ]
    )
    multi_idx = pd.MultiIndex.from_product(
        [["macro avg", "weighted avg"], ["precision", "recall", "f1-score"]],
        names=["avg_type", "score_type"],
    )
    classifier_strategies_metrics_df = pd.concat(
        [pd.DataFrame(columns=strategy_multi_idx), pd.DataFrame(columns=multi_idx)]
    )

    for i, (export_job_id, data_dict) in enumerate(exported_jobs_data_and_metrics_dict.items()):
        classifier_strategies_metrics_df.loc[i, ("Strategy", "region")] = data_dict["region"]
        classifier_strategies_metrics_df.loc[i, ("Strategy", "use_scn_features")] = data_dict[
            "use_scn_features"
        ]
        classifier_strategies_metrics_df.loc[i, ("Strategy", "learning_objs_limit")] = data_dict[
            "learning_objs_limit"
        ]
        classifier_strategies_metrics_df.loc[i, ("Strategy", "strategy_numb")] = data_dict[
            "strategy_numb"
        ]
        classifier_strategies_metrics_df.loc[i, ("Strategy", "strategy_label")] = data_dict[
            "strategy_label"
        ]

        for avg_type, score_dict in list(data_dict["classif_report"].items())[-2:]:
            for key, item in score_dict.items():
                if key == "support":
                    classifier_strategies_metrics_df.loc[i, ("", "support")] = item
                else:
                    classifier_strategies_metrics_df.loc[
                        i, (avg_type.replace(" (corr)", ""), key)
                    ] = item

    return classifier_strategies_metrics_df
