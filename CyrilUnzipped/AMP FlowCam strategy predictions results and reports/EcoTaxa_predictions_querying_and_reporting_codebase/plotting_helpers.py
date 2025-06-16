from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd  # type: ignore
import seaborn as sns  # type: ignore
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text
from matplotlib.backends.backend_pdf import PdfPages

tab20_colormap = plt.cm.get_cmap("tab20", 20)  # 20 distinct colors
tab20b_colormap = plt.cm.get_cmap("tab20b", 20)  # 20 distinct colors

# Create a dictionary mapping from the 'tab20' colormap
tab20_colors = {f"color_{i}": tab20_colormap(i) for i in range(tab20_colormap.N)}
tab20b_colors = {f"color_{i}b": tab20b_colormap(i) for i in range(tab20b_colormap.N)}
tab40_colors = {**tab20_colors, **tab20b_colors}


def compute_dual_top_class_relat_abund_dfs(
    df: pd.DataFrame,
    taxonomic_col1: str,
    taxonomic_col2: str,
    sample_col: str,
    number_of_largest: int = 11,
) -> tuple[list[str], pd.DataFrame, pd.DataFrame]:
    """
    Computes relative abundance DataFrames for the top classes in two taxonomic columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing sample and taxonomic information.
    taxonomic_col1 : str
        Name of the first taxonomic column.
    taxonomic_col2 : str
        Name of the second taxonomic column.
    sample_col : str
        Name of the sample identifier column.
    number_of_largest : int, optional
        Number of top classes to consider (default is 11).

    Returns
    -------
    top_combined_classes : list of str
        List of combined top classes from both taxonomic columns.
    top_classes_relative_abund_df1 : pandas.DataFrame
        Relative abundance DataFrame for taxonomic_col1.
    top_classes_relative_abund_df2 : pandas.DataFrame
        Relative abundance DataFrame for taxonomic_col2.
    """
    count_df1 = df.groupby([sample_col, taxonomic_col1]).size().unstack(fill_value=0)
    top_classes_df1 = count_df1.sum(axis=0).nlargest(number_of_largest)

    count_df2 = df.groupby([sample_col, taxonomic_col2]).size().unstack(fill_value=0)
    top_classes_df2 = count_df2.sum(axis=0).nlargest(number_of_largest)

    top_combined_classes = get_combined_top_classes(
        number_of_largest, top_classes_df1, top_classes_df2
    )

    # Group all other classes under "others"
    count_top_classes_df1 = reshape_df_top_classes_and_others(count_df1, top_combined_classes)
    count_top_classes_df2 = reshape_df_top_classes_and_others(count_df2, top_combined_classes)

    # Order alphabetically the classes columns
    count_top_classes_df1 = count_top_classes_df1.sort_index(axis=1)
    count_top_classes_df2 = count_top_classes_df2.sort_index(axis=1)

    # Convert counts to relative abundance
    top_classes_relative_abund_df1 = get_top_classes_relative_abund_df(count_top_classes_df1)
    top_classes_relative_abund_df2 = get_top_classes_relative_abund_df(count_top_classes_df2)
    return (
        top_combined_classes,
        top_classes_relative_abund_df1,
        top_classes_relative_abund_df2,
    )


def get_predicted_taxonomic_redistribution_df(
    taxo_category_mapping_df: pd.DataFrame,
    discarded_classes_prediction_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Computes a contingency table of predicted taxonomic redistributions and adds percentage and annotation columns.

    Parameters
    ----------
    taxo_category_mapping_df : pandas.DataFrame
        DataFrame containing taxonomic category mappings.
    discarded_classes_prediction_df : pandas.DataFrame
        DataFrame with predictions for discarded classes.

    Returns
    -------
    predicted_taxonomic_redistribution_df : pandas.DataFrame
        DataFrame with contingency counts, percentages, annotation strings, and mappings.
    """
    predicted_taxonomic_redistribution_df = pd.crosstab(
        discarded_classes_prediction_df["object_newname"],
        discarded_classes_prediction_df["predicted_newName"],
    )
    predicted_taxonomic_redistribution_df = predicted_taxonomic_redistribution_df.T
    predicted_taxonomic_redistribution_df.reset_index(inplace=True)
    predicted_taxonomic_redistribution_df.columns.name = None
    predicted_taxonomic_redistribution_df.set_index("predicted_newName", inplace=True)
    total_counts = predicted_taxonomic_redistribution_df.sum()
    predicted_taxonomic_redistribution_df.loc["total_counts"] = total_counts
    for col in predicted_taxonomic_redistribution_df.columns:
        predicted_taxonomic_redistribution_df.loc[:, f"{col} pct"] = (
            predicted_taxonomic_redistribution_df.loc[:, col].divide(
                predicted_taxonomic_redistribution_df.loc["total_counts", col]
            )
            * 100
        )
        percentages_str = predicted_taxonomic_redistribution_df.loc[:, f"{col} pct"].map(
            lambda x: f"{x: .1f}" if x > 0 else f"{x: .1f}"
        )
        contingency_str = predicted_taxonomic_redistribution_df.loc[:, col].map(
            lambda x: f"(n={x})" if x > 0 else ""
        )
        predicted_taxonomic_redistribution_df.loc[:, f"{col} annotations"] = (
            percentages_str + r"\% " + "\\small" + contingency_str
        )
    predicted_taxonomic_redistribution_df = predicted_taxonomic_redistribution_df.merge(
        taxo_category_mapping_df.set_index("newName").loc[:, "Copepoda_mapping"],
        how="left",
        left_index=True,
        right_index=True,
    )

    return predicted_taxonomic_redistribution_df


def get_top_classes_relative_abund_df(count_top_classes_df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts a count DataFrame of top classes to relative abundance per row.

    Parameters
    ----------
    count_top_classes_df : pandas.DataFrame
        DataFrame with counts of top classes (columns) per sample (rows).

    Returns
    -------
    top_classes_relative_abundance_df : pandas.DataFrame
        DataFrame with relative abundances (fractions) for each class per sample.
    """
    top_classes_relative_abundance_df = count_top_classes_df.divide(
        count_top_classes_df.sum(axis=1), axis=0
    )
    relative_abundance1_col_idx = top_classes_relative_abundance_df.columns.to_list()
    relative_abundance1_col_idx.remove("Others")
    relative_abundance1_col_idx.sort(reverse=True)
    top_classes_relative_abundance_df = top_classes_relative_abundance_df[
        ["Others"] + relative_abundance1_col_idx
    ]
    return top_classes_relative_abundance_df


def reshape_df_top_classes_and_others(
    count_df: pd.DataFrame, top_combined_classes: list[str]
) -> pd.DataFrame:
    """
    Reshapes a count DataFrame to include only top classes and aggregates all others into an 'Others' column.

    Parameters
    ----------
    count_df : pandas.DataFrame
        DataFrame with class counts per sample.
    top_combined_classes : list of str
        List of top class names to retain.

    Returns
    -------
    reshaped_df : pandas.DataFrame
        DataFrame with columns for each top class and an 'Others' column.
    """
    other_columns = [col for col in count_df.columns if col not in top_combined_classes]
    count_df["Others"] = count_df.loc[:, other_columns].sum(axis=1)
    df_top_classes = [
        categ for categ in top_combined_classes if categ in count_df.columns.to_list()
    ]
    count_df = count_df[df_top_classes + ["Others"]]
    return count_df


def get_combined_top_classes(
    number_of_largest: int,
    top_classes_df1: pd.Series,
    top_classes_df2: pd.Series,
) -> list[str]:
    """
    Combines the top classes from two DataFrames and returns the most frequent ones.

    Parameters
    ----------
    number_of_largest : int
        Number of top classes to return.
    top_classes_df1 : pandas.Series
        Series of top classes from the first DataFrame.
    top_classes_df2 : pandas.Series
        Series of top classes from the second DataFrame.

    Returns
    -------
    top_combined_classes : list of str
        List of the most frequent combined top classes.
    """
    combined_count_df = pd.concat([top_classes_df1, top_classes_df2], axis=0)

    combined_count_df.sort_values()
    top_combined_classes = combined_count_df.index.drop_duplicates().to_list()[:number_of_largest]
    top_combined_classes = (
        top_combined_classes[:number_of_largest]
        if len(top_combined_classes) > number_of_largest
        else top_combined_classes
    )

    return top_combined_classes


def round_percent_number(value: float) -> str:
    """
    Converts a float value to a percentage string for LaTeX formatting.

    Parameters
    ----------
    value : float
        The value to convert, expected to be a fraction (e.g., 0.05 for 5%).

    Returns
    -------
    percent_str : str
        String representation of the percentage (e.g., '5\\%').
        Returns '<1\\%' if value is between 0 and 0.01, and an empty string if value is 0 or negative.
    """
    if value >= 0.01:
        value_str = f"{value * 100:.0f}"
        return value_str + r"\%"
    elif 0.01 > value > 0.0:
        return r"<1\%"
    else:
        return ""


def get_df_index_in_str_list(df: pd.DataFrame, str_list: list[str]) -> set[str]:
    """
    Finds the intersection between a DataFrame's index and a list of strings.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame whose index will be compared.
    str_list : list of str
        List of strings to compare against the DataFrame index.

    Returns
    -------
    common_elements : set of str
        Set of elements present in both the DataFrame index and the string list.
    """
    columns_set = set(df.index)
    strings_set = set(str_list)

    # Find the intersection
    common_elements = columns_set.intersection(strings_set)
    return common_elements


def built_classif_figure_by_process(
    conf_matrix: np.ndarray,
    classif_report: dict[str, dict[str, float]],
    fig_params_dict: dict,
    predicted_x_axis_labels: Union[list[str], pd.Series],
    actual_y_axis_labels: Union[list[str], pd.Series],
    frequency_y_report_labels: list[str],
) -> Figure:
    """
    Build a classification figure with confusion matrix and classification report.

    Parameters
    ----------
    conf_matrix : np.ndarray
        Confusion matrix as a 2D numpy array.
    classif_report : dict of str to dict of str to float
        Classification report as a nested dictionary.
    fig_params_dict : dict
        Dictionary of figure parameters (e.g., title, extra classes, etc.).
    predicted_x_axis_labels : list of str or pd.Series
        Labels for the predicted values (x-axis).
    actual_y_axis_labels : list of str or pd.Series
        Labels for the actual values (y-axis).
    frequency_y_report_labels : list of str
        Frequency labels for the y-axis of the report.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The resulting matplotlib Figure object.
    """
    fig = plt.figure(figsize=(27, 20))
    gs = GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[20, 1])

    ax0 = fig.add_subplot(gs[0, 0])

    ax1 = fig.add_subplot(gs[:, 1])

    fig.suptitle(
        fig_params_dict["figure_title"],
        fontsize=26,
        fontweight="bold",
    )
    plot_classif_conf_matrix(
        conf_matrix,
        ax0,
        predicted_x_axis_labels,
        actual_y_axis_labels,
        fig_params_dict["extra_training_only_classes"],
    )

    plot_classif_report_matrix(
        classif_report,
        fig_params_dict["learning_objs_limit"],
        ax1,
        frequency_y_report_labels,
        fig_params_dict["extra_training_only_classes"],
    )

    scn_objs_strat_numb_text_box = (
        f"Use of SCN features: {'Yes' if fig_params_dict['use_scn_features'] else 'No'}\n"
        f"Max learning objects: "
        f"{'Maximum' if fig_params_dict['learning_objs_limit'] == 999999 else fig_params_dict['learning_objs_limit']} "
        f"objects/class\nStrategy N° {fig_params_dict['strategy_numb']}"
    )

    fig.text(
        0.02,
        0.92,
        scn_objs_strat_numb_text_box,
        fontsize=21,
        bbox=dict(facecolor="gray", alpha=0.2),
    )

    return fig


def plot_classif_conf_matrix(
    conf_matrix: np.ndarray,
    axes: Axes,
    predicted_x_axis_labels: Union[list[str], pd.Series],
    actual_y_axis_labels: Union[list[str], pd.Series],
    extra_training_only_classes: Union[list[str], pd.Series],
) -> None:
    """
    Plot a confusion matrix heatmap with special highlighting for extra training-only classes.

    Parameters
    ----------
    conf_matrix : np.ndarray
        Confusion matrix as a 2D numpy array.
    axes : matplotlib.axes.Axes
        Axes object to plot on.
    predicted_x_axis_labels : list of str or pd.Series
        Labels for the predicted values (x-axis).
    actual_y_axis_labels : list of str or pd.Series
        Labels for the actual values (y-axis).
    extra_training_only_classes : list of str or pd.Series
        List of extra classes only present in training.

    Returns
    -------
    None
    """
    if extra_training_only_classes:
        conf_matrix = conf_matrix[: -len(extra_training_only_classes)]
        actual_y_axis_labels = actual_y_axis_labels[: -len(extra_training_only_classes)]
        mask = np.full(conf_matrix.shape, True)
        mask[:, -len(extra_training_only_classes) :] = False  # For column
    else:
        mask = np.full(conf_matrix.shape, True)

    annot_labels = pd.DataFrame(conf_matrix)
    annot_labels = annot_labels.map(round_percent_number)  # type: ignore
    conf_matrix_axes: plt.Axes = sns.heatmap(
        conf_matrix,
        mask=~mask,
        ax=axes,
        annot=annot_labels,
        fmt="",
        cmap="GnBu",
        cbar=False,
        linewidths=0.5,
        linecolor="#ced5b4",
        annot_kws={"fontsize": 14},
    )
    if extra_training_only_classes:
        sns.heatmap(
            conf_matrix,
            mask=mask,
            ax=axes,
            annot=annot_labels,
            fmt="",
            cmap="Reds",
            cbar=False,
            linewidths=0.5,
            linecolor="#ced5b4",
            annot_kws={"fontsize": 14},
        )
    conf_matrix_axes.set_title(
        "Confusion Matrix - In percent of Actual Value",
        fontdict={"fontsize": 20, "fontweight": "bold"},
        y=1.025,
    )
    conf_matrix_axes.set_xlabel(
        "\n Predicted Values\n", fontdict={"fontsize": 18, "fontweight": "bold"}
    )
    conf_matrix_axes.set_ylabel("Actual Values\n", fontdict={"fontsize": 18, "fontweight": "bold"})
    conf_matrix_axes.xaxis.set_ticklabels(
        predicted_x_axis_labels,
        fontsize=15,
        rotation=-45,
        ha="left",
        rotation_mode="anchor",
    )
    conf_matrix_axes.yaxis.set_ticklabels(
        actual_y_axis_labels,
        fontsize=15,
        rotation=0,
    )
    for _, spine in conf_matrix_axes.spines.items():
        spine.set_visible(True)
        spine.set_color("white")

    if extra_training_only_classes:

        n_rows, n_cols = conf_matrix.shape

        # Highlight the last three rows and last three columns (indexing starts from 0)
        # highlight_rows_start = n_rows - len(extra_training_only_classes)
        highlight_cols_start = n_cols - len(extra_training_only_classes)

        # Define soft colors using hexadecimal color codes
        soft_red = "#fbb4ae"  # Hex for a light shade of red

        # Add patch for the column
        rect_col = patches.Rectangle(
            (highlight_cols_start, 0),
            width=n_rows,
            height=conf_matrix.shape[0],
            linewidth=2,
            edgecolor=soft_red,
            facecolor="none",
        )

        conf_matrix_axes.add_patch(rect_col)

        # Create custom legend items
        red_patch = Patch(color="red", label="Extra\ntraining\nclasses")
        conf_matrix_axes.legend(
            handles=[red_patch],
            # loc="upper left",
            bbox_to_anchor=(0.95, -0.2, 0.1, 1),
            bbox_transform=conf_matrix_axes.transAxes,
            loc="lower left",
            borderaxespad=0.0,
        )

        for label in conf_matrix_axes.get_xticklabels():
            if (
                label.get_text() in extra_training_only_classes
            ):  # Specify the label you want to change
                label.set_color("red")  # Set the color


def plot_classif_report_matrix(
    classif_report: dict[str, dict[str, float]],
    learning_objs_limit: int,
    axes: Axes,
    frequency_y_report_labels: list[str],
    extra_training_only_classes: Union[list[str], pd.Series],
    cbar=True,
    yshared=False,
) -> None:
    """
    Plot a classification report as a heatmap matrix.

    Parameters
    ----------
    classif_report : dict of str to dict of str to float
        Classification report as a nested dictionary.
    learning_objs_limit : int
        Maximum number of learning objects per class.
    axes : matplotlib.axes.Axes
        Axes object to plot on.
    frequency_y_report_labels : list of str
        Frequency labels for the y-axis.
    extra_training_only_classes : list of str or pd.Series
        List of extra classes only present in training.
    cbar : bool, optional
        Whether to display the colorbar (default is True).
    yshared : bool, optional
        Whether to share y-axis ticks (default is False).

    Returns
    -------
    None
    """
    if extra_training_only_classes:
        bottom_metrics_list = [
            "macro avg (corr)",
            "weighted avg",
        ]
    else:
        bottom_metrics_list = ["macro avg", "weighted avg"]

    len_bottom_metrics = len(bottom_metrics_list)

    df_report = pd.DataFrame(classif_report).transpose().drop("support", axis=1)
    classif_report_axes: plt.Axes = sns.heatmap(
        df_report,
        ax=axes,
        annot=True,
        fmt=".2f",
        cmap="GnBu",
        cbar=cbar,
        annot_kws={"fontsize": 14},
    )
    classif_report_axes.set_title(
        (
            f"Classification Report Matrix\nmax "
            f"{learning_objs_limit if learning_objs_limit < 999999 else 'available'} "
            f"learning objects per class"
        ),
        fontdict={"fontsize": 18, "fontweight": "bold"},
        y=1.025,
    )
    classif_report_axes.xaxis.set_ticklabels(
        classif_report_axes.get_xmajorticklabels(), fontsize=16
    )
    classif_report_axes.tick_params(
        axis="x", bottom=True, top=True, labelbottom=True, labeltop=True
    )
    classif_report_axes.yaxis.set_ticklabels(classif_report_axes.get_ymajorticklabels(), fontsize=7)
    tick_label: Text
    if classif_report_axes.axes:
        for i, tick_label in enumerate(classif_report_axes.axes.get_yticklabels()):
            tick_text: str = tick_label.get_text()
            if tick_text in bottom_metrics_list:
                tick_label.set_color("#d66f00")
                tick_label.set_fontweight("bold")
                tick_label.set_fontsize(14)
                for j in range(3 * i, 3 * i + 3):
                    annotation: Text = classif_report_axes.axes.texts[j]
                    annotation.set_fontweight("heavy")
                    annotation.set_color("#000000")
                    annotation.set_fontsize(18)
            elif tick_text in extra_training_only_classes:
                tick_label.set_color("red")
                for j in range(3 * i, 3 * i + 3):
                    annotation = classif_report_axes.axes.texts[j]
                    annotation.set_text("-")
                    annotation.set_fontweight("heavy")
                    annotation.set_color("#000000")
                    annotation.set_fontsize(18)
            elif tick_text == "accuracy":
                classif_report_axes.axes.texts[3 * i].set_text("")
                classif_report_axes.axes.texts[3 * i + 1].set_text("model's scores")
                classif_report_axes.axes.texts[3 * i + 1].set_fontsize(18)
                classif_report_axes.axes.texts[3 * i + 1].set_fontweight("heavy")
                classif_report_axes.axes.texts[3 * i + 1].set_horizontalalignment("center")
                classif_report_axes.axes.texts[3 * i + 2].set_text("")

        frequency_y_report_labels = [
            (
                "$\\fontfamily{iwona}\\selectfont\\textbf{"
                + lbl_str.split("\n")[0]
                + "}$\n$\\fontfamily{iwona}\\selectfont\\textrm{"
                + lbl_str.split("\n")[1]
                + "}$"
                if ("\n" in lbl_str)
                else ("$\\fontfamily{iwona}\\selectfont\\textbf{" + lbl_str + "}$")
            )
            for lbl_str in frequency_y_report_labels
        ]

        classif_report_axes_label = [
            ("$\\fontfamily{iwona}\\selectfont\\textbf{" + str(lbl_str.get_text()) + "}$")
            for lbl_str in (classif_report_axes.get_ymajorticklabels()[-len_bottom_metrics:])
        ]

        frequency_y_report_labels = (
            frequency_y_report_labels + classif_report_axes_label  # type: ignore
        )

        if yshared:
            classif_report_axes.axes.set_yticks([])

        classif_report_axes.yaxis.set_ticklabels(
            frequency_y_report_labels,
            fontsize=14,
        )

        for tick_label in classif_report_axes.get_ymajorticklabels()[-len_bottom_metrics:]:
            tick_label.set_fontweight("bold")
            tick_label.set_fontsize(18)


def plot_dual_relative_abundance(
    top_combined_classes: list[str],
    top_classes_relative_abundance1: pd.DataFrame,
    top_classes_relative_abundance2: pd.DataFrame,
    redistributed: bool,
) -> Figure:
    """
    Plot the relative abundance of top taxonomic classes for
    two datasets per sample using stacked bar chart.

    Parameters
    ----------
    top_combined_classes : list of str
        List of top combined class names.
    top_classes_relative_abundance1 : pandas.DataFrame
        Relative abundance DataFrame for the first dataset (e.g., validation).
    top_classes_relative_abundance2 : pandas.DataFrame
        Relative abundance DataFrame for the second dataset (e.g., prediction).
    redistributed : bool
        Whether the data has been redistributed.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The resulting matplotlib Figure object.
    """
    colormap = dict(zip(top_combined_classes + ["Others"], tab40_colors.values()))

    ax: Axes
    fig, ax = plt.subplots(figsize=(15, 10))
    _ = top_classes_relative_abundance2.plot(
        kind="bar",
        stacked=True,
        color=[colormap[x] for x in top_classes_relative_abundance2.columns],
        ax=ax,
        position=0,
        width=0.4,
        edgecolor="#d6d6d6",
        align="center",
    )

    bar2 = top_classes_relative_abundance1.plot(
        kind="bar",
        stacked=True,
        color=[colormap[x] for x in top_classes_relative_abundance1.columns],
        ax=ax,
        position=1,
        width=0.4,
        edgecolor="#d6d6d6",
        align="center",
    )

    for i, container in enumerate(bar2.containers):
        rect: Rectangle
        for rect in container:
            rect.set_x(rect.get_x() + 0.05)
            rect.set_width(0.3)

    xmin, xmax = ax.get_xlim()

    ax.set_xlim(xmin, xmax + 0.2)

    # Adding annotations 'Val' and 'Pred'
    for i in range(len(bar2.get_xticks())):
        ax.text(i - 0.2, 1, "Val", ha="center", va="bottom")
        ax.text(i + 0.2, 1, "Pred", ha="center", va="bottom")

    ax.set_ylabel("Relative Abundance")
    ax.set_xlabel("Sample Short ID")
    plt.title(
        (
            f"Relative Abundance of Top Taxonomic "
            f"Instances per Sample{' (Redistributed)' if redistributed else ''}"
        )
    )
    plt.xticks(rotation=0)

    # Create a single legend for both plots
    handles, labels = ax.get_legend_handles_labels()
    labels.reverse()
    handles.reverse()
    by_label = dict(zip(labels, handles))
    plt.legend(
        by_label.values(),
        by_label.keys(),
        title="Most abundant classes (In Val or Pred)",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    return fig


def plot_discarded_categ_redistributed_predictions_heatmap(
    predicted_taxonomic_redistribution_df: pd.DataFrame,
    discarded_sples_and_train_categ_pred_data: pd.DataFrame,
    train_kept_train_only_newName_classes_list: list[str] = [],
) -> Figure:
    """
    Plot a heatmap showing the predicted redistribution of discarded taxonomic categories.

    Parameters
    ----------
    predicted_taxonomic_redistribution_df : pandas.DataFrame
        DataFrame with contingency counts, percentages, and annotations for predicted redistributions.
    discarded_sples_and_train_categ_pred_data : pandas.DataFrame
        DataFrame with information about discarded samples and their categories.
    train_kept_train_only_newName_classes_list : list of str, optional
        List of extra training-only class names to highlight (default is empty list).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The resulting matplotlib Figure object.
    """
    discarded_newname_categories_list = discarded_sples_and_train_categ_pred_data.loc[
        :, "newName"
    ].to_list()

    original_cat_count_label = discarded_sples_and_train_categ_pred_data.loc[
        :, "category_frequency_labels"
    ].to_list()

    fig_width = (
        10
        if len(discarded_newname_categories_list) < 7
        else (10 + len(discarded_newname_categories_list) - 6)
    )
    fig = plt.figure(figsize=(fig_width, 8))
    fig.subplots_adjust(top=0.9)
    ax: plt.Axes = sns.heatmap(
        predicted_taxonomic_redistribution_df.drop("total_counts").loc[
            :, [f"{categ_name} pct" for categ_name in discarded_newname_categories_list]
        ],
        annot=predicted_taxonomic_redistribution_df.drop("total_counts").loc[
            :, [f"{categ_name} annotations" for categ_name in discarded_newname_categories_list]
        ],
        fmt="",
        cmap="YlGnBu",
        cbar=True,
        xticklabels=original_cat_count_label,
    )

    plt.title("Predictions of discarded taxa from training", pad=15)
    plt.xlabel("Actual discarded Taxa")
    plt.ylabel("Predicted Taxa")
    plt.xticks(
        rotation=-25,
        ha="left",
        rotation_mode="anchor",
    )
    for label in ax.get_yticklabels():
        if (
            label.get_text() in train_kept_train_only_newName_classes_list
        ):  # Specify the label you want to change
            label.set_color("red")  # Set the color
    for label in ax.get_xticklabels():
        label.set_fontsize(12)

    # Check if any discarded taxa is predicted as
    # extra category only present in training set
    common_elements = get_df_index_in_str_list(
        predicted_taxonomic_redistribution_df, train_kept_train_only_newName_classes_list
    )

    if len(common_elements) > 0:
        red_patch = Patch(color="red", label="Extra training\nclasses")
        ax.legend(
            handles=[red_patch],
            fontsize=10,
            loc="upper right",
            bbox_to_anchor=(-0.15, -1.05, 0.15, 1),
            bbox_transform=ax.transAxes,
            borderaxespad=0.8,
        )

    return fig


def create_classification_figure(job_data: dict) -> Figure:
    """
    Create the classification figure for retained samples.

    Parameters
    ----------
    job_data : dict
        Dictionary containing job data from extract_job_data

    Returns
    -------
    matplotlib.figure.Figure
        The generated classification figure
    """
    with plt.style.context("fivethirtyeight"):
        fig = built_classif_figure_by_process(
            job_data["confusion_matrix"],
            job_data["classif_report"],
            job_data["conf_mtrx_fig_params"],
            job_data["retained_samples_data"]["newName"],
            job_data["retained_samples_data"]["newName"],
            job_data["retained_samples_data"]["category_frequency_labels"].to_list(),
        )
        plt.tight_layout()
        return fig


def create_discarded_categories_heatmap(
    job_data: dict, taxo_category_mapping_df: pd.DataFrame
) -> Optional[Figure]:
    """
    Create heatmap showing redistribution of discarded categories.

    Parameters
    ----------
    job_data : dict
        Dictionary containing job data from extract_job_data
    taxo_category_mapping_df : pandas.DataFrame
        DataFrame mapping taxonomy categories

    Returns
    -------
    matplotlib.figure.Figure or None
        The generated heatmap figure, or None if no discarded classes
    """
    discarded_classes_df = job_data["discarded_classes_df"]
    if discarded_classes_df.empty:
        return None

    # Get discarded samples data
    discarded_samples_data = job_data["samples_data"].loc[
        ~job_data["samples_data"]["learned_category"]
    ]

    # Get taxonomic redistribution data
    predicted_discarded_taxo_redistribution_df = get_predicted_taxonomic_redistribution_df(
        taxo_category_mapping_df,
        discarded_classes_df,
    )

    # Create heatmap figure
    with plt.style.context("fivethirtyeight"):
        fig = plot_discarded_categ_redistributed_predictions_heatmap(
            predicted_discarded_taxo_redistribution_df,
            discarded_samples_data,
            job_data["retained_train_only_classes"],
        )
        plt.tight_layout()
        return fig


def create_relative_abundance_figure(
    prediction_df: pd.DataFrame, redistributed: bool = False
) -> Optional[Figure]:
    """
    Create relative abundance figure for original or redistributed predictions.

    Parameters
    ----------
    prediction_df : pandas.DataFrame
        DataFrame containing prediction data
    redistributed : bool, default=False
        Whether to use redistributed column names

    Returns
    -------
    matplotlib.figure.Figure or None
        The generated relative abundance figure or None if required columns don't exist
    """
    if redistributed:
        col1, col2 = "object_newname_redistr", "predicted_newName_redistr"
        # Check if redistributed columns exist
        if col1 not in prediction_df.columns or col2 not in prediction_df.columns:
            print(f"Warning: Redistributed columns {col1} and/or {col2} not found in DataFrame")
            return None
    else:
        col1, col2 = "object_newname", "predicted_newName"

    # Compute relative abundance data
    top_combined, top_rel_abund1, top_rel_abund2 = compute_dual_top_class_relat_abund_dfs(
        prediction_df, col1, col2, "short_sample_id"
    )

    # Create relative abundance figure
    with plt.style.context("seaborn-v0_8-pastel"):
        fig = plot_dual_relative_abundance(
            top_combined, top_rel_abund1, top_rel_abund2, redistributed
        )
        plt.tight_layout()
        return fig


def create_region_scn_maxObj_comp_strategy_metrics_figs(
    metrics_df: pd.DataFrame, results_path: Path, show_figure: bool = False
) -> None:
    """
    Generate and save comparative plots of classifier metrics (precision, recall, f1-score)
    across different classifier strategies, SCN feature usage, and learning object limits for each region.

    For each unique combination of region, SCN feature usage, and learning object limit, this function:
      - Filters the metrics DataFrame for the relevant subset.
      - Plots the macro and weighted averages of precision, recall, and f1-score for each strategy.
      - Saves the resulting figure as a PDF in a region-specific subfolder.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        MultiIndex DataFrame containing classifier metrics, including strategy, region, SCN usage, and object limits.
    results_path : Path
        Base directory where the generated figures will be saved, organized by region.
    show_figure : bool, optional
        If True, displays each figure interactively; otherwise, figures are only saved to disk.

    Returns
    -------
    None
    """
    # Extract unique parameter values
    df_no_top_level = metrics_df.droplevel(level=0, axis=1)
    regions = list(df_no_top_level.loc[:, "region"].unique())
    use_scn_cases = list(df_no_top_level.loc[:, "use_scn_features"].unique())
    learning_objs_limits = list(df_no_top_level.loc[:, "learning_objs_limit"].unique())

    # Retain only necessary columns
    metrics_df = metrics_df[["Strategy", "macro avg", "weighted avg"]]

    # Create plots for each parameter combination
    for region in regions:
        # Create region's figure directory
        region_figures_path: Path = results_path / region / "figures"
        region_figures_path.mkdir(parents=True, exist_ok=True)

        # Filter data for this region
        df_region = metrics_df.loc[df_no_top_level["region"] == region, :]

        for scn_use in use_scn_cases:
            # Filter for SCN feature usage
            df_region_scn = df_region.loc[df_no_top_level["use_scn_features"] == scn_use, :]

            for max_objs in learning_objs_limits:
                # Filter for learning objects limit
                df_filtered = df_region_scn.loc[
                    df_no_top_level["learning_objs_limit"] == max_objs, :
                ]

                # Skip if no data for this combination
                if df_filtered.empty:
                    continue

                # Create the figure
                fig, axes = plot_region_scn_maxObj_comp_metrics_fig(
                    region, scn_use, max_objs, df_filtered
                )

                # Format axes and save figure
                save_region_scn_maxObj_comp_metrics_fig_as_pdf(
                    fig, region_figures_path, region, scn_use, max_objs, show_figure
                )
    if not show_figure:
        plt.close("all")


def plot_region_scn_maxObj_comp_metrics_fig(
    region: str, scn_use: bool, max_objs: int, data: pd.DataFrame
) -> tuple:
    """
    Create a matplotlib figure with three vertically stacked subplots, each showing
    macro and weighted averages for precision, recall, and f1-score across classifier strategies.

    Parameters
    ----------
    region : str
        Name of the region (used for figure width and titles).
    data : pd.DataFrame
        DataFrame containing strategy labels and corresponding macro/weighted metrics.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created matplotlib Figure object.
    axes : np.ndarray
        Array of Axes objects for the three subplots (precision, recall, f1-score).
    """
    # Determine figure width based on region
    fig_width = 20 if region == "Gulf" else 13

    # Create figure and axes with the fivethirtyeight style
    with plt.style.context("fivethirtyeight"):
        fig, axes = plt.subplots(3, 1, figsize=(fig_width, 8), sharex=True)

        # Plot metrics for each score type
        for i, score_type in enumerate(["precision", "recall", "f1-score"]):
            # Get data
            x = data[("Strategy", "strategy_label")].to_list()
            y1 = data[("macro avg", score_type)].to_list()
            y2 = data[("weighted avg", score_type)].to_list()

            # Create plots
            axes[i].plot(x, y1, marker="o", label="macro avg")
            axes[i].plot(x, y2, marker="o", label="weighted avg")

            axes[i].set_title(score_type)
            axes[i].set_ylabel(f"{score_type}")
            axes[i].set_ylim((0, 1))
            fig_title = (
                f"region: {region} - use of SCN features: {'Yes' if scn_use else 'No'} - "
                f"max learning objects/class: {'maximum available' if max_objs == 999999 else max_objs}"
            )
            fig.suptitle(fig_title, fontsize=20)
            plt.xticks(
                fontsize=10,
                rotation=-35,
                ha="left",
                va="top",
                rotation_mode="anchor",
            )
            for j, label in enumerate(axes[i].get_xticklabels()):
                if j % 2 == 1:
                    label.set_color("dimgray")

        # Add legend to middle plot
        axes[1].legend(["macro avg", "weighted avg"], loc="center left", bbox_to_anchor=(1.02, 0.5))

        plt.tight_layout()

        return fig, axes


def save_region_scn_maxObj_comp_metrics_fig_as_pdf(
    fig,
    output_path: Path,
    region: str,
    scn_use: bool,
    max_objs: int,
    show_figure: bool = False,
) -> None:
    """
    Save a matplotlib figure as a PDF file in the specified output directory, with a filename
    encoding the region, SCN feature usage, and learning object limit.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    output_path : Path
        Directory where the PDF will be saved.
    region : str
        Name of the region (used in the filename).
    scn_use : bool
        Whether SCN features were used (used in the filename).
    max_objs : int
        Maximum number of learning objects per class (used in the filename).
    show_figure : bool, optional
        If True, displays the figure interactively before saving.

    Returns
    -------
    None
    """
    if show_figure:
        plt.show()

    # Create output filename
    filename = (
        f"{region}-{'With SCN' if scn_use else 'Without SCN'}"
        f"-{max_objs if max_objs < 999999 else 'Max'} Objects"
        f"-classifier_strategy_detailed_scores_comp.pdf"
    )
    output_filepath = output_path / filename

    # Save figure to PDF
    with PdfPages(output_filepath) as pdf_file:
        pdf_file.savefig(fig)


def create_f1_scores_strategy_comp_by_scn_maxObj_figs(
    metrics_df: pd.DataFrame,
    results_folderpath: Path,
    score_type: str = "f1-score",
    show_figure: bool = False,
) -> None:
    """
    Plot F1-scores for classifiers across different regions, SCN feature usage, and object limits.

    Parameters
    ----------
    metrics_df : pandas.DataFrame
        DataFrame containing classifier metrics. Must have a MultiIndex with
        ('Strategy', 'region'), ('Strategy', 'use_scn_features'), and
        ('Strategy', 'learning_objs_limit') as levels, and columns including
        ('Strategy', 'strategy_label'), ('macro avg', 'f1-score'), and
        ('weighted avg', 'f1-score').
    results_folderpath : pathlib.Path
        Path to the folder where the generated plots (PDF files) will be saved.
        Plots are saved in subdirectories by region.
    score_type : str, optional
        The score metric to plot on the y-axis. Default is "f1-score".
        Should match the column name in the metrics DataFrame.
    show_figure : bool, optional
        If True, display the generated figures interactively. If False (default),
        figures are closed after saving.

    Returns
    -------
    None
        This function does not return any value. It saves plots as PDF files
        in the specified results folder.

    Notes
    -----
    - For each region, a grid of subplots is created, with rows corresponding to
      different object limits and columns to SCN feature usage.
    - Each subplot shows the macro and weighted average F1-scores for all strategies.
    - Legends are added to the figure, and plots are saved as PDF files in
      region-specific subfolders under `results_folderpath`.
    """
    # Set up indexed dataframe for easier access
    metrics_df.sort_index(inplace=True)
    indexed_metrics_df = metrics_df.set_index(
        [
            ("Strategy", "region"),
            ("Strategy", "use_scn_features"),
            ("Strategy", "learning_objs_limit"),
        ]
    )

    # Get unique values for each dimension
    regions = list(indexed_metrics_df.index.get_level_values(0).unique())
    scn_usage = list(indexed_metrics_df.index.get_level_values(1).unique())
    max_objs = list(indexed_metrics_df.index.get_level_values(2).unique())

    # Create plots for each region
    for region in regions:
        nrows, ncols = len(max_objs), len(scn_usage)
        fig_width = 30 if region == "Gulf" else 20

        with plt.style.context("fivethirtyeight"):
            # Create figure and axes
            fig, axes = plt.subplots(
                nrows, ncols, figsize=(fig_width, 8), sharex=True, constrained_layout=True
            )

            # Convert axes to 2D array regardless of dimensions
            if nrows == 1 and ncols == 1:
                # Single subplot case
                axes = np.array([[axes]])
            elif nrows == 1:
                # Single row case
                axes = np.array([axes])
            elif ncols == 1:
                # Single column case
                axes = np.array([ax.reshape(-1, 1) for ax in axes])

            # Set main figure title
            fig_title = (
                f"Comparison of classifier's F1-Score for the {region} region"
                " with different learning strategies"
            )
            fig.suptitle(fig_title, fontsize=20)

            # Create subplots for each combination of parameters
            for j, use_scn in enumerate(scn_usage):
                for i, max_obj in enumerate(max_objs):
                    create_f1_scores_strategy_comp_subplot(
                        axes[i, j], indexed_metrics_df, region, use_scn, max_obj, score_type
                    )

            # Add legend
            fig.legend(
                ["macro avg", "weighted avg"], loc="center left", bbox_to_anchor=(0.90, 0.92)
            )

            # Finalize and save plot
            plt.tight_layout()
            if show_figure:
                plt.show()

            # Save to PDF
            save_f1_scores_strategy_comp_fig_to_pdf(fig, region, results_folderpath)

            if not show_figure:
                plt.close("all")


def create_f1_scores_strategy_comp_subplot(
    ax, indexed_metrics_df, region, use_scn, max_objs, score_type
):
    """
    Create a single subplot displaying F1-score data for a specific combination of parameters.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object on which to plot the data.
    indexed_metrics_df : pandas.DataFrame
        MultiIndexed DataFrame containing classifier metrics, indexed by region, SCN feature usage,
        and object limit.
    region : str
        The region for which the subplot is being created.
    use_scn : bool
        Indicates whether SCN features are used (True) or not (False).
    max_objs : int
        The maximum number of objects per class used for training. If set to 999999, it is treated as 'maximum'.
    score_type : str
        The type of score to plot on the y-axis (e.g., "f1-score", "precision", "recall").
    """
    try:
        # Get data for this combination of parameters
        data = indexed_metrics_df.loc[
            (region, use_scn, max_objs),
            [
                ("Strategy", "strategy_label"),
                ("Strategy", "strategy_numb"),
                ("macro avg", "f1-score"),
                ("weighted avg", "f1-score"),
            ],
        ]

        # Only plot if data exists and is not empty
        if not data.empty:
            ax.plot(
                data[("Strategy", "strategy_label")],
                data[("macro avg", score_type)],
                marker="o",
                label="macro avg",
            )
            ax.plot(
                data[("Strategy", "strategy_label")],
                data[("weighted avg", score_type)],
                marker="o",
                label="weighted avg",
            )
    except KeyError:
        # Index doesn't exist - skip plotting for this combination
        pass

    # Set title and formatting regardless of whether data exists
    ax.set_title(
        f"Learning with {'maximum' if max_objs == 999999 else max_objs} objs/class - "
        f"{'Not u' if not use_scn else 'U'}sing SCN features",
        fontsize=12,
    )
    ax.set_ylabel(f"{score_type}")
    ax.set_ylim((0, 1))

    # Format x-axis labels
    for k, label in enumerate(ax.get_xticklabels()):
        if k % 2 == 1:
            label.set_color("dimgray")
        label.set_size(10)
        label.set_rotation(-35)
        label.set_rotation_mode("anchor")
        label.set_ha("left")
        label.set_va("top")


def save_f1_scores_strategy_comp_fig_to_pdf(
    fig: Figure, region: str, results_folderpath: Path
) -> None:
    """
    Save a matplotlib figure as a PDF file in a region-specific directory.

    Parameters
    ----------
    fig : Figure
        The matplotlib Figure object to be saved.
    region : str
        The name of the region, used to construct the output directory and filename.
    results_folderpath : Path
        The base directory where results are stored. The figure will be saved under
        `<results_folderpath>/<region>/figures/`.

    Returns
    -------
    None
        This function does not return anything. It saves the figure as a PDF file.
    """
    filename = f"{region}-classifier_learning_strategy_f1scores_comp.pdf"
    figures_path = results_folderpath / region / "figures"
    figures_path.mkdir(parents=True, exist_ok=True)

    pdf_filepath = figures_path / filename
    with PdfPages(pdf_filepath) as pdf_file:
        pdf_file.savefig(fig)
