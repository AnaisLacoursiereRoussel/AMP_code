import typing

from ecotaxa_categories_name_id_mapping import (batch_map_category_name_to_id,
                                                map_category_name_to_id)
from ecotaxa_py_client.models.export_req import ExportReq
from ecotaxa_py_client.models.prediction_req import PredictionReq
from ecotaxa_py_client.models.project_filters import ProjectFilters


def prediction_request_factory(
    project_id: int,
    source_project_ids: list[int],
    learning_limit: int,
    features: list[str],
    source_taxo_categories: list[int],
    use_scn: bool,
    pre_mapping: dict[str, int] = {},
) -> PredictionReq:
    """
    Create a PredictionReq object with the specified parameters.

    Parameters
    ----------
    project_id : int
        The ID of the project for which predictions are to be made.
    source_project_ids : list of int
        List of source project IDs to use for training.
    learning_limit : int
        Maximum number of learning objects to use.
    features : list of str
        List of feature names to use for prediction.
    source_taxo_categories : list of int
        List of taxonomic category IDs to use for learning.
    use_scn : bool
        Whether to use SCN (special category normalization) flags.
    pre_mapping : dict of str to int, optional
        Mapping from replaced category IDs (as strings) to replacing category IDs.

    Returns
    -------
    PredictionReq
        Configured PredictionReq object.
    """
    prediction_request = PredictionReq(
        project_id=project_id,
        source_project_ids=source_project_ids,
        learning_limit=learning_limit,
        features=features,
        categories=source_taxo_categories,
        use_scn=use_scn,
        pre_mapping=pre_mapping,
    )
    return prediction_request


def export_request_factory(
    project_id: int,
) -> ExportReq:
    """
    Create an ExportReq object with the specified project ID and default export settings.

    Parameters
    ----------
    project_id : int
        The ID of the project to export.

    Returns
    -------
    ExportReq
        Configured ExportReq object with default parameters.
    """
    export_request = ExportReq(
        project_id=project_id,
        exp_type="TSV",
        use_latin1=True,
        tsv_entities="OS",
        split_by="",
        coma_as_separator=False,
        format_dates_times=False,
        with_images=False,
        with_internal_ids=True,
        only_first_image=False,
        sum_subtotal="",
        pre_mapping={},
        formulae={},
        out_to_ftp=False,
    )
    return export_request


def project_filters_factory(
    taxo: str = "",
    taxochild: str = "",
) -> ProjectFilters:
    """
    Create a ProjectFilters object with the specified taxonomic filters and default values for other filters.

    Parameters
    ----------
    taxo : str, optional
        Taxonomic filter (default is an empty string).
    taxochild : str, optional
        Child taxonomic filter (default is an empty string).

    Returns
    -------
    ProjectFilters
        Configured ProjectFilters object.
    """
    project_filters = ProjectFilters(
        taxo=taxo,
        taxochild=taxochild,
        statusfilter="",
        depthmin="",
        depthmax="",
        samples="",
        instrum="",
        daytime="",
        month="",
        fromdate="",
        todate="",
        fromtime="",
        totime="",
        inverttime="",
        validfromdate="",
        validtodate="",
        freenum="",
        freenumst="",
        freenumend="",
        freetxt="",
        freetxtval="",
        filt_annot="",
        filt_last_annot="",
    )
    return project_filters


def prediction_scenario_factory(
    strategy_settings_dict: dict[str, typing.Any],
) -> dict[str, typing.Any]:
    """
    Construct a prediction scenario dictionary with all required request and filter objects.

    This function processes the strategy settings, maps category names to IDs,
    handles taxonomic category name mappings, and creates all necessary request and filter
    objects for a prediction scenario.

    Parameters
    ----------
    strategy_settings_dict : dict of str to Any
        Dictionary containing strategy settings, including:
            - "training_classes_list": list of str, category names for training
            - "taxo_category_names_mapping_dict": dict of str to str, mapping of replaced to replacing category names
            - "project_id": int, project ID
            - "source_project_ids": list of int, source project IDs
            - "learning_objs_limit": int, learning object limit
            - "use_scn_flags": bool, use SCN flags
            - "features": list of str, feature names

    Returns
    -------
    dict of str to Any
        Dictionary containing:
            - "prediction_scenario": updated strategy_settings_dict
            - "prediction_request": PredictionReq object
            - "prediction_project_filters": ProjectFilters object
            - "export_request": ExportReq object
            - "export_project_filters": ProjectFilters object

    Raises
    ------
    KeyError
        If a replaced or replacing category name cannot be mapped to a valid ID.
    """
    learning_category_id_list = batch_map_category_name_to_id(
        strategy_settings_dict["training_classes_list"]
    )
    strategy_settings_dict["learning_category_id_list"] = learning_category_id_list

    if strategy_settings_dict["taxo_category_names_mapping_dict"]:
        taxo_category_ids_mapping: dict[str, int] = {}
        for replaced_name, replacing_name in strategy_settings_dict["taxo_category_names_mapping_dict"].items():
            replaced_id = map_category_name_to_id(replaced_name)
            replacing_id = map_category_name_to_id(replacing_name)
            if (replaced_id + replacing_id) > 1:
                taxo_category_ids_mapping[str(replaced_id)] = replacing_id
            else:
                raise KeyError()
    else:
        taxo_category_ids_mapping = {}
    strategy_settings_dict["taxo_category_ids_mapping"] = taxo_category_ids_mapping

    prediction_request = prediction_request_factory(
        project_id=strategy_settings_dict["project_id"],
        source_project_ids=strategy_settings_dict["source_project_ids"],
        learning_limit=strategy_settings_dict["learning_objs_limit"],
        source_taxo_categories=strategy_settings_dict["learning_category_id_list"],
        use_scn=strategy_settings_dict["use_scn_flags"],
        features=strategy_settings_dict["features"],
        pre_mapping=strategy_settings_dict["taxo_category_ids_mapping"],
    )
    prediction_project_filters = project_filters_factory()
    export_request = export_request_factory(project_id=strategy_settings_dict["project_id"])
    export_project_filters = project_filters_factory()

    return {
        "prediction_scenario": strategy_settings_dict,
        "prediction_request": prediction_request,
        "prediction_project_filters": prediction_project_filters,
        "export_request": export_request,
        "export_project_filters": export_project_filters,
    }
