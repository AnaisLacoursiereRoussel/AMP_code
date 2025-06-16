"""
Provides access to functions to request prediction, job status and job file
to and from the EcoTaxa API.
"""

import json
import pathlib
import shutil
import time
from typing import Any, Optional

import ecotaxa_py_client  # type: ignore
from ecotaxa_py_client.api import authentification_api, jobs_api, objects_api
from ecotaxa_py_client.models.body_export_object_set_object_set_export_post import (
    BodyExportObjectSetObjectSetExportPost,
)
from ecotaxa_py_client.models.body_predict_object_set_object_set_predict_post import (
    BodyPredictObjectSetObjectSetPredictPost,
)
from ecotaxa_py_client.models.export_req import ExportReq
from ecotaxa_py_client.models.login_req import LoginReq
from ecotaxa_py_client.models.prediction_req import PredictionReq
from ecotaxa_py_client.models.project_filters import ProjectFilters
from predictions_scenarios_factory import prediction_scenario_factory

predict_export_folderpath = pathlib.Path(r"exported_classification_results")

FLOWCAM_IMAGES_FEATURES = [
    "fre.abd_area",
    "fre.abd_diameter",
    "fre.abd_volume",
    "fre.biovolume_cylind",
    "fre.biovolume_Pspheroid",
    "fre.biovolume_sphere",
    "fre.circle_fit",
    "fre.circularity_hu",
    "fre.compactness",
    "fre.convex_perimeter",
    "fre.edge_gradient",
    "fre.elongation",
    "fre.esd_diameter",
    "fre.esd_volume",
    "fre.fd_diameter",
    "fre.feret_max_angle",
    "fre.feret_min_angle",
    "fre.filled_area",
    "fre.geodesic_aspectratio",
    "fre.geodesic_length",
    "fre.geodesic_thickness",
    "fre.image_h",
    "fre.image_w",
    "fre.intensity",
    "fre.length",
    "fre.perimeter",
    "fre.raw_area",
    "fre.raw_convex_hull_area",
    "fre.raw_convex_perimeter",
    "fre.raw_feret_max",
    "fre.raw_feret_mean",
    "fre.raw_feret_min",
    "fre.raw_filled_area",
    "fre.raw_legendre_major",
    "fre.raw_legendre_minor",
    "fre.raw_perimeter",
    "fre.roughness",
    "fre.sigma_intensity",
    "fre.sum_intensity",
    "fre.symmetry",
    "fre.width",
]


FLOWCAM_IMAGES_FEATURES = [feat.lower() for feat in FLOWCAM_IMAGES_FEATURES]


def get_yes_or_no_input(prompt: str) -> str:
    """
    Prompt the user for a yes or no response.

    Parameters
    ----------
    prompt : str
        The message displayed to the user.

    Returns
    -------
    str
        The user's response, either 'y' or 'n' (lowercase).

    Notes
    -----
    The function will continue to prompt until a valid response is given.
    """
    while True:
        response = input(prompt).lower()  # Convert to lowercase to accept 'Y' or 'N' as well
        if response in ["y", "n"]:
            return response
        else:
            print("Please enter 'y' for yes or 'n' for no.")


def setup_api_configuration_token() -> None:
    """
    Set up the EcoTaxa API configuration by authenticating the user and retrieving an access token.

    Prompts the user for EcoTaxa credentials, performs login, and sets the access token
    in the global configuration object.

    Returns
    -------
    None

    Raises
    ------
    ecotaxa_py_client.ApiException
        If authentication fails or the API call encounters an error.
    """
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
        # Create an instance of the API class
        api_instance = authentification_api.AuthentificationApi(api_client)
        # ask user for credentials
        print("Please enter your EcoTaxa credentials to login:")
        username = input("Username: ")
        password = input("Password: ")
        # setup your credentials
        login_req = LoginReq(
            password=password,
            username=username,
        )

        # example passing only required values which don't have defaults set
        try:
            # Login
            api_response_token = api_instance.login(login_req)
            print("api_response_token", api_response_token)
            # set YOUR_ACCESS_TOKEN configuration
            configuration.access_token = api_response_token
        except ecotaxa_py_client.ApiException as e:
            print("Exception when calling AuthentificationApi->login: %s\n" % e)


def predict_and_record_objects_zipfile(
    prediction_scenario: dict[Any, Any],
    prediction_request: PredictionReq,
    prediction_project_filters: ProjectFilters,
    export_request: ExportReq,
    export_project_filters: ProjectFilters,
    predict_export_folderpath: pathlib.Path = predict_export_folderpath,
) -> tuple[Optional[str], dict[Any, Any]]:
    """
    Request prediction on a project with a set of prediction parameters,
    followed by an export request and a download of the export zipfile.

    Parameters
    ----------
    prediction_scenario : dict
        Dictionary containing scenario-specific parameters.
    prediction_request : PredictionReq
        The prediction request object for the API.
    prediction_project_filters : ProjectFilters
        Filters to apply to the prediction request.
    export_request : ExportReq
        The export request object for the API.
    export_project_filters : ProjectFilters
        Filters to apply to the export request.
    predict_export_folderpath : pathlib.Path, optional
        Path to the folder where the exported zipfile will be saved (default is global variable).

    Returns
    -------
    tuple
        A tuple containing:
        - The file path to the exported zipfile (str) if successful, otherwise None.
        - The prediction scenario dictionary.

    Notes
    -----
    This function handles the full workflow: prediction, export, job status polling, and file download.
    """
    # Enter a context with an instance of the API client
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
        # Create an instance of the API class
        api_instance = objects_api.ObjectsApi(api_client)
        body_predict_object_set_object_set_predict_post = BodyPredictObjectSetObjectSetPredictPost(
            filters=prediction_project_filters,
            request=prediction_request,
        )

        try:
            # Predict Object Set
            prediction_api_response = api_instance.predict_object_set(
                body_predict_object_set_object_set_predict_post
            )
            prediction_errors = prediction_api_response.errors
            prediction_warnings = prediction_api_response.warnings
            prediction_job_id = prediction_api_response.job_id
            print(
                "Prediction request server response:",
                f"prediction_errors: {prediction_errors}",
                f"prediction_job_id: {prediction_job_id}",
                f"prediction_warnings: {prediction_warnings}",
                sep="\n",
            )
        except ecotaxa_py_client.ApiException as e:
            print("Exception when calling ObjectsApi->predict_object_set: %s\n" % e)

        if prediction_job_id:
            # Create an instance of the API class
            prediction_job_state = request_job_state(api_client, prediction_job_id)
        else:
            prediction_job_state = "E"

        if prediction_job_state == "F":
            # Create an instance of the API class
            api_instance = objects_api.ObjectsApi(api_client)
            body_export_object_set_object_set_export_post = BodyExportObjectSetObjectSetExportPost(
                filters=export_project_filters,
                request=export_request,
            )

            # example passing only required values which don't have defaults set
            try:
                # Export Object Set
                export_api_response = api_instance.export_object_set(
                    body_export_object_set_object_set_export_post
                )
                export_errors = export_api_response.errors
                export_warnings = export_api_response.warnings
                export_job_id = export_api_response.job_id
                print(
                    "Export query server response",
                    f"export_errors: {export_errors}",
                    f"export_job_id: {export_job_id}",
                    f"export_warnings: {export_warnings}",
                    sep="\n",
                )
            except ecotaxa_py_client.ApiException as e:
                print("Exception when calling ObjectsApi->export_object_set: %s\n" % e)

        if export_job_id:
            # Create an instance of the API class
            export_job_state = request_job_state(api_client, export_job_id)
        else:
            export_job_state = "E"

        if export_job_state == "F":
            # Create an instance of the JobsApi class
            jobs_api_instance = jobs_api.JobsApi(api_client)
            job_id = export_job_id
            if job_id is None:
                raise ValueError("export_job_id is None, cannot retrieve job file.")
            # example passing only required values which don't have defaults set
            try:
                # Get Job File
                get_export_file_api_response = jobs_api_instance.get_job_file_without_preload_content(
                    job_id
                )
                # get_export_file_api_response = jobs_api_instance.get_job_file(
                #     job_id
                # )
                print(get_export_file_api_response)

                prediction_job_filepath = (
                    predict_export_folderpath
                    / f"exportjob{export_job_id}-predictjob{prediction_job_id}.zip"
                )
                with open(
                    prediction_job_filepath,
                    "wb",
                ) as out_file:
                    shutil.copyfileobj(get_export_file_api_response, out_file)
                return str(prediction_job_filepath), prediction_scenario
            except ecotaxa_py_client.ApiException as e:
                print("Exception when calling JobsApi->get_job_file: %s\n" % e)

        return None, prediction_scenario


def request_job_state(api_client: ecotaxa_py_client.ApiClient, job_id: int) -> Optional[str]:
    """
    Poll the EcoTaxa API for the status of a job until it finishes, errors, or times out.

    Parameters
    ----------
    api_client : ecotaxa_py_client.ApiClient
        The API client instance to use for requests.
    job_id : int
        The ID of the job to poll.

    Returns
    -------
    Optional[str]
        The final state of the job ('F' for finished, 'E' for error, 'A' for question, or a timeout message).

    Notes
    -----
    Will prompt the user if the job takes longer than 10 minutes.
    """
    api_instance = jobs_api.JobsApi(api_client)

    iteration_number = 0
    while True:
        try:
            # Get Job
            get_job_api_response = api_instance.get_job(job_id)
            job_state: Optional[str] = get_job_api_response.state
            job_step = get_job_api_response.step
            job_progress_pct = get_job_api_response.progress_pct
            print(
                "Job state request server response:",
                f"job_state: {job_state}",
                f"job_step: {job_step}",
                f"job_progress_pct: {job_progress_pct}",
                sep="\n",
            )
            if iteration_number > 50:
                print("Time out, the job took more than 10 min")
                job_state = (
                    f"Timed out. last job state = {job_state}, progress = {job_progress_pct}"
                )
                response = get_yes_or_no_input("Stop process of querying EcoTaxa (Y/N)?")
                if response.lower() == "y":
                    break
                elif response.lower() == "n":
                    continue
            elif job_state == "F":
                print(f"job# {job_id} has finished")
                break
            elif job_state == "E":
                print(f"An error as occured and stopped job# {job_id}")
                break
            elif job_state == "A":
                print(f"Question asked and stopped job# {job_id}")
                print(get_job_api_response.question)
                break
            else:
                time.sleep(15)
                iteration_number += 1

        except ecotaxa_py_client.ApiException as e:
            print("Exception when calling JobsApi->get_job: %s\n" % e)
            break
    return job_state


def queryEcotaxaAPI(
    strategy_settings_dict: dict[str, Any],
    learning_objs_limit: int,
    features: list[str],
    use_scn_flags: bool,
) -> None:
    """
    Run a prediction and export scenario on the EcoTaxa API using the provided settings.

    Parameters
    ----------
    strategy_settings_dict : dict
        Dictionary containing scenario strategy settings.
    learning_objs_limit : int
        Maximum number of learning objects to use.
    features : list of str
        List of feature names to include in the prediction.
    use_scn_flags : bool
        Whether to use scenario-specific flags.

    Returns
    -------
    None

    Notes
    -----
    The function saves the exported zipfile and a JSON file with scenario and result details.
    """
    strategy_settings_dict["learning_objs_limit"] = learning_objs_limit
    strategy_settings_dict["features"] = features
    strategy_settings_dict["use_scn_flags"] = use_scn_flags

    prediction_scenario = prediction_scenario_factory(
        strategy_settings_dict=strategy_settings_dict,
    )
    prediction_results: tuple[Optional[str], dict[str, str]] = predict_and_record_objects_zipfile(
        **prediction_scenario
    )
    # prediction_scenario_list.append(prediction_results)
    print("exported zip and appended scenario")

    prediction_jobs_data_json_filepath = (
        (f"{prediction_results[0].strip('.zip')}_settings.json")
        if prediction_results[0] is not None
        else None
    )
    print(prediction_jobs_data_json_filepath)

    if prediction_jobs_data_json_filepath:
        with open(prediction_jobs_data_json_filepath, "w", encoding="utf-8") as f:
            json.dump(prediction_results, f, ensure_ascii=False, indent=4)


if __name__ == "__main__":
    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api")
    # In case of ssl error
    configuration.verify_ssl = True

    # load prediction scenario JSON file
    scenarios_details_json_filepath = (
        pathlib.Path("scenarios_details_JSON") / "scenarios_implementation_details.json"
    )

    with open(scenarios_details_json_filepath, "r", encoding="utf-8") as f:
        scenarios_details_json: dict = json.load(f)

    region_list = ["Gulf", "NL 2020", "NL 2021", "PA"]
    use_scn_features_list = [False, True]
    learning_objs_limit_list = [200, 5000, 20000, 999999]

    setup_api_configuration_token()
    n = 0

    for use_scn in use_scn_features_list:
        for max_objs in learning_objs_limit_list:
            for query_setting in list(scenarios_details_json.values()):
                n += 1
                print(n, "of", len(scenarios_details_json), "scenarios")
                queryEcotaxaAPI(
                    strategy_settings_dict=query_setting,
                    learning_objs_limit=max_objs,
                    features=FLOWCAM_IMAGES_FEATURES,
                    use_scn_flags=use_scn,
                )
