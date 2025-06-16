# %% [markdown]
# ## Build Prediction Scenario Settings from Excel file configuration details into a JSON file

# %% [markdown]
# Import necessary libraries for data manipulation, file system operations, and type hinting.

# %%
import json
import pathlib
from typing import Union

import numpy as np
import pandas as pd

# %% [markdown]
# Define the path to the Excel configuration file. Then, load the 'strategies' sheet,
# which outlines the different prediction strategies. Subsequently, iterate through a
# predefined list of regional titles, loading each corresponding sheet from the Excel file.
# These sheets contain specific implementation details for each region and are stored in a
# dictionary of pandas DataFrames.

# %%
strategies_config_excel_filepath = r"regional_prediction_strategies.xlsx"

# Loading the strategies definition sheet
strategies_config_df = pd.read_excel(strategies_config_excel_filepath, sheet_name="strategies")

# Loading each regional implementation details sheet as a dataframe into a dictionary
regional_title_list = ["Gulf", "NL 2020", "NL 2021", "PA"]

region_implementation_config_dict: dict[str, pd.DataFrame] = {}

for region in regional_title_list:
    region_df = pd.read_excel(strategies_config_excel_filepath, sheet_name=region)
    region_implementation_config_dict[region] = region_df

# %% [markdown]
# For each strategy in the `strategies_config_df` build the variables needed
# to query the EcoTaxa API, using implementation details found in the
# corresponding region's dataframe for the corresponding strategy number.<br>
# -> Record the scenarios details into a JSON file

# %% [markdown]
# Define column labels for accessing data within the regional strategy DataFrames and a mapping for
# columns in the main `strategies_config_df`. This facilitates easier and more readable access to
# specific data points later in the script.

# %%
region_strategy_df_col_labels = [
    "newName",
    "stage",
    "ecotaxa_category",
    "selected_samples_classes_count",
    "region_training_classes_count",
    "all_regions_training_classes_count",
    "sample_presence",
    "regional_training_presence",
    "global_training_presence",
    "Broad category",
    "Category type",
    "comments",
]

strategy_cfg_col_map = {
    "Region": 0,
    "project_id": 1,
    "strategy_number": 2,
    "Broad classes strategy": 3,
    "source_project_ids": 4,
    "special case": 5,
    "Training class set": 6,
    "Select sples classes present": 7,
    "Extra regional training only classes": 8,
    "Low training threshold": 9,
    "Region classes learning-title": 10,
}
# %% [markdown]
# Iterate through each strategy defined in `strategies_config_df`. For every strategy,
# extract and structure its configuration details, including regional specifics, project IDs,
# and class lists for training and discarding. This information is compiled into a dictionary
# for each scenario, which is then stored in a main dictionary `region_strategy_scenarios_dict`
# keyed by a tuple of region and strategy number.

# %%
region_strategy_scenarios_dict: dict[str, dict[str, Union[str, int, float, dict, list]]] = {}
for strategy_record in strategies_config_df.to_records(index=False):
    scenario_data_dict: dict[str, Union[str, int, float, dict, list]] = {}

    # Get the strategy dataframe with the relevant data only
    region_df = region_implementation_config_dict[strategy_record[strategy_cfg_col_map["Region"]]]

    # Retrieve specific strategy df from the regional strategy df,
    # discarding taxonomic categories not present in either selected samples
    # or regional training sets.
    strategy_df = region_df.loc[
        (region_df["sample_presence"] == "yes")
        | (region_df["regional_training_presence"] == "yes"),
        region_strategy_df_col_labels + [strategy_record[strategy_cfg_col_map["strategy_number"]]],
    ]

    scenario_data_dict["region"] = str(strategy_record[strategy_cfg_col_map["Region"]])
    scenario_data_dict["project_id"] = int(strategy_record[strategy_cfg_col_map["project_id"]])
    scenario_data_dict["source_project_ids"] = [
        int(src_id)
        for src_id in str(strategy_record[strategy_cfg_col_map["source_project_ids"]]).split(", ")
    ]

    training_set_used_string = (
        strategy_record[strategy_cfg_col_map["Region"]]
        if strategy_record[strategy_cfg_col_map["Training class set"]] == "Regional only"
        else "all regions"
    )

    scenario_title = (
        f"{strategy_record[strategy_cfg_col_map['Region']]} Selected Samples prediction using "
        f"{training_set_used_string} training set,\n"
        f"{strategy_record[strategy_cfg_col_map['Region classes learning-title']]},\n"
        f"{strategy_record[strategy_cfg_col_map['Broad classes strategy']]}"
    )
    scenario_data_dict["scenario_title"] = scenario_title
    scenario_data_dict["taxo_category_names_mapping_dict"] = {}

    # Get the list of discarded training classes from the region's strategy
    # implementation dataframe
    discarded_training_classes_list = strategy_df.loc[
        strategy_df[strategy_record[strategy_cfg_col_map["strategy_number"]]] == "Discarded",
        "ecotaxa_category",
    ].tolist()

    scenario_data_dict["discarded_training_classes_list"] = discarded_training_classes_list
    # Get the list of training classes from the region's strategy
    # implementation dataframe
    training_classes_list = strategy_df.loc[
        strategy_df[strategy_record[strategy_cfg_col_map["strategy_number"]]] == "kept",
        "ecotaxa_category",
    ].tolist()

    scenario_data_dict["training_classes_list"] = training_classes_list
    scenario_data_dict["strategy_df"] = strategy_df.to_json()
    scenario_data_dict["strategy_scenario_record"] = {
        col_name: record_field
        for col_name, record_field in zip(list(strategies_config_df.columns), strategy_record)
    }
    region_strategy_scenarios_dict[
        (
            f"({strategy_record[strategy_cfg_col_map['Region']]}, "
            f"{strategy_record[strategy_cfg_col_map['strategy_number']]})"
        )
    ] = scenario_data_dict

# %% [markdown]
# Specify the output file path for the JSON file that will store the compiled scenario details.
# A custom JSON encoder class, `npEncoder`, is defined to handle specific numpy integer types (int32, int64)
# by converting them to standard Python integers, ensuring compatibility with JSON serialization. The compiled
# `region_strategy_scenarios_dict` is then written to the specified JSON file using this custom encoder.

# %%
scenarios_details_json_filepath = pathlib.Path("scenarios_details_JSON") / "scenarios_implementation_details.json"


class npEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


with open(scenarios_details_json_filepath, "w", encoding="utf-8") as f:
    json.dump(region_strategy_scenarios_dict, f, cls=npEncoder)


# %% [markdown]
# Define the path to the previously created JSON file containing scenario implementation details.
# Open and read this JSON file, loading its contents into a dictionary named `dic`. This step is
# typically for verification or further processing of the stored scenario configurations.

# %%
scenarios_details_json_filepath = (
    pathlib.Path("scenarios_details_JSON") / "scenarios_implementation_details.json"
)

with open(scenarios_details_json_filepath, "r", encoding="utf-8") as f:
    dic = json.load(f)


# %% [markdown]
# As a test, retrieve and display the details for the first strategy of the "Gulf" region from the loaded
# `dic` dictionary. This helps in verifying that the data has been correctly processed and stored.

# %%
# test with 1st gulf strategy
test_1st_gulf_strat_dict = dic["(Gulf, 1)"]
test_1st_gulf_strat_dict
