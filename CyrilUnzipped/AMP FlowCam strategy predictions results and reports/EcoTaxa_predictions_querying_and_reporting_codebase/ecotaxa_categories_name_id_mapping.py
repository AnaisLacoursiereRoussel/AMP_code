import typing
import pandas as pd
from typing import cast

# Load EcoTaxa taxonomic categories reference data
taxo_categories_ecotaxa_df = pd.read_csv(
    r"ecotaxa_server_taxoexport_20231012_145829.tsv",
    sep="\t",
    usecols=["id", "display_name"],
)


def map_category_name_to_id(category_name: str) -> int:
    """
    Map a taxonomic category display name to its corresponding EcoTaxa category ID.

    Parameters
    ----------
    category_name : str
        The display name of the taxonomic category.

    Returns
    -------
    int
        The ID number of the taxonomic category if found, otherwise 0.

    Raises
    ------
    None
        Prints a message and returns 0 if the category name is not found.
    """
    try:
        category_id = cast(
            int, taxo_categories_ecotaxa_df.set_index("display_name").loc[category_name, "id"]
        )
        return int(category_id)
    except KeyError:
        print(f"no taxo category found with the name : {category_name}")
        return 0


def map_id_to_category_name(category_id: int) -> typing.Optional[str]:
    """
    Map an EcoTaxa category ID to its corresponding taxonomic category display name.

    Parameters
    ----------
    category_id : int
        The ID number of the taxonomic category.

    Returns
    -------
    str or None
        The display name of the taxonomic category if found, otherwise None.

    Raises
    ------
    None
        Prints a message and returns None if the category ID is not found.
    """
    try:
        category_name = str(
            taxo_categories_ecotaxa_df.set_index("id").loc[category_id, "display_name"]
        )
        return category_name
    except KeyError:
        print(f"no taxo category found with the ID number : {category_id}")
        return None


def batch_map_category_name_to_id(category_names_list: list[str]) -> list[int]:
    """
    Map a list of taxonomic category display names to their corresponding EcoTaxa category IDs.

    Parameters
    ----------
    category_names_list : list of str
        List of display names of the taxonomic categories.

    Returns
    -------
    list of int
        List of ID numbers corresponding to the provided category names.

    Raises
    ------
    KeyError
        If any category name in the list is not found.
    """
    category_ids_list: list[int] = []
    for category_name in category_names_list:
        category_id = map_category_name_to_id(category_name)
        if category_id > 0:
            category_ids_list.append(category_id)
        else:
            raise KeyError()
    return category_ids_list


def batch_map_id_to_category_name(category_ids_list: list[int]) -> list[str]:
    """
    Map a list of EcoTaxa category IDs to their corresponding taxonomic category display names.

    Parameters
    ----------
    category_ids_list : list of int
        List of ID numbers of the taxonomic categories.

    Returns
    -------
    list of str
        List of display names corresponding to the provided category IDs.

    Raises
    ------
    KeyError
        If any category ID in the list is not found.
    """
    category_names_list: list[str] = []
    for category_id in category_ids_list:
        category_name = map_id_to_category_name(category_id)
        if category_name:
            category_names_list.append(category_name)
        else:
            raise KeyError()
    return category_names_list
