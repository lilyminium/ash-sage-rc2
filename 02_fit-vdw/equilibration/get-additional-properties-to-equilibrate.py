"""
This script selects amine properties from an existing dataset and checks which ones
are not yet equilibrated. It saves the non-equilibrated properties to a new dataset file
for further processing.
"""

import click

import pandas as pd

from eveq.storage.storage import LocalStoredEquilibrationData

from openff.evaluator.datasets.datasets import PhysicalPropertyDataSet


@click.command()
@click.option(
    "--existing-storage-path",
    "-esp",
    "existing_storage_path",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    default="../../data/stored_data",
    help="Path to the existing storage directory.",
)
@click.option(
    "--input-dataset-path",
    "-i",
    "input_dataset_path",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    default="../../01_download-data/physprop/final/intermediate/renamed-filtered.csv",
    help=(
        "Path to the intermediate CSV file containing filtered properties. "
        "This file is from the ash-sage-rc2 dataset."
    )
)
@click.option(
    "--output-name",
    "-o",
    "output_name",
    type=str,
    default="not-equilibrated",
    help="Name for the output dataset containing non-equilibrated properties.",
)
def main(
    existing_storage_path: str = "../../data/stored_data",
    input_dataset_path: str = "../../01_download-data/physprop/final/intermediate/renamed-filtered.csv",
    output_name: str = "not-equilibrated",
):
    storage = LocalStoredEquilibrationData(existing_storage_path)
    print(f"Number of objects in storage: {len(storage._cached_retrieved_objects)}")

    # load existing intermediate filtered set
    df = pd.read_csv(
        input_dataset_path,
        index_col=0
    )
    df["Id"] = df["Id"].astype(str)

    print(
        f"{len(df)} properties found"
    )

    dataset = PhysicalPropertyDataSet.from_pandas(df)

    not_equilibrated = []
    for physical_property in dataset.properties:
        if not storage.contains_all_property_boxes(physical_property):
            not_equilibrated.append(physical_property)

    print(
        f"{len(not_equilibrated)} properties not equilibrated"
    )

    not_equilibrated_dataset = PhysicalPropertyDataSet()
    not_equilibrated_dataset.add_properties(*not_equilibrated)
    with open(f"{output_name}.json", "w") as f:
        f.write(not_equilibrated_dataset.json())

    not_equilibrated_dataset_df = not_equilibrated_dataset.to_pandas()
    not_equilibrated_dataset_df.to_csv(f"{output_name}.csv")


if __name__ == "__main__":
    main()
