"""
This script remaps properties from an old physical property dataset to a new one by index.

This script is a helper utility to copy previously-benchmarked properties from the rc1 fit
to the rc2 fit to avoid re-computing data.
"""

import json
import pathlib
import sys
import click
from loguru import logger

import openff.evaluator
from openff.evaluator.datasets import PhysicalPropertyDataSet

logger.remove()
logger.add(sys.stdout)

@click.command(help=__doc__)
@click.option(
    "--old-dataset-file",
    "-old",
    default="../../../ash-sage-rc1/01_curate-physprop/validation/final/output/validation-set.json",
    help=(
        "The JSON file containing the old dataset. "
        "This should be a valid PhysicalPropertyDataSet JSON file."
    ),
)
@click.option(
    "--new-dataset-file",
    "-new",
    default="../../01_curate-physprop/validation/final/output/validation-set.json",
    help=(
        "The JSON file containing the new dataset. "
        "This should be a valid PhysicalPropertyDataSet JSON file."
    ),
)
@click.option(
    "--output-file",
    "-o",
    default="mappings/old-validation-to-new-validation.json",
    help=(
        "The JSON file to write the output mapping to."
    ),
)
def main(
    old_dataset_file: str,
    new_dataset_file: str,
    output_file: str,
):
    # this depends on a branch of main with property hashing
    logger.info(f"Using Evaluator version {openff.evaluator.__version__}")

    old_dataset = PhysicalPropertyDataSet.from_json(old_dataset_file)
    new_dataset = PhysicalPropertyDataSet.from_json(new_dataset_file)

    logger.info(f"Loaded old dataset with {len(old_dataset.properties)} properties.")
    logger.info(f"Loaded new dataset with {len(new_dataset.properties)} properties.")

    # start by hashing properties
    old_hash_to_index_map: dict[int, int] = {}
    for i, prop in enumerate(old_dataset.properties):
        old_hash_to_index_map[prop.get_property_hash()] = i


    old_to_new_index_map: dict[int, int] = {}
    for j, prop in enumerate(new_dataset.properties):
        prop_hash = prop.get_property_hash()
        if prop_hash in old_hash_to_index_map:
            i = old_hash_to_index_map[prop_hash]
            old_to_new_index_map[i] = j

    logger.info(f"Found {len(old_to_new_index_map)} matching properties by hash.")

    # save to file
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(old_to_new_index_map, f, indent=4)
    logger.info(f"Saved mapping to {output_file}")



if __name__ == "__main__":
    main()
