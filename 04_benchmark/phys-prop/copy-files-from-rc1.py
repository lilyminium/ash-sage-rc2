"""
This script copies physical property files from an old data set to a new data set
based on a provided mapping of old property indices to new property indices.

This is intended to be used with the output of `remap-from-rc1-results.py`.
"""

import json
import pathlib
import shutil
import sys

import click
from loguru import logger

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--input-directory",
    "-i",
    default="../../../ash-sage/04_benchmark/phys-prop/thermoml-validation/validation",
    help=(
        "The directory containing the files to copy. "
        "This should contain subdirectories rep-0, rep-1, etc. with force field subdirs."
    ),
)
@click.option(
    "--output-directory",
    "-o",
    default="validation",
    help=(
        "The directory to copy the files to. "
        "This will contain subdirectories rep-0, rep-1, etc. with force field subdirs."
    ),
)
@click.option(
    "--mapping-file",
    "-m",
    default="mappings/old-validation-to-new-validation.json",
    help=(
        "The JSON file containing the mapping from old property indices to new property indices. "
        "This should be a dictionary mapping old indices (int) to new indices (int)."
    ),
)
def main(
    input_directory: str = "../../../ash-sage/04_benchmark/phys-prop/thermoml-validation/validation",
    output_directory: str = "validation",
    mapping_file: str = "mappings/old-validation-to-new-validation.json",
):
    
    # we follow a dataset / rep / property_index structure
    with open(mapping_file, "r") as f:
        old_to_new_mapping = {
            int(k): v
            for k, v in json.load(f).items()
        }

    logger.info(f"Loaded mapping for {len(old_to_new_mapping)} properties.")

    input_directory = pathlib.Path(input_directory)
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    # copy files
    for old_index, new_index in old_to_new_mapping.items():
        old_files = list(input_directory.glob(f"rep-*/*/prop-{old_index:04d}.*"))
        for old_file in old_files:
            forcefield = old_file.parent.name  # ff-x
            rep_num = old_file.parent.parent.name  # rep-x
            new_dir = output_directory / rep_num / forcefield
            new_dir.mkdir(parents=True, exist_ok=True)
            new_file = new_dir / f"prop-{new_index:04d}{old_file.suffix}"
            # copy2 preserves metadata
            shutil.copy2(old_file, new_file)
            logger.info(f"Copied {old_file} to {new_file}")

if __name__ == "__main__":
    main()
