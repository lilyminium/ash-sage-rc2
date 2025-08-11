"""
This script renames the property IDs in a PhysicalPropertyDataSet CSV file
to a more compact format based on the property type and a hash of the property attributes.
It is intended to be run after the initial filtering of ThermoML data (just to save time).
"""
import hashlib
import json
import logging
import sys

import click
import tqdm

import pandas as pd
from openff.evaluator.utils.serialization import TypedJSONEncoder
from openff.evaluator.datasets.datasets import PhysicalPropertyDataSet

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

def _get_raw_hash(physical_property) -> int:
    """
    Get the raw hash of a property based on its attributes.

    This method serializes the property attributes into a JSON string,
    sorts the keys, and computes a SHA-256 hash of the resulting string.
    The hash is then converted to an integer.

    Note: unlike `_get_hash`, this method does not truncate the hash value,
    so the number can be quite large.
    """

    type_ = type(physical_property)
    clsname = f"{type_.__module__}.{type_.__qualname__}"

    obj = {
        "type": clsname,
        "substance": physical_property.substance,
        "phase": physical_property.phase,
        "thermodynamic_state": physical_property.thermodynamic_state,
        "value": physical_property.value,
        "uncertainty": physical_property.uncertainty,
        "source": physical_property.source,
        "metadata": physical_property.metadata,
    }
    serialized = json.dumps(obj, sort_keys=True, cls=TypedJSONEncoder)
    return int(hashlib.sha256(serialized.encode("utf-8")).hexdigest(), 16)

def _get_hash(physical_property) -> int:
    """
    Returns a hash of the property based on attributes that are expected to
    have a meaningful value for the property. Hashes will change based on:

    - the property type
    - the value and uncertainty of the property
    - thermodynamic state
    - phase
    - substance
    - source
    - metadata

    The hash value will not depend on:
    - the id of the property (which is expected to be unique)
    - the gradients

    Note: as with the default __hash__ method in Python,
    the hash is truncated to the size of a Py_ssize_t, which is
    platform-dependent.

    Returns
    -------
    int
        The hash value of the property.
    """
    # hash() truncates the value returned from an object’s custom __hash__()
    # method to the size of a Py_ssize_t.
    # see https://docs.python.org/3/library/functions.html#hash
    # and https://github.com/python/cpython/blob/main/Include/cpython/pyhash.h#L8-L17
    
    raw_hash = _get_raw_hash(physical_property)

    # here we mimic the Python hash function for ease of comparison
    # and uses Mersenne primes for truncation
    if sys.hash_info.width == 64:
        mod = (1 << 61) - 1  # Mersenne prime for 64-bit hash
    else:
        mod = (1 << 31) - 1

    return raw_hash % mod


@click.command()
@click.option(
    "--input-file",
    "-i",
    default="output/initial-filtered.csv",
    help=(
        "The CSV file containing existing parsed ThermoML data. "
        "This should be a valid PhysicalPropertyDataSet CSV file."
    )
)
@click.option(
    "--output-file",
    "-o",
    default="intermediate/renamed-filtered.csv",
    help=(
        "The CSV file to save the renamed properties to. "
        "This will be a valid PhysicalPropertyDataSet CSV file."
    ),
)
def main(
    input_file: str,
    output_file: str,
):
    df = pd.read_csv(input_file, index_col=0)
    df["Id"] = df["Id"].astype(str)  # Ensure Id is string type
    dataset = PhysicalPropertyDataSet.from_pandas(df)
    logger.info(f"Loaded {len(dataset)} properties from {input_file}")

    SHORT_NAMES = {
        "Density": "dens",
        "EnthalpyOfMixing": "dhmix",
    }

    all_properties = {}

    for physical_property in tqdm.tqdm(dataset.properties):
        # Set the hash attribute to the computed hash
        hash_value = _get_hash(physical_property)
        clsname = SHORT_NAMES[type(physical_property).__name__]
        prop_id = f"{clsname}_{hash_value}"
        physical_property.id = prop_id
        all_properties[hash_value] = physical_property

    logger.info(f"Found {len(all_properties)} unique properties by hash")

    sorted_properties = [
        all_properties[hash_value] for hash_value in sorted(all_properties)
    ]
    new_dataset = PhysicalPropertyDataSet()
    new_dataset.add_properties(*sorted_properties)
    new_df = new_dataset.to_pandas()
    new_df.to_csv(output_file, index=False)
    logger.info(f"Saved renamed properties to {output_file}")



if __name__ == "__main__":
    main()
