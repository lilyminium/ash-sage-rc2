"""
This script collates training data from ForceBalance runs for vdW optimization.
It reads the training data from the `optimize.tmp` directory, which contains results from
intermediate iterations of the ForceBalance optimization.
The script aggregates the data per iteration and saves it to a CSV file.

Note: this hardcodes some expected paths.
For example, it assumes that the training set is in `targets/phys-prop/training-set.json`
and that the results are in `optimize.tmp/phys-prop/iter*/results.json`.
It also expects 15 iterations of results.
"""

import pathlib
import sys

import click
import tqdm
import pandas as pd

from loguru import logger
from openff.evaluator.client.client import RequestResult
from openff.evaluator.datasets.datasets import PhysicalPropertyDataSet

logger.remove()
logger.add(sys.stdout)

@click.command()
@click.option(
    "--input-directory",
    "-i",
    default="../refit",
    help=(
        "The directory containing the ForceBalance run for vdW optimization"
    ),
)
@click.option(
    "--output-file",
    "-o",
    default="output/training-per-iteration.csv",
    help="The output CSV file to save the training data per iteration",
)
def main(
    input_directory: str = "../refit",
    output_file: str = "output/training-per-iteration.csv",
):
    input_directory = pathlib.Path(input_directory)

    reference_dataset = input_directory / "targets/phys-prop/training-set.json"
    reference = PhysicalPropertyDataSet.from_json(str(reference_dataset))

    opt_directory = input_directory / "optimize.tmp"
    results_files = sorted(opt_directory.glob("phys-prop/iter*/results.json"))
    assert len(results_files) == 16, (
        f"Expected 16 iterations, found {len(results_files)}. "
        "Check the optimize.tmp directory for results."
    )

    # collect data for first two columns, Property type and reference value
    data_over_iterations: dict[str, list[str, float]] = {}
    for prop in reference.properties:
        data_over_iterations[prop.id] = [
            type(prop).__name__, prop.value.m  # we assume we know the units
        ]

    # iterate over results files and collect data for each iteration
    iter_cols: list[str] = []
    for result_file in tqdm.tqdm(results_files):
        iter_cols.append(result_file.parent.name)
        result = RequestResult.from_json(result_file)
        dataset = result.estimated_properties
        for prop in dataset.properties:
            data_over_iterations[prop.id].append(prop.value.m)

    df = pd.DataFrame.from_dict(
        data_over_iterations,
        orient="index",
        columns=["Property type", "Reference"] + iter_cols
    )
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(output_file)
    logger.info(f"Results saved to {output_file}")


if __name__ == "__main__":
    main()
