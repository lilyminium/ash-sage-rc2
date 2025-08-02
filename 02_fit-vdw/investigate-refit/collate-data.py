import pathlib
import logging

import click
import tqdm
import pandas as pd

from openff.evaluator.client.client import RequestResult
from openff.evaluator.datasets.datasets import PhysicalPropertyDataSet


logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

@click.command()
@click.option(
    "--input-directory",
    "-i",
    default="../refit",
    help="The directory containing the results files from the refit",
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

    data_over_iterations = {}
    for prop in reference.properties:
        data_over_iterations[prop.id] = [
            type(prop).__name__, prop.value.m
        ]

    iter_cols = []
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
