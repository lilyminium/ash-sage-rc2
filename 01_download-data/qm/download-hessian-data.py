import pathlib
import click
import tqdm

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds

import qcportal as ptl
from openff.qcsubmit.utils.utils import portal_client_manager
from openff.qcsubmit.results import (
    BasicResultCollection,
    OptimizationResult,
    OptimizationResultCollection
)

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"


@click.command()
@click.option(
    "--input-directory",
    "-i",
    default="../02_curate-data/data/tables/optimization",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help=(
        "Directory containing the input optimization data. "
    ),
)
@click.option(
    "--output-directory",
    "-o",
    default="qm-data",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    help=(
        "Directory to write the output optimization data. "
    ),
)
@click.option(
    "--exclude-bad-ids",
    "-x",
    default="../02_curate-data/bad-qcarchive_ids.dat",
    type=click.Path(exists=True, dir_okay=False),
    help=(
        "File containing a list of bad ids to exclude from the optimization data. "
        "This is a file with one id per line."
    ),
)
def main(
    input_directory: str,
    output_directory: str,
    exclude_bad_ids: str = "../02_curate-data/bad-qcarchive_ids.dat"
):
    dataset = ds.dataset(input_directory)
    print(f"Loaded {dataset.count_rows()} records")
    
    # filter out bad ids
    exclude_ids = set()
    if exclude_bad_ids:
        with open(exclude_bad_ids, "r") as f:
            exclude_ids = set(
                [int(line.strip()) for line in f.readlines()]
            )
        print(f"Excluding {len(exclude_ids)} bad ids")

    subset = dataset.filter(
        ~pc.field("id").isin(exclude_ids)
    )
    print(f"Filtered to {subset.count_rows()} records")

    df = subset.to_table(columns=["id", "inchi_key", "cmiles", "energy"]).to_pandas()

    # get lowest energy ids per inchi_key
    entries = []
    for inchi_key, subdf in tqdm.tqdm(
        df.groupby("inchi_key"),
        desc="Finding lowest energy id",
    ):
        # get the lowest energy id
        index = np.argmin(subdf.energy.values)
        qca_id = subdf.id.values[index]
        cmiles = subdf.cmiles.values[index]
        result = OptimizationResult(
            type="optimization",
            record_id=qca_id,
            cmiles=cmiles,
            inchi_key=inchi_key
        )
        entries.append(result)

    optimization_collection = OptimizationResultCollection(
        type="OptimizationResultCollection",
        entries={
            QCFRACTAL_URL: entries
        }
    )

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    output_opt_file = output_directory / "optimization_results.json"

    with open(output_opt_file, "w") as f:
        f.write(optimization_collection.json(indent=4))
    
    print(f"Working with {optimization_collection.n_results} results")
    print(f"Found {optimization_collection.n_molecules} molecules")
    print(f"Wrote optimization results to {output_opt_file}")

    # convert to hessian
    with portal_client_manager(
        lambda x: ptl.PortalClient(x, cache_dir="../02_curate-data")
    ):
        hessian_set = optimization_collection.to_basic_result_collection(
            driver="hessian"
        )
    
    print(f"Found {hessian_set.n_results} hessian calculations")
    print(f"Found {hessian_set.n_molecules} hessian molecules")

    output_hessian_file = output_directory / "hessian_results.json"
    with open(output_hessian_file, "w") as f:
        f.write(hessian_set.json(indent=4))
    print(f"Wrote hessian results to {output_hessian_file}")


if __name__ == "__main__":
    main()
