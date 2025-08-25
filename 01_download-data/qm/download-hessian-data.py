"""
This script downloads Hessian data of optimzation results from the MolSSI QCArchive.
It filters out bad IDs, and saves the results in a specified output directory.

The output files are:
- `optimization_results.json`: Contains a combined OptimizationResultCollection with all existing optimization results.
- `hessian_results.json`: Contains a combined BasicResultCollection with all existing Hessian results
"""

import pathlib
import click
import tqdm
import logging

import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds

import qcportal as ptl
from openff.qcsubmit.utils.utils import portal_client_manager
from openff.qcsubmit.results import (
    BasicResultCollection,
    OptimizationResult,
    OptimizationResultCollection
)

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
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
        "This should be downloaded from download-qcdata-tables.py."
    ),
)
@click.option(
    "--output-directory",
    "-o",
    default="qm-data",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    help=(
        "Directory to write the output data. "
        "`optimization_results.json` contains a combined OptimizationResultCollection "
        "with all existing optimization results. "
        "`hessian_results.json` contains a combined BasicResultCollection "
        "with all existing Hessian results."
    ),
)
@click.option(
    "--exclude-bad-ids",
    "-x",
    default="../02_curate-data/bad-qcarchive_ids.dat",
    type=click.Path(exists=True, dir_okay=False),
    help=(
        "File containing a list of bad ids to exclude from the optimization data. "
        "This is a file with one id (int) per line."
    ),
)
def main(
    input_directory: str,
    output_directory: str,
    exclude_bad_ids: str = "../02_curate-data/bad-qcarchive_ids.dat"
):
    dataset = ds.dataset(input_directory)
    logger.info(f"Loaded {dataset.count_rows()} records")
    
    # filter out bad ids
    exclude_ids = set()
    if exclude_bad_ids:
        with open(exclude_bad_ids, "r") as f:
            exclude_ids = set(
                [int(line.strip()) for line in f.readlines()]
            )
        logger.info(f"Excluding {len(exclude_ids)} bad ids")

    subset = dataset.filter(
        ~pc.field("id").isin(exclude_ids)
    )
    logger.info(f"Filtered to {subset.count_rows()} records")

    df = subset.to_table(columns=["id", "inchi_key", "cmiles", "energy"]).to_pandas()
    assert len(df) > 0, "No records found after filtering"

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
    
    logger.info(f"Working with {optimization_collection.n_results} results")
    logger.info(f"Found {optimization_collection.n_molecules} molecules")
    logger.info(f"Wrote optimization results to {output_opt_file}")

    # convert to hessian
    with portal_client_manager(
        lambda x: ptl.PortalClient(x, cache_dir=".")
    ):
        hessian_set = optimization_collection.to_basic_result_collection(
            driver="hessian"
        )
    
    logger.info(f"Found {hessian_set.n_results} hessian calculations")
    logger.info(f"Found {hessian_set.n_molecules} hessian molecules")

    output_hessian_file = output_directory / "hessian_results.json"
    with open(output_hessian_file, "w") as f:
        f.write(hessian_set.json(indent=4))
    logger.info(f"Wrote hessian results to {output_hessian_file}")


if __name__ == "__main__":
    main()
