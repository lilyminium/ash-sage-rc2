"""
This script calculates RMSD and TFD values for molecular conformations
as optimized by a force field from QM.
It uses a single data directory as input, assuming that QM data
is in a subdirectory named "qm" and force field data
is in a subdirectory named after the force field name.
It saves the results in a specified output directory.
Data is batch-processed with Dask and is written to work on a SLURM cluster.


The format of the output files is Parquet, with each file containing
the following columns:
- qcarchive_id (int): The QCArchive ID of the molecule.
- rmsd (float): The RMSD value between the QM and force field coordinates (A)
- tfd (float): The TFD value between the QM and force field coordinates (dimensionless)
- method (str): The FF name
"""

import pathlib
import sys
import typing
import click
import tqdm
import time

from loguru import logger
from click_option_group import optgroup

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq


from openff.toolkit import Molecule, ForceField

from yammbs.analysis import get_rmsd, get_tfd

logger.remove()
logger.add(sys.stdout)

def batch_get_rmsd(
    qcarchive_ids: list[int],
    qm_directory: str = None,
    ff_directory: str = None,
) -> list[dict]:
    """
    Calculate RMSD and TFD for a batch of QCArchive IDs.

    Parameters
    ----------
    qcarchive_ids : list[int]
        List of QCArchive IDs to process.
    qm_directory : str
        Input directory containing QM data files.
    ff_directory : str
        Input directory containing force field data files.
    
    Returns
    -------
    list[dict]
        List of dictionaries containing the results for each QCArchive ID.
        Each dictionary contains:
        - qcarchive_id (int): The QCArchive ID of the molecule.
        - rmsd (float): The RMSD value between the QM and force field coordinates (A)
        - tfd (float): The TFD value between the QM and force field coordinates
        - method (str): The FF name
    """
    qm_dataset = ds.dataset(qm_directory)
    qm_subset = qm_dataset.filter(
        pc.field("qcarchive_id").isin(qcarchive_ids)
    )
    qm_rows = qm_subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "coordinates"]
    ).to_pylist()
    qm_coordinates_by_qcarchive_id = {
        row["qcarchive_id"]: np.array(row["coordinates"]).reshape((-1, 3))
        for row in qm_rows
    }

    ff_subset = ds.dataset(ff_directory).filter(
        pc.field("qcarchive_id").isin(qcarchive_ids)
    )
    ff_rows = ff_subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "coordinates", "method"]
    ).to_pylist()
    

    entries = []
    for row in tqdm.tqdm(ff_rows):
        try:
            qm_coordinates = qm_coordinates_by_qcarchive_id[row["qcarchive_id"]]
        except KeyError:
            logger.warning(f"Skipping record {row['qcarchive_id']}: QM coordinates not found")
            continue
        mol = Molecule.from_mapped_smiles(
            row["mapped_smiles"],
            allow_undefined_stereo=True
        )
        ff_coordinates = np.array(row["coordinates"]).reshape((-1, 3))
        rmsd = get_rmsd(mol, qm_coordinates, ff_coordinates)
        try:
            tfd = get_tfd(mol, qm_coordinates, ff_coordinates)
        except Exception as e:
            logger.info(e)
            tfd = np.nan
        entry = {
            "qcarchive_id": row["qcarchive_id"],
            "rmsd": rmsd,
            "tfd": tfd,
            "method": row["method"],
        }
        entries.append(entry)
    return entries



@click.command()
@click.option(
    "--forcefield",
    "-ff",
    "forcefield",
    default="openff_unconstrained-2.2.1.offxml",
    help="The force field to select, which should have results in the data directory.",
)
@click.option(
    "--data",
    "-d",
    "data_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data",
    help=(
        "Directory containing data files. "
        "This directory should contain a subdirectory named 'qm' with QM data, "
        "and a subdirectory named after the force field with force field data."
    ),
)
@click.option(
    "--rmsd",
    "-r",
    "rmsd_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="rmsd",
    help="Directory to write RMSD/TFD output files to",
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-workers",
    help="The number of workers to distribute the labelling across. Use -1 to request "
    "one worker per batch.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.option(
    "--worker-type",
    help="The type of worker to distribute the labelling across.",
    type=click.Choice(["lsf", "local", "slurm"]),
    default="local",
    show_default=True,
)
@optgroup.option(
    "--batch-size",
    help="The number of molecules to processes at once on a particular worker.",
    type=int,
    default=500,
    show_default=True,
)
@optgroup.group("Cluster configuration", help="Options to configure cluster workers.")
@optgroup.option(
    "--memory",
    help="The amount of memory (GB) to request per queue worker.",
    type=int,
    default=3,
    show_default=True,
)
@optgroup.option(
    "--walltime",
    help="The maximum wall-clock hours to request per queue worker.",
    type=int,
    default=2,
    show_default=True,
)
@optgroup.option(
    "--queue",
    help="The SLURM queue to submit workers to.",
    type=str,
    default="cpuqueue",
    show_default=True,
)
@optgroup.option(
    "--conda-environment",
    help="The conda environment that SLURM workers should run using.",
    type=str,
)
def main(
    forcefield: str = "openff_unconstrained-2.2.1.offxml",
    data_directory: str = "data",
    rmsd_directory: str = "rmsd",
    worker_type: typing.Literal["slurm", "local"] = "local",
    queue: str = "free",
    conda_environment: str = "ib-dev",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed

    logger.info(f"{time.ctime()} - Starting batch optimization")
    start_time = time.time()

    # load input data QM
    data_directory = pathlib.Path(data_directory)
    qm_directory = data_directory / "qm"
    qm_dataset = ds.dataset(qm_directory)
    logger.info(f"Loaded {qm_dataset.count_rows()} rows from {qm_directory}")

    # load input FF data
    ff_name = pathlib.Path(forcefield).stem
    ff_directory = data_directory / ff_name
    ff_dataset = ds.dataset(ff_directory)
    logger.info(f"Loaded {ff_dataset.count_rows()} rows from {ff_directory}")

    rmsd_directory = pathlib.Path(rmsd_directory) / ff_name
    rmsd_directory.mkdir(parents=True, exist_ok=True)
    
    input_qcarchive_ids = ff_dataset.to_table(
        columns=["qcarchive_id"]
    ).to_pydict()["qcarchive_id"]
    input_qcarchive_ids = set(input_qcarchive_ids)
    logger.info(f"Loaded {len(input_qcarchive_ids)} rows to process")

    output_dataset = ds.dataset(rmsd_directory)
    n_files = 0
    # check if we have any existing files to avoid re-processing data
    if output_dataset.count_rows():
        existing_qcarchive_ids = output_dataset.to_table(
            columns=["qcarchive_id"]
        ).to_pydict()["qcarchive_id"]
        logger.info(f"Loaded {len(existing_qcarchive_ids)} rows from {rmsd_directory}")
        input_qcarchive_ids -= set(existing_qcarchive_ids)
        logger.info(f"Filtered to {len(input_qcarchive_ids)} new rows to process")
        n_files  = len(output_dataset.files)

    input_qcarchive_ids = sorted(input_qcarchive_ids)

    with batch_distributed(
        input_qcarchive_ids,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_get_rmsd,
            qm_directory=str(qm_directory.resolve()),
            ff_directory=str(ff_directory.resolve()),
        ))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Calculating RMSD",
        ):
            entries = future.result()
            table = pa.Table.from_pylist(entries)
            table_file = rmsd_directory / f"batch-{n_files:04d}.parquet"
            pq.write_table(table, table_file)
            logger.info(f"Wrote {len(entries)} entries to {table_file}")
            n_files += 1

    logger.info(f"{time.ctime()} - Finished batch RMSD")
    elapsed_time = time.time() - start_time
    logger.info(f"Elapsed time: {elapsed_time / 60:.2f} min")
    logger.info("Done!")


if __name__ == "__main__":
    main()
