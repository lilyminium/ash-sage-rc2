"""
This script benchmarks the optimization of molecules using OpenMM and a specified force field.
It reads a set of molecules from a QM dataset, optimizes them using the provided force field
and saves the optimized coordinates and energies.

It works off a single data directory, assumes the QM data is in a subdirectory called `qm`,
and saves the results in a subdirectory named after the force field used for optimization.
The output files are saved in Parquet format for efficient storage and retrieval.
Data is batch-processed with Dask and is written to work on a SLURM cluster.
\b
The format of the output files is:
- `qcarchive_id` (int): The QCArchive ID of the molecule.
- `cmiles` (str): The canonical SMILES of the molecule.
- `mapped_smiles` (str): The mapped SMILES of the molecule.
- `coordinates` (list[float]): The optimized coordinates of the molecule in Angstroms.
- `energy` (float): The optimized energy of the molecule in kcal/mol.
- `method` (str): The name of the force field used for optimization.
- `dataset` (str): The dataset the molecule belongs to.
"""

import pathlib
import logging
import typing
import click
import sys
import tqdm
import time

from click_option_group import optgroup

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

import openmm
import openmm.app
import openmm.unit
from openff.units import unit
from openff.toolkit import Molecule, ForceField

logger = logging.getLogger(__name__)
logging.basicConfig(
   level=logging.INFO,
   stream=sys.stdout,
   format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def optimize_single(
    row,
    forcefield: ForceField,
    ff_name: str,
) -> typing.Optional[dict[str, typing.Any]]:
    """
    Optimize a single molecule using the provided force field.
    Returns a dictionary with the optimized coordinates and energy.

    Parameters
    ----------
    row : dict[str, typing.Any]
        A dictionary containing the molecule data, including:
        - "qcarchive_id": The QCArchive ID of the molecule.
        - "cmiles": The canonical SMILES of the molecule.
        - "mapped_smiles": The mapped SMILES of the molecule.
        - "coordinates": The initial coordinates of the molecule.
        - "dataset": The dataset the molecule belongs to.
    forcefield : ForceField
        The force field to use for optimization.
    ff_name : str
        The name of the force field, used for labeling the output.

    Returns
    -------
    dict[str, typing.Any] or None
        A dictionary containing the optimized coordinates and energy, or None if optimization fails.
    """
    mol = Molecule.from_mapped_smiles(
        row["mapped_smiles"],
        allow_undefined_stereo=True
    )
    positions = np.array(row["coordinates"]).reshape((-1, 3))
    try:
        system = forcefield.create_interchange(mol.to_topology()).to_openmm(
            combine_nonbonded_forces=False,
        )
    except Exception as e:
        logger.warning(f"Skipping record {row['qcarchive_id']}: {e}")
        return None
    
    context = openmm.Context(
        system,
        openmm.VerletIntegrator(0.1 * openmm.unit.femtoseconds),
        openmm.Platform.getPlatformByName("Reference"),
    )

    context.setPositions(
        (positions * openmm.unit.angstrom).in_units_of(openmm.unit.nanometer),
    )
    openmm.LocalEnergyMinimizer.minimize(
        context=context,
        tolerance=10,
        maxIterations=0,
    )
    state = context.getState(getPositions=True, getEnergy=True)
    coordinates = state.getPositions(asNumpy=True).value_in_unit(openmm.unit.angstrom)
    energy = state.getPotentialEnergy().value_in_unit(openmm.unit.kilocalorie_per_mole)

    return {
        "qcarchive_id": row["qcarchive_id"],
        "cmiles": row["cmiles"],
        "mapped_smiles": row["mapped_smiles"],
        "coordinates": coordinates.flatten().tolist(),
        "energy": energy,
        "method": ff_name,
        "dataset": row["dataset"],
    }


def batch_optimize(
    qcarchive_ids: list[str],
    qm_directory: str,
    forcefield_path: str,
) -> list[dict[str, typing.Any]]:
    """
    Optimize a batch of molecules using the provided force field.

    Parameters
    ----------
    qcarchive_ids : list[str]
        A list of QCArchive IDs to process.
    qm_directory : str
        The directory containing the input QM data.
        This will get loaded as a dataset.
    forcefield_path : str
        The path to the force field file to use for optimization.

    Returns
    -------
    list[dict[str, typing.Any]]
        A list of dictionaries containing the optimized coordinates and energy for each molecule.
        Each dictionary contains:
        - "qcarchive_id": The QCArchive ID of the molecule.
        - "cmiles": The canonical SMILES of the molecule.
        - "mapped_smiles": The mapped SMILES of the molecule.
        - "coordinates": The optimized coordinates (A) of the molecule.
        - "energy": The optimized energy (kcal/mol) of the molecule.
        - "method": The name of the force field used for optimization.
        - "dataset": The dataset the molecule belongs to.

    """
    dataset = ds.dataset(qm_directory)
    subset = dataset.filter(
        pc.field("qcarchive_id").isin(qcarchive_ids)
    )
    forcefield = ForceField(forcefield_path)
    ff_name = pathlib.Path(forcefield_path).stem
    rows = subset.to_table().to_pylist()

    entries = []
    for row in tqdm.tqdm(rows):
        entry = optimize_single(row, forcefield, ff_name)
        if entry is not None:
            entries.append(entry)
    return entries



@click.command(help=__doc__)
@click.option(
    "--forcefield",
    "-ff",
    "forcefield",
    default="openff_unconstrained-2.2.1.offxml",
    help="Force field to use for labeling",
)
@click.option(
    "--data",
    "-d",
    "data_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data",
    help=(
        "Directory containing input QM data files. "
        "This should contain a 'qm' subdirectory with QM data files."
    ),
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

    # load input data
    data_directory = pathlib.Path(data_directory)
    input_directory = data_directory / "qm"
    input_dataset = ds.dataset(input_directory)
    logger.info(f"Loaded {input_dataset.count_rows()} rows from {input_directory}")
    
    # filter just to get QCArchive IDs to iterate over
    input_qcarchive_ids = input_dataset.to_table(
        columns=["qcarchive_id"]
    ).to_pydict()["qcarchive_id"]
    input_qcarchive_ids = set(input_qcarchive_ids)

    # work out output directory, which will be based on the FF name
    ff_name = pathlib.Path(forcefield).stem
    output_directory = data_directory / ff_name
    output_directory.mkdir(parents=True, exist_ok=True)
    output_dataset = ds.dataset(output_directory)
    n_files = 0
    # load any existing files to avoid re-processing
    if output_dataset.count_rows():
        existing_qcarchive_ids = output_dataset.to_table(
            columns=["qcarchive_id"]
        ).to_pydict()["qcarchive_id"]
        logger.info(f"Loaded {len(existing_qcarchive_ids)} rows from {output_directory}")
        input_qcarchive_ids -= set(existing_qcarchive_ids)
        logger.info(f"Filtered to {len(input_qcarchive_ids)} new rows to process")
        n_files  = len(output_dataset.files)

    input_qcarchive_ids = sorted(input_qcarchive_ids)

    # batch compute to save time
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
            batch_optimize,
            forcefield_path=forcefield,
            qm_directory=str(input_directory.resolve()),
        ))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Optimizing batches",
        ):
            entries = future.result()
            table = pa.Table.from_pylist(entries)
            table_file = output_directory / f"batch-{n_files:04d}.parquet"
            pq.write_table(table, table_file)
            logger.info(f"Wrote {len(entries)} entries to {table_file}")
            n_files += 1

    logger.info(f"{time.ctime()} - Finished batch optimization")
    elapsed_time = time.time() - start_time
    logger.info(f"Elapsed time: {elapsed_time / 60:.2f} min")
    logger.info("Done!")


if __name__ == "__main__":
    main()
