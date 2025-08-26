"""
This script benchmarks the optimization of torsiondrives using OpenMM and a specified force field.
It reads a set of molecules from a QM dataset, optimizes them using the provided force field
and saves the optimized coordinates and energies.

It works off a single data directory, assumes the QM data is in a subdirectory called `qm`,
and saves the results in a subdirectory named after the force field used for optimization.
The output files are saved in Parquet format for efficient storage and retrieval.
Data is batch-processed with Dask and is written to work on a SLURM cluster.
\b
The format of the output files is:
- `torsiondrive_id` (int): The ID of the torsion drive.
- `cmiles` (str): The canonical SMILES of the molecule.
- `mapped_smiles` (str): The mapped SMILES of the molecule.
- `dihedral_indices` (list[int]): The indices of the atoms involved in the dihedral.
- `angle` (float): The dihedral angle in degrees.
- `coordinates` (list[float]): The optimized coordinates of the molecule in Angstrom.
- `energy` (float): The optimized energy of the molecule in kcal/mol.
- `method` (str): The name of the force field used for optimization.
- `dataset` (str): The name of the dataset the molecule belongs to.
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

import openmm
import openmm.app
import openmm.unit
from openff.units import unit
from openff.toolkit import Molecule, ForceField

from openff.interchange.operations.minimize import _DEFAULT_ENERGY_MINIMIZATION_TOLERANCE

logger.remove()
logger.add(sys.stdout)


def optimize_single(
    row,
    forcefield: ForceField,
    ff_name: str,
):
    """
    Optimize a single row of data using the provided force field.
    The dihedral atoms are constrained; other atoms have a 1 kcal/mol/Å² restraint applied.

    Parameters
    ----------
    row : dict
        A dictionary containing the data for a single torsion drive.
    forcefield : ForceField
        The force field to use for optimization.
    ff_name : str
        The name of the force field, used for labeling.

    Returns
    -------
    dict[str, typing.Any]
        A dictionary containing the optimized data, including coordinates and energy.
    """
    restrain_k = 1.0

    mol = Molecule.from_mapped_smiles(
        row["mapped_smiles"],
        allow_undefined_stereo=True
    )
    positions = np.array(row["coordinates"]).reshape((-1, 3))
    mol.add_conformer(positions * openmm.unit.angstrom)

    interchange = forcefield.create_interchange(mol.to_topology())
    additional_forces = []
    dihedral_indices = list(map(int, row["dihedral_indices"]))

    # apply 1kcal/mol/Å² restraints to all atoms except the dihedral indices
    if True:
        restraint_force = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint_force.addGlobalParameter(
            "k",
            restrain_k * openmm.unit.kilocalorie_per_mole / openmm.unit.angstrom**2,
        )
        for parameter in ("x0", "y0", "z0"):
            restraint_force.addPerParticleParameter(parameter)


        positions = interchange.positions.to("nanometer")

        for atom_index in range(mol.n_atoms):
            if atom_index in dihedral_indices:
                continue

            particle_index = restraint_force.addParticle(atom_index)
            restraint_force.setParticleParameters(
                particle_index,
                atom_index,
                [x.to_openmm() for x in positions[atom_index]],
            )
        additional_forces = [restraint_force]

    simulation = interchange.to_openmm_simulation(
        openmm.LangevinMiddleIntegrator(
            293.15 * openmm.unit.kelvin,
            1.0 / openmm.unit.picosecond,
            2.0 * openmm.unit.femtosecond,
        ),
        combine_nonbonded_forces=True,
        additional_forces=additional_forces,
    )
    simulation.context.computeVirtualSites()

    # set the dihedral indices to zero mass
    for index in dihedral_indices:
        simulation.system.setParticleMass(index, 0.0)

    simulation.minimizeEnergy(
        tolerance=_DEFAULT_ENERGY_MINIMIZATION_TOLERANCE.to_openmm(),
        maxIterations=10_000,
    )

    state = simulation.context.getState(getPositions=True, getEnergy=True)
    coordinates = state.getPositions(asNumpy=True).value_in_unit(openmm.unit.angstrom)
    energy = state.getPotentialEnergy().value_in_unit(openmm.unit.kilocalorie_per_mole)

    return {
        "torsiondrive_id": row["torsiondrive_id"],
        "cmiles": row["cmiles"],
        "mapped_smiles": row["mapped_smiles"],
        "dihedral_indices": row["dihedral_indices"],
        "angle": row["angle"],
        "coordinates": coordinates.flatten().tolist(),
        "energy": energy,
        "method": ff_name,
        "dataset": row["dataset"],
    }


def batch_optimize(
    qcarchive_ids: list[str],
    qm_directory: str,
    forcefield_path: str,
):
    """
    Batch optimize a list of QCArchive IDs using the specified force field.

    Parameters
    ----------
    qcarchive_ids : list[str]
        A list of QCArchive IDs to process.
    qm_directory : str
        The directory containing the QM data files.
    forcefield_path : str
        The path to the force field file to use for optimization.

    Returns
    -------
    list[dict[str, typing.Any]]
        A list of dictionaries containing the optimized data for each QCArchive ID.
        Each dictionary contains the following keys:
            - "torsiondrive_id": The QCArchive ID.
            - "cmiles": The canonical SMILES representation of the molecule.
            - "mapped_smiles": The mapped SMILES representation of the molecule.
            - "dihedral_indices": The indices of the dihedral atoms.
            - "angle": The dihedral angle.
            - "coordinates": The optimized coordinates of the molecule.
            - "energy": The optimized energy of the molecule.
            - "method": The name of the force field used for optimization.
            - "dataset": The name of the dataset the molecule belongs to.
    """
    dataset = ds.dataset(qm_directory)
    subset = dataset.filter(
        pc.field("torsiondrive_id").isin(qcarchive_ids)
    )
    forcefield = ForceField(forcefield_path)
    ff_name = pathlib.Path(forcefield_path).stem
    rows = subset.to_table().to_pylist()

    entries = []
    for row in tqdm.tqdm(rows):
        try:
            entry = optimize_single(row, forcefield, ff_name)
        except:
            continue
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
    help="Directory containing data files",
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
    
    # filter just to get torsiondrive IDs to iterate over
    input_qcarchive_ids = input_dataset.to_table(
        columns=["torsiondrive_id"]
    ).to_pydict()["torsiondrive_id"]
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
            columns=["torsiondrive_id"]
        ).to_pydict()["torsiondrive_id"]
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
