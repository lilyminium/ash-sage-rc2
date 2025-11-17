import collections
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

from yammbs.analysis import get_rmsd, get_tfd

logger.remove()
logger.add(sys.stdout)


def batch_get_rmsd(
    batch_mapped_smiles: list[str],
    qm_directory: str = None,
    ff_directory: str = None,
):
    qm_dataset = ds.dataset(qm_directory)
    qm_subset = qm_dataset.filter(
        pc.field("mapped_smiles").isin(batch_mapped_smiles)
    )
    qm_rows = qm_subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "coordinates", "energy"]
    ).to_pylist()
    qm_coordinates_by_qcarchive_id = {
        row["qcarchive_id"]: np.array(row["coordinates"]).reshape((-1, 3))
        for row in qm_rows
    }
    qcarchive_id_by_mapped_smiles = collections.defaultdict(list)
    for row in qm_rows:
        qcarchive_id_by_mapped_smiles[row["mapped_smiles"]].append(row["qcarchive_id"])

    qm_energy = {
        row["qcarchive_id"]: row["energy"]
        for row in qm_rows
    }
        

    ff_subset = ds.dataset(ff_directory).filter(
        pc.field("mapped_smiles").isin(batch_mapped_smiles)
    )
    ff_df = ff_subset.to_table(
        columns=["qcarchive_id", "cmiles", "mapped_smiles", "coordinates", "method", "energy"]
    ).to_pandas()
    
    # calculate all-to-all RMSD
    # and for each force field, take the *lowest* RMSD match to QM.
    # store the initial QCArchive ID and match QCArchive ID.
    rows = []
    for ff, ff_subdf in ff_df.groupby("method"):
        for mapped_smiles, smiles_df in tqdm.tqdm(ff_subdf.groupby("mapped_smiles")):
            qcarchive_ids = qcarchive_id_by_mapped_smiles[mapped_smiles]
            mol = Molecule.from_mapped_smiles(mapped_smiles, allow_undefined_stereo=True)
            inchi = mol.to_inchi(fixed_hydrogens=True)

            matched_qcarchive_id = -1
            current_rmsd = np.inf
            for _, row in smiles_df.iterrows():
                ff_coordinates = np.array(row["coordinates"]).reshape((-1, 3))
                assert row["qcarchive_id"] in qcarchive_ids
                for qca_id in qcarchive_ids:
                    qm_coordinates = qm_coordinates_by_qcarchive_id[qca_id]
                    rmsd = get_rmsd(mol, qm_coordinates, ff_coordinates)
                    print(f"Comparing FF {row['qcarchive_id']} to QM {qca_id}: RMSD = {rmsd:.4f} Å")
                    if rmsd < current_rmsd:
                        matched_qcarchive_id = qca_id
                        current_rmsd = rmsd
                        try:
                            # this sometimes errors out because RDKit looks for a bond that isn't there?
                            current_tfd = get_tfd(mol, qm_coordinates, ff_coordinates)
                        except Exception as e:
                            logger.info(f"Could not compute TFD for {inchi}: {e}")
                            current_tfd = np.nan
            
                entry = {
                    "mapped_smiles": mapped_smiles,
                    "cmiles": row["cmiles"],
                    "inchi": inchi,
                    "ff_qcarchive_id": row["qcarchive_id"],
                    "ff_energy": row["energy"],
                    "rmsd": current_rmsd,
                    "tfd": current_tfd,
                    "qm_qcarchive_id": matched_qcarchive_id,
                    "qm_energy": qm_energy[matched_qcarchive_id],
                    "method": ff,
                }
                rows.append(entry)

    return rows

@click.command()
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
@click.option(
    "--rmsd",
    "-r",
    "rmsd_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="rmsd",
    help="Directory to write RMSD files to",
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

    print(f"{time.ctime()} - Starting batch optimization")
    start_time = time.time()

    data_directory = pathlib.Path(data_directory)
    qm_directory = data_directory / "qm"
    qm_dataset = ds.dataset(qm_directory)
    print(f"Loaded {qm_dataset.count_rows()} rows from {qm_directory}")

    ff_name = pathlib.Path(forcefield).stem
    ff_directory = data_directory / ff_name
    ff_dataset = ds.dataset(ff_directory)
    print(f"Loaded {ff_dataset.count_rows()} rows from {ff_directory}")

    rmsd_directory = pathlib.Path(rmsd_directory) / ff_name
    rmsd_directory.mkdir(parents=True, exist_ok=True)
    
    input_mapped_smiles = ff_dataset.to_table(
        columns=["mapped_smiles"]
    ).to_pydict()["mapped_smiles"]
    input_mapped_smiles = set(input_mapped_smiles)
    print(f"Loaded {len(input_mapped_smiles)} mapped smiles to process")

    output_dataset = ds.dataset(rmsd_directory)
    n_files = 0
    if output_dataset.count_rows():
        existing_mapped_smiles = output_dataset.to_table(
            columns=["mapped_smiles"]
        ).to_pydict()["mapped_smiles"]
        print(f"Loaded {len(existing_mapped_smiles)} mapped smiles from {rmsd_directory}")
        input_mapped_smiles -= set(existing_mapped_smiles)
        print(f"Filtered to {len(input_mapped_smiles)} new mapped smiles to process")
        n_files  = len(output_dataset.files)

    # input_mapped_smiles = sorted(input_mapped_smiles)
    input_mapped_smiles = [
        "[H:42][C:23]1([C:22]([C:21]([C:20]([C:19]([C:26]([C:25]([C:24]1([H:44])[H:45])([H:46])[H:47])([H:48])[H:49])([H:35])[N:18]([H:34])[c:10]2[c:2]([c:3]([c:8]([c:16]([c:11]2[S:12](=[O:13])(=[O:15])[N:14]([H:32])[H:33])[F:17])[F:9])[S:4][C:5]([H:27])([H:28])[C:6]([H:29])([H:30])[O:7][H:31])[F:1])([H:36])[H:37])([H:38])[H:39])([H:40])[H:41])[H:43]"
    ]

    with batch_distributed(
        input_mapped_smiles,
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
            print(f"Wrote {len(entries)} entries to {table_file}")
            n_files += 1

    print(f"{time.ctime()} - Finished batch RMSD")
    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time / 60:.2f} min")
    print("Done!")


if __name__ == "__main__":
    main()
