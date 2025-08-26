import collections
import pathlib
import logging
import typing
import click
import tqdm
import json
import time

from click_option_group import optgroup

import numpy as np
import pyarrow as pa
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

import MDAnalysis as mda
from openff.units import unit
from openff.toolkit import Molecule

import seaborn as sns
from matplotlib import pyplot as plt

def batch_filter(
    mapped_smiles: list[str],
    pattern: str = None,
    topology_group: str = None,
) -> dict[str, list[str]]:

    pattern_matches = {}
    for smiles in tqdm.tqdm(mapped_smiles):
        mol = Molecule.from_mapped_smiles(
            smiles,
            allow_undefined_stereo=True
        )
        matches = mol.chemical_environment_matches(pattern)
        mol_matches = []
        if matches:
            for match in matches:
                # assume pattern was ordered sensibly
                if topology_group == "ImproperTorsions":
                    central = match[1]
                    others = sorted([match[0], match[2], match[3]])
                    match = [others[0], central, others[1], others[2]]
                elif match[0] > match[-1]:
                    match = match[::-1]
                mol_matches.append(tuple(match))
            pattern_matches[smiles] = mol_matches
    return pattern_matches

def compute_bond_difference(mm_value, qm_value):
    return mm_value - qm_value

def compute_angle_difference(mm_value, qm_value):
    difference = mm_value - qm_value
    if difference > 180:
        difference -= 360
    if difference < -180:
        difference += 360
    return difference


def get_unique_values(dataset, column_name: str = "mapped_smiles"):
    """
    Get the unique values in a column of a dataset in a memory efficient way.
    """
    scanner = dataset.scanner(columns=[column_name])
    unique_values = set()
    for batch in scanner.to_batches():
        col = batch.column(column_name)
        unique_values.update(col.to_pylist())
    return unique_values

                    

@click.command()
@click.option(
    "--input-directory",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="topology-values",
    help="Directory containing data files",
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="topology-comparisons",
    help="Directory to write the output files to",
)
@click.option(
    "--input-file",
    "-f",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default="comparison-patterns.json",
    help="JSON file containing the patterns to search for.",
)
@click.option(
    "--index",
    "-n",
    type=int,
    default=0,
    help="Index of the pattern to search for in the input file.",
    show_default=True,
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
    input_file: str,
    index: int,
    input_directory: str = "topology-values",
    output_directory: str = "topology-comparisons",
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

    print(f"{time.ctime()} - Starting batch filter")
    start_time = time.time()

    with open(input_file, "r") as f:
        contents = json.load(f)

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    kwargs = contents[index]
    pattern = kwargs["pattern"]
    topology_group = kwargs["topology_group"]
    output_stem = str(output_directory / kwargs["output_stem"])
    
    

    print(f"Searching for {pattern} and saving to {output_stem}")

    input_directory = pathlib.Path(input_directory)
    input_dataset = ds.dataset(input_directory)
    print(f"Loaded {input_dataset.count_rows()} rows from {input_directory}")

    unique_smiles = get_unique_values(input_dataset, column_name="mapped_smiles")
    print(f"Found {len(unique_smiles)} unique SMILES in {input_directory}")
    unique_smiles = sorted(unique_smiles)

    filtered_smiles_indices = {}

    with batch_distributed(
        unique_smiles,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_filter,
            pattern=pattern,
            topology_group=topology_group,
        ))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Filtering topology batches",
        ):
            filtered_smiles_indices.update(future.result())

    print(f"Found {len(filtered_smiles_indices)} matches for {pattern} in {topology_group}")

    initial_expression = (
        (pc.field("mapped_smiles").isin(list(filtered_smiles_indices)))
        & (pc.field("topology") == topology_group)
    )


    subset = input_dataset.filter(initial_expression)
    qm_subset = subset.filter(pc.field("method") == "qm")
    qm_rows = qm_subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "atom_indices", "value"]
    ).to_pylist()
    all_qm_values = collections.defaultdict(dict)
    for qm_row in tqdm.tqdm(qm_rows, desc="Getting QM values"):
        atom_indices = tuple(qm_row["atom_indices"])
        mapped_smiles = qm_row["mapped_smiles"]
        if atom_indices in filtered_smiles_indices[mapped_smiles]:
            all_qm_values[qm_row["qcarchive_id"]][atom_indices] = qm_row["value"]

    ff_subset = subset.filter(pc.field("method") != "qm")
    ff_rows = ff_subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "atom_indices", "value", "method"]
    ).to_pylist()

    diff_func = compute_bond_difference if topology_group == "Bonds" else compute_angle_difference

    all_output_entries = []
    for ff_row in tqdm.tqdm(ff_rows, desc="Getting FF values"):
        atom_indices = tuple(ff_row["atom_indices"])
        qcarchive_id = ff_row["qcarchive_id"]
        if atom_indices in all_qm_values[qcarchive_id]:
            qm_value = all_qm_values[qcarchive_id][atom_indices]
            entry = {
                "qcarchive_id": qcarchive_id,
                "mapped_smiles": ff_row["mapped_smiles"],
                "atom_indices": list(atom_indices),
                "mm_value": ff_row["value"],
                "qm_value": qm_value,
                "difference": diff_func(ff_row["value"], qm_value),
                "method": ff_row["method"].replace("_unconstrained", ""),
            }
            all_output_entries.append(entry)

    df = pd.DataFrame(all_output_entries)
    output_csv = output_stem + ".csv"
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} to {output_csv}")
    

    # plot
    ax = sns.boxplot(
        data=df,
        x="difference",
        y="method",
    )
    output_png = output_stem + ".png"
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    print(f"Saved plot to {output_png}")
    

    print(f"{time.ctime()} - Finished getting pattern values")
    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time / 60:.2f} min")
    print("Done!")


if __name__ == "__main__":
    main()
