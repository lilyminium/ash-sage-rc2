"""
Select torsiondrive data for valence parameters.

Molecules are first sorted by size, filtered so that they can be parameterized by the force field,
and then filtered for chemical diversity. The final set of molecules is selected to aim for a
minimum number of records for each parameter.
"""
import functools
import json
import logging
import multiprocessing

from openff.toolkit import Molecule, ForceField

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import SimDivFilters
from rdkit.Chem import rdFingerprintGenerator
from openff.qcsubmit.results.results import TorsionDriveResult, TorsionDriveResultCollection

import qcportal as ptl
import tqdm
import click

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def select_by_size(
    df: pd.DataFrame,
    n_to_select: int = 1,
):
    cmiles_pool = df.cmiles.values
    torsion_ids = df.id.values
    rdmols = [
        Chem.MolFromSmiles(cmiles) for cmiles in cmiles_pool
    ]
    mws = [
        rdMolDescriptors.CalcExactMolWt(rdmol)
        for rdmol in rdmols
    ]
    sorted_indices = np.argsort(mws)
    selected_indices = sorted_indices[:n_to_select]
    selected_ids = [
        torsion_ids[i]
        for i in selected_indices
    ]
    return selected_ids


def select_by_chemical_diversity(
    df: pd.DataFrame,
    n_to_select: int = 1,
):
    if n_to_select == 1:
        return [df.id.values[0]]
    
    cmiles_pool = df.cmiles.values
    torsion_ids = df.id.values

    rdmols = [
        Chem.MolFromSmiles(cmiles) for cmiles in cmiles_pool
    ]
    mpfgen = rdFingerprintGenerator.GetMorganGenerator()
    fingerprints = [
        mpfgen.GetFingerprint(rdmol)
        for rdmol in rdmols
    ]
    mmp = SimDivFilters.MaxMinPicker()
    pick_size = min([n_to_select, len(rdmols)])
    picked_indices = list(
        mmp.LazyBitVectorPick(
            fingerprints,
            poolSize=len(rdmols),
            pickSize=pick_size,
        )
    )
    return [torsion_ids[i] for i in picked_indices]


def cmiles_to_inchi(cmiles: str) -> str:
    return Molecule.from_mapped_smiles(
        cmiles,
        allow_undefined_stereo=True
    ).to_inchikey(fixed_hydrogens=True)



def add_torsions_to_selected_dataset(
    torsion_id: int,
    overall_df: pd.DataFrame,
    parameter_counts: dict[str, int],
    selected_ids: set[int],
):
    # update counts of other parameters in central bond
    if torsion_id in selected_ids:
        return
    selected_ids.add(torsion_id)
    subdf = overall_df[
        overall_df.id == torsion_id
    ]
    assert len(subdf) == 1
    for row in subdf.all_parameter_ids.values:
        for pid in row:
            parameter_counts[pid] += 1
    

def remove_torsion_from_selected_dataset(
    torsion_id: int,
    overall_df: pd.DataFrame,
    parameter_counts: dict[str, int],
):
    subdf = overall_df[
        overall_df.id == torsion_id
    ]
    assert len(subdf) == 1
    for row in subdf.all_parameter_ids.values:
        for pid in row:
            parameter_counts[pid] -= 1


# @functools.cache
def can_parameterize_cmiles(cmiles: str, forcefield: ForceField) -> bool:
    try:
        mol = Molecule.from_mapped_smiles(
            cmiles,
            allow_undefined_stereo=True
        )
        mol.assign_partial_charges("zeros")
        forcefield.create_interchange(
            mol.to_topology(),
            charge_from_molecules=[mol]
        )
    except Exception as e:
        logger.info(e)
        return (False, cmiles)
    return (True, cmiles)

def filter_for_ff(
    cmiles_list: list[str],
    forcefield: ForceField,
    n_processes: int = 4
):
    """
    Filter the cmiles list for those that are compatible with the forcefield.
    """
    filtered_cmiles = []
    with multiprocessing.Pool(n_processes) as pool:
        results = list(
            tqdm.tqdm(
                pool.imap(
                    functools.partial(can_parameterize_cmiles, forcefield=forcefield),
                    cmiles_list,
                ),
                total=len(cmiles_list),
                desc="Filtering cmiles",
            )
        )
    filtered_cmiles = [
        cmiles
        for can_param, cmiles in results
        if can_param
    ]
    # for cmiles in tqdm.tqdm(cmiles_list, desc="Check FF parameterization"):
    #     if can_parameterize_cmiles(cmiles, forcefield):
    #         filtered_cmiles.append(cmiles)
    return filtered_cmiles


def filter_for_exclude_smarts(
    cmiles: list[str],
    patterns: list[str],
) -> list[str]:
    successful = []
    for cmi in cmiles:
        mol = Molecule.from_smiles(cmi, allow_undefined_stereo=True)
        if any(
            mol.chemical_environment_matches(pattern)
            for pattern in patterns
        ):
            continue
        else:
            successful.append(cmi)
    return successful


@click.command()
@click.option(
    "--input-directory",
    "-i",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="parameters/torsions",
    help="Directory containing the input tables."
)
@click.option(
    "--table-directory",
    "-t",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data/tables/torsiondrive",
    help="Directory containing the input tables."
)
@click.option(
    "--forcefield",
    "-ff",
    type=str,
    help="Path to the forcefield file.",
    default="openff_unconstrained-2.2.1.offxml",
)
@click.option(
    "--n-records",
    "-n",
    type=int,
    default=10,
    help="Number of records to select for each parameter.",
)
@click.option(
    "--output-count-file",
    "-oc",
    type=str,
    default="torsion-counts.json",
    help="Path to the output JSON file with counts.",
)
@click.option(
    "--output-file",
    "-o",
    type=str,
    default="torsiondrives.json",
    help="Path to the output JSON file.",
)
@click.option(
    "--exclude-qcarchive-file",
    "-xq",
    "exclude_qcarchive_files",
    type=str,
    multiple=True,
    default=["bad-qcarchive_ids.dat"],
    help="Path to the file containing bad qcarchive ids.",
)
@click.option(
    "--exclude-cmiles-file",
    "-xc",
    "exclude_cmiles_files",
    type=str,
    multiple=True,
    default=["failed-charge-cmiles.dat"],
    help="Path to the file containing bad cmiles.",
)
@click.option(
    "--exclude-smarts-file",
    "-xs",
    "exclude_smarts_files",
    type=str,
    multiple=True,
    default=[],
    help="Path to the file containing bad smarts.",
)
@click.option(
    "--exclude-dataset-file",
    "-xd",
    "exclude_dataset_files",
    type=str,
    multiple=True,
    default=[],
    help="Path to the file containing dataset entries to exclude, e.g. existing training data.",
)
@click.option(
    "--exclude-dataset-names",
    "-xdn",
    "exclude_dataset_names",
    type=str,
    multiple=True,
    default=[],
    help="Names of datasets to exclude from the selection, e.g. 'OpenFF Industry Benchmark Season 1 v1.2'.",
)
@click.option(
    "--n-processes",
    "-np",
    type=int,
    default=1,
    help="Number of processes to use for parallel processing.",
)
def main(
    input_directory: str = "parameters/torsions",
    table_directory: str = "data/tables/torsiondrive",
    forcefield: str = "openff_unconstrained-2.2.1.offxml",
    n_records: int = 10,
    output_count_file: str = "torsion-counts.json",
    output_file: str = "torsiondrives.json",
    exclude_qcarchive_files: list[str] = ["bad-qcarchive_ids.dat"],
    exclude_cmiles_files: list[str] = ["failed-charge-cmiles.dat"],
    exclude_smarts_files: list[str] = [],
    exclude_dataset_files: list[str] = [],
    exclude_dataset_names: list[str] = [],
    n_processes: int = 1,
):
    # load force field and parameters
    ff = ForceField(forcefield)
    all_parameter_ids = []
    handler = ff.get_parameter_handler("ProperTorsions")
    for parameter in handler.parameters:
        all_parameter_ids.append(parameter.id)

    logger.info(f"Loaded {len(all_parameter_ids)} torsion ids")
    logger.info(f"Looking for {n_records} each")

    # Load exclude file
    exclude_ids = set()
    if exclude_qcarchive_files:
        for exclude_file in exclude_qcarchive_files:
            with open(exclude_file, "r") as f:
                exclude_ids |= set([int(x.strip()) for x in f.readlines()])
        logger.info(f"Loaded {len(exclude_ids)} exclude ids from {exclude_file}")
    
    if exclude_dataset_files:
        for exclude_file in exclude_dataset_files:
            with open(exclude_file, "r") as f:
                contents = json.load(f)
            for entry_list in contents["entries"].values():
                for entry in entry_list:
                    exclude_ids.add(int(entry["record_id"]))
            logger.info(f"Loaded {len(exclude_ids)} exclude ids from {exclude_file}")

    exclude_cmiles = set()
    if exclude_cmiles_files:
        for exclude_file in exclude_cmiles_files:
            with open(exclude_file, "r") as f:
                exclude_cmiles |= set([x.strip() for x in f.readlines()])
        logger.info(f"Loaded {len(exclude_cmiles)} exclude cmiles from {exclude_file}")

    exclude_smarts = set()
    if exclude_smarts_files:
        for exclude_file in exclude_smarts_files:
            with open(exclude_file, "r") as f:
                exclude_smarts |= set([x.strip() for x in f.readlines()])
        logger.info(f"Loaded {len(exclude_smarts)} exclude smarts from {exclude_file}")

    dataset = ds.dataset(input_directory)
    df = dataset.to_table().to_pandas()
    logger.info(f"Loaded dataset of {len(df)} records")
    # Filter out exclude cmiles
    df = df[~df.cmiles.isin(exclude_cmiles)]
    logger.info(f"Filtered dataset to {len(df)} records after excluding bad cmiles")

    # Filter out exclude ids
    if exclude_ids:
        df = df[
            ~df.id.isin(exclude_ids)
        ]
        logger.info(f"Filtered dataset to {len(df)} records after excluding bad QCA IDs")

    # filter out cmiles that cannot be parameterized
    cmiles_list = df.cmiles.unique()
    cmiles_list = filter_for_exclude_smarts(
        cmiles_list,
        exclude_smarts,
    )
    filtered_cmiles = filter_for_ff(
        cmiles_list,
        forcefield=ff,
        n_processes=n_processes
    )
    # filtered_cmiles = cmiles_list
    logger.info(f"Filtered cmiles to {len(filtered_cmiles)} CMILES")
    df = df[
        df.cmiles.isin(filtered_cmiles)
    ]

    table_dataset = ds.dataset(table_directory)
    logger.info(f"Loading {table_dataset.count_rows()} records from {table_directory}")
    table_subset = table_dataset.filter(
        ~pc.field("id").isin(exclude_ids)
    )
    logger.info(f"Filtered to {table_subset.count_rows()} records after removing excluded ids")
    table_subset = table_subset.filter(
        ~pc.field("dataset_name").isin(exclude_dataset_names)
    )
    logger.info(f"Filtered to {table_subset.count_rows()} records after removing excluded dataset names")

    parameter_to_ids = {
        k: subdf.id.unique()
        for k, subdf in df.groupby(by="parameter_id")
    }
    
    parameter_counts = {
        k: 0
        for k in all_parameter_ids
    }

    selected_ids = set()
    passed_parameters = set()

    
    # first pass: pick all molecules that contain rare parameters
    for parameter_id, torsion_ids in tqdm.tqdm(
        parameter_to_ids.items(),
        desc="First pass"
    ):
        if len(torsion_ids) <= n_records:

            for torsion_id in torsion_ids:
                add_torsions_to_selected_dataset(
                    torsion_id,
                    df,
                    parameter_counts,
                    selected_ids
                )
            passed_parameters.add(parameter_id)
            logger.info(f"Selected {len(torsion_ids)} for {parameter_id}")
    
    # second pass: pick remaining torsions
    sorted_parameters_by_rarity = sorted(
        parameter_to_ids.items(),
        key=lambda x: len(x[1])
    )
    for parameter_id, _ in tqdm.tqdm(
        sorted_parameters_by_rarity,
        desc="Second pass"
    ):
        if parameter_id in passed_parameters:
            continue
        n_remaining = n_records - parameter_counts[parameter_id]
        if n_remaining <= 0:
            continue
        parameter_torsions = parameter_to_ids.get(parameter_id, [])
        parameter_torsions = set(parameter_torsions) - selected_ids

        if len(parameter_torsions) <= n_remaining:
            additional_ids = parameter_torsions
        else:
            subdf = df[
                df.id.isin(parameter_torsions)
            ]
            # first pick smaller molecules
            smallest_ids = select_by_size(
                subdf,
                n_to_select=max([20, n_records * 5]),
            )
            # then pick diverse molecules
            subdf = subdf[
                subdf.id.isin(smallest_ids)
            ]
            
            # sort subdf by size / index of smallest_ids
            subdf = subdf.set_index("id", drop=False)
            subdf = subdf.loc[
                list(map(int, smallest_ids)),
                :,
            ]
            additional_ids = select_by_chemical_diversity(
                subdf,
                n_to_select=n_remaining,
            )

        for torsion_id in additional_ids:
            add_torsions_to_selected_dataset(
                torsion_id,
                df,
                parameter_counts,
                selected_ids
            )
        logger.info(f"Selected {len(additional_ids)} for {parameter_id}")


    logger.info(f"Selected {len(selected_ids)} torsions")

    with open(output_count_file, "w") as f:
        json.dump(parameter_counts, f, indent=4)


    # logger.info all parameters with low counts
    logger.info("Parameters with 0 counts:")
    for parameter_id, count in parameter_counts.items():
        if not count:
            logger.info(parameter_id)
    logger.info("Parameters with 1-5 counts:")
    for parameter_id, count in parameter_counts.items():
        if count < 5 and count:
            logger.info(parameter_id)

    torsiondrive_results = []
    torsiondrive_df = df[df.id.isin(selected_ids)]
    assert len(torsiondrive_df) == len(selected_ids)

    for _, row in torsiondrive_df.iterrows():
        cmiles = row.cmiles
        torsiondrive_results.append(
            TorsionDriveResult(
                type="torsion",
                record_id=row["id"],
                cmiles=row["cmiles"],
                inchi_key=cmiles_to_inchi(cmiles),
            )
        )

    torsiondrive_collection = TorsionDriveResultCollection(
        type="TorsionDriveResultCollection",
        entries={
            QCFRACTAL_URL: torsiondrive_results
        }
    )

    with open(output_file, "w") as f:
        f.write(torsiondrive_collection.json(indent=4))
    
    logger.info(f"Wrote {len(torsiondrive_results)} records to {output_file}")



if __name__ == "__main__":
    main()
