"""
Select optimization data for valence parameters.

Molecules are first sorted by size, filtered so that they can be parameterized by the force field,
and then filtered for chemical diversity. The final set of molecules is selected to aim for a
minimum number of records for each parameter.

The process is as follows:

- first, all CMILES are labelled with the parameter IDs of matching valence parameters (in a separate script in `label/`)
- this is loaded as input into the script
- for all parameter IDs where there are fewer than `n_records` CMILES, all CMILES are selected
  following a filtering process for any bad CMILES, QCA IDs, and SMARTS patterns
- for parameter IDs with more than `n_records` CMILES:
  - CMILES are first sorted and filtered by size to select the smallest molecules
  - they are filtered for bad QCA IDs, bad CMILES, and bad SMARTS patterns
  - they are selected for chemical diversity using RDKit's MaxMinPicker
- for all selected CMILES in the dataset, the n_conformer lowest energy conformers are selected
- the output is saved into the output_file as a `OptimizationResultCollection`.

"""
import collections
import functools
import logging
import multiprocessing
import json

from openff.toolkit import Molecule, ForceField

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import SimDivFilters
from rdkit.Chem import rdFingerprintGenerator
from openff.qcsubmit.results.results import OptimizationResult, OptimizationResultCollection

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
    cmiles_pool: list[str],
    n_to_select: int = 1,
) -> list[str]:
    """
    Select the smallest molecules from a pool of canonical SMILES, by molecular weight.

    Parameters
    ----------
    cmiles_pool : list[str]
        List of canonical SMILES strings to select from.
    n_to_select : int, optional
        Number of molecules to select, by default 1.
    
    Returns
    -------
    list[str]
        List of selected (smallest) canonical SMILES strings.
    """
    cmiles_pool = list(cmiles_pool)
    rdmols = [
        Chem.MolFromSmiles(cmiles) for cmiles in cmiles_pool
    ]
    mws = [
        rdMolDescriptors.CalcExactMolWt(rdmol)
        for rdmol in rdmols
    ]
    sorted_indices = np.argsort(mws)
    selected_indices = sorted_indices[:n_to_select]
    selected_cmiles = [
        cmiles_pool[i]
        for i in selected_indices
    ]
    return selected_cmiles


def select_by_chemical_diversity(
    cmiles_pool: list[str],
    n_to_select: int = 1,
):
    """
    Select a diverse set of molecules from a pool of canonical SMILES,
    using RDKit's MaxMinPicker and the Tanimoto distance between Morgan fingerprints.

    Parameters
    ----------
    cmiles_pool : list[str]
        List of canonical SMILES strings to select from.
    n_to_select : int, optional
        Number of molecules to select, by default 1.
    
    Returns
    -------
    list[str]
        List of selected (diverse) canonical SMILES strings.
    """
    if n_to_select == 1:
        return [cmiles_pool[0]]
    
    cmiles_pool = list(cmiles_pool)

    rdmols = [
        Chem.MolFromSmiles(cmiles) for cmiles in cmiles_pool
    ]
    mpfgen = rdFingerprintGenerator.GetMorganGenerator()
    fingerprints = [
        mpfgen.GetFingerprint(rdmol)
        for rdmol in rdmols
    ]
    # this by default will use Tanimoto distance
    mmp = SimDivFilters.MaxMinPicker()
    pick_size = min([n_to_select, len(rdmols)])
    picked_indices = list(
        mmp.LazyBitVectorPick(
            fingerprints,
            poolSize=len(rdmols),
            pickSize=pick_size,
        )
    )
    return [cmiles_pool[i] for i in picked_indices]


def cmiles_to_inchi(cmiles: str) -> str:
    """
    Convert a canonical SMILES string to an InChIKey, with fixed hydrogens.

    Parameters
    ----------
    cmiles : str
        Canonical SMILES string to convert.
    
    Returns
    -------
    str
        InChIKey representation of the molecule.
    """
    return Molecule.from_mapped_smiles(
        cmiles,
        allow_undefined_stereo=True
    ).to_inchikey(fixed_hydrogens=True)



def add_cmiles_to_selected_dataset(
    cmiles: str,
    overall_df: pd.DataFrame,
    parameter_counts: dict[str, int],
    selected_cmiles: set[str],
):
    """
    Add a cmiles to the selected dataset, updating the parameter counts.
    Note, this is a side-effecting helper function that modifies the `selected_cmiles` and `parameter_counts` dictionaries.

    Parameters
    ----------
    cmiles : str
        The canonical SMILES string to add.
    overall_df : pd.DataFrame
        The DataFrame containing the overall dataset.
    parameter_counts : dict[str, int]
        A dictionary mapping parameter IDs to their counts.
    selected_cmiles : set[str]
        A set of selected canonical SMILES strings.
    """
    # update counts of other parameters in selected cmiles
    if cmiles in selected_cmiles:
        return
    selected_cmiles.add(cmiles)
    cmiles_parameters = overall_df[
        overall_df.cmiles == cmiles
    ].parameter_id.unique()
    for cmiles_parameter in cmiles_parameters:
        parameter_counts[cmiles_parameter] += 1
    

def remove_cmiles_from_selected_dataset(
    cmiles: str,
    overall_df: pd.DataFrame,
    parameter_counts: dict[str, int],
):
    """
    Remove a cmiles from the selected dataset, updating the parameter counts.
    Note, this is a side-effecting helper function that modifies the `parameter_counts` dictionary.
    It does not modify the `selected_cmiles` set.

    Parameters
    ----------
    cmiles : str
        The canonical SMILES string to remove.
    overall_df : pd.DataFrame
        The DataFrame containing the overall dataset.
    parameter_counts : dict[str, int]
        A dictionary mapping parameter IDs to their counts.
    """
    cmiles_parameters = overall_df[
        overall_df.cmiles == cmiles
    ].parameter_id.unique()
    for cmiles_parameter in cmiles_parameters:
        parameter_counts[cmiles_parameter] -= 1
    
@functools.cache
def can_parameterize_cmiles(cmiles: str, forcefield: ForceField) -> tuple[bool, str]:
    """
    Check if a canonical SMILES can be parameterized by the given force field.

    Parameters
    ----------
    cmiles : str
        The canonical SMILES string to check.
    forcefield : ForceField
        The force field to use for parameterization.
    
    Returns
    -------
    tuple[bool, str]
        A tuple where the first element is a boolean indicating if the cmiles can be parameterized,
        and the second element is the cmiles string.
    """
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
        # logger.info(e)
        return (False, cmiles)
    return (True, cmiles)


def filter_for_ff(
    cmiles_list: list[str],
    forcefield: ForceField,
    n_processes: int = 4
) -> list[str]:
    """
    Filter the cmiles list for those that can be parameterized by the forcefield.

    Parameters
    ----------
    cmiles_list : list[str]
        List of canonical SMILES strings to filter.
    forcefield : ForceField
        The force field to use for filtering.
    n_processes : int, optional
        Number of processes to use for parallel processing, by default 4.

    Returns
    -------
    list[str]
        List of canonical SMILES strings that can be parameterized by the force field.
    """
    filtered_cmiles = []
    with multiprocessing.Pool(n_processes) as pool:
        results = list(
            pool.imap(
                functools.partial(can_parameterize_cmiles, forcefield=forcefield),
                cmiles_list,
            ),
        )
    filtered_cmiles = [
        cmiles
        for can_param, cmiles in results
        if can_param
    ]
    return filtered_cmiles


def filter_cmiles_for_good_qca_ids(
    table_subset: ds.Dataset,
    cmiles: list[str],
) -> set[str]:
    """
    Filter the cmiles list for those that have QCArchive IDs to include in the provided table subset.

    Parameters
    ----------
    table_subset : ds.Dataset
        The dataset to filter against, containing QCArchive IDs.
    cmiles : list[str]
        List of canonical SMILES strings to filter.

    Returns
    -------
    set[str]
        Set of canonical SMILES strings that have QCArchive IDs in the provided dataset.
    """
    if not len(cmiles):
        return set()
    
    subset = table_subset.filter(
        pc.field("cmiles").isin(cmiles)
    )
    if not subset.count_rows():
        return set()
    return set(
        subset.to_table(columns=["cmiles"]).to_pydict()["cmiles"]
    )


def filter_for_exclude_smarts(
    cmiles: list[str],
    patterns: list[str],
) -> list[str]:
    """
    Filter the cmiles list for those that do not match any of the provided SMARTS patterns
    
    Parameters
    ----------
    cmiles : list[str]
        List of canonical SMILES strings to filter.
    patterns : list[str]
        List of SMARTS patterns to exclude.

    Returns
    -------
    list[str]
        List of canonical SMILES strings that do not match any of the SMARTS patterns.
    """
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
    default="parameters/valence",
    help="Directory containing the input tables."
)
@click.option(
    "--table-directory",
    "-t",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data/tables/optimization",
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
    default=100,
    help="Number of records to select for each parameter.",
)
@click.option(
    "--n-conformers",
    "-nc",
    type=int,
    default=1,
    help="Number of conformers to select for each cmiles.",
)
@click.option(
    "--output-count-file",
    "-oc",
    type=str,
    default="valence-counts.json",
    help="Path to the output JSON file with counts.",
)
@click.option(
    "--output-file",
    "-o",
    type=str,
    default="optimizations.json",
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
    input_directory: str = "parameters/valence",
    table_directory: str = "data/tables/optimization",
    forcefield: str = "openff_unconstrained-2.2.1.offxml",
    n_records: int = 100,
    n_conformers: int = 1,
    output_count_file: str = "valence-counts.json",
    output_file: str = "optimizations.json",
    exclude_dataset_names: list[str] = [],
    exclude_qcarchive_files: list[str] = [],
    exclude_cmiles_files: list[str] = [],
    exclude_smarts_files: list[str] = [],
    exclude_dataset_files: list[str] = [],
    n_processes: int = 1
):
    # load force field and parameters
    ff = ForceField(forcefield)
    parameter_types = ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"]
    all_parameter_ids = []
    for parameter_type in parameter_types:
        handler = ff.get_parameter_handler(parameter_type)
        for parameter in handler.parameters:
            all_parameter_ids.append(parameter.id)

    logger.info(f"Loaded {len(all_parameter_ids)} parameter ids")
    logger.info(f"Looking for {n_records} each")

    # Load files with QCA IDs to exclude
    exclude_ids = set()
    if exclude_qcarchive_files:
        for exclude_file in exclude_qcarchive_files:
            with open(exclude_file, "r") as f:
                exclude_ids |= set([int(x.strip()) for x in f.readlines()])
        logger.info(f"Loaded {len(exclude_ids)} exclude ids from {exclude_file}")

    # Load files with datasets (.json) to exclude
    if exclude_dataset_files:
        for exclude_file in exclude_dataset_files:
            with open(exclude_file, "r") as f:
                contents = json.load(f)
            for entry_list in contents["entries"].values():
                for entry in entry_list:
                    exclude_ids.add(int(entry["record_id"]))
            logger.info(f"Loaded {len(exclude_ids)} exclude ids from {exclude_file}")

    # load CMILES to exclude
    exclude_cmiles = set()
    if exclude_cmiles_files:
        for exclude_file in exclude_cmiles_files:
            with open(exclude_file, "r") as f:
                exclude_cmiles |= set([x.strip() for x in f.readlines()])
        logger.info(f"Loaded {len(exclude_cmiles)} exclude cmiles from {exclude_file}")

    # load SMARTS to exclude
    exclude_smarts = set()
    if exclude_smarts_files:
        for exclude_file in exclude_smarts_files:
            with open(exclude_file, "r") as f:
                exclude_smarts |= set([x.strip() for x in f.readlines()])
        logger.info(f"Loaded {len(exclude_smarts)} exclude smarts from {exclude_file}")

    # actually load input labelled datasets matching parameter IDs to cmiles
    dataset = ds.dataset(input_directory)
    df = dataset.to_table().to_pandas()
    logger.info(f"Loaded {len(df)} CMILES from {input_directory}")

    # filter out CMILES to exclude
    df = df[~df.cmiles.isin(exclude_cmiles)]
    logger.info(f"Filtered to {len(df)} CMILES by removing excluded cmiles")

    # load originally downloaded data tables with QCA Ids, etc
    table_dataset = ds.dataset(table_directory)
    logger.info(f"Loading {table_dataset.count_rows()} records from {table_directory}")

    # exclude QCA IDs -- this is independent from `df` and CMILES
    # since some conformers might just be bad
    table_subset = table_dataset.filter(
        ~pc.field("id").isin(exclude_ids)
    )
    logger.info(f"Filtered to {table_subset.count_rows()} records after removing excluded ids")

    # exclude by dataset names
    table_subset = table_subset.filter(
        ~pc.field("dataset_name").isin(exclude_dataset_names)
    )
    logger.info(f"Filtered to {table_subset.count_rows()} records after removing excluded dataset names")
    
    
    # === begin filtering CMILES ===

    # map parameter ids to cmiles
    # this is a dictionary of all cmiles (value) matching a parameter id (key)
    parameter_to_cmiles: dict[str, set[str]] = {
        k: set(subdf.cmiles.unique())
        for k, subdf in df.groupby(by="parameter_id")
    }

    # set up a counting dictionary for how many cmiles we have per parameter
    parameter_counts: dict[str, int] = {
        k: 0
        for k in all_parameter_ids
    }

    # set up a set for selected dataset
    selected_cmiles: set[str] = set()
    passed_parameters: set[str] = set()

    # first pass: pick all molecules that contain rare parameters
    # by just taking all CMILES for parameters with fewer than n_records
    for parameter_id, parameter_cmiles in tqdm.tqdm(
        parameter_to_cmiles.items(),
        desc="First pass - selecting rare parameters"
    ):
        # select all CMILES where there are fewer than n_records
        if len(parameter_cmiles) <= n_records:
            # filter cmiles for good qca ids
            parameter_cmiles = filter_cmiles_for_good_qca_ids(
                table_subset,
                parameter_cmiles
            )
            # exclude any matching SMARTS
            parameter_cmiles = filter_for_exclude_smarts(
                parameter_cmiles,
                exclude_smarts,
            )
            # ensure FF can parameterize
            parameter_cmiles = filter_for_ff(
                parameter_cmiles,
                ff,
                n_processes=n_processes
            )
            logger.info(f"Selected {len(parameter_cmiles)} cmiles for {parameter_id}")
            for cmiles in parameter_cmiles:
                add_cmiles_to_selected_dataset(
                    cmiles,
                    df,
                    parameter_counts,
                    selected_cmiles
                )

            passed_parameters.add(parameter_id)

    logger.info(f"Processed {len(passed_parameters)} parameters in first pass with {len(selected_cmiles)} cmiles")
    
    # second pass: pick remaining cmiles
    sorted_parameters_by_rarity = sorted(
        parameter_to_cmiles.items(),
        key=lambda x: len(x[1])
    )
    for parameter_id, _ in tqdm.tqdm(
        sorted_parameters_by_rarity,
        desc="Second pass"
    ):
        if parameter_id in passed_parameters:
            continue
        # if we already have enough for this parameter, skip
        n_remaining = n_records - parameter_counts[parameter_id]
        if n_remaining <= 0:
            logger.info(f"Skipping {parameter_id} as we already have enough ({parameter_counts[parameter_id]})")
            continue
        parameter_cmiles = parameter_to_cmiles.get(parameter_id, [])
        parameter_cmiles = set(parameter_cmiles) - selected_cmiles

        n_original = len(parameter_cmiles)

        parameter_cmiles = filter_cmiles_for_good_qca_ids(
            table_subset,
            parameter_cmiles
        )
        logger.info(f"Filtered {n_original} to {len(parameter_cmiles)} cmiles after filtering for good qca ids for {parameter_id}")

        # take all remaining if that's all that's left
        if len(parameter_cmiles) <= n_remaining:
            additional_cmiles = parameter_cmiles
        else:

            logger.info(f"Filtering {len(parameter_cmiles)} cmiles for {parameter_id}")

            parameter_cmiles = select_by_size(
                parameter_cmiles,
                # choose a reasonable number to keep
                # do this first to save time
                n_to_select=max([100, n_records * 5]),
            )
            logger.info(f"Filtered to {len(parameter_cmiles)} cmiles after size filtering for {parameter_id}")

            parameter_cmiles = filter_for_exclude_smarts(
                parameter_cmiles,
                exclude_smarts,
            )
            logger.info(f"Filtered to {len(parameter_cmiles)} cmiles after excluding smarts for {parameter_id}")

            parameter_cmiles = filter_for_ff(
                parameter_cmiles,
                ff,
                n_processes=n_processes
            )
            logger.info(f"Filtered to {len(parameter_cmiles)} cmiles after checking for FF compatibility for {parameter_id}")


            if not parameter_cmiles:
                logger.info(f"No cmiles left for {parameter_id}")
                additional_cmiles = set()
            else:
                logger.info(f"Selecting {n_remaining} cmiles from {len(parameter_cmiles)} for {parameter_id}")
                # then pick diverse molecules
                additional_cmiles = select_by_chemical_diversity(
                    parameter_cmiles,
                    n_to_select=n_remaining,
                )
                logger.info(f"Filtered to {len(additional_cmiles)} cmiles after chemical diversity filtering for {parameter_id}")

        logger.info(f"Selected {len(additional_cmiles)} cmiles for {parameter_id}")
        for cmiles in additional_cmiles:
            add_cmiles_to_selected_dataset(
                cmiles,
                df,
                parameter_counts,
                selected_cmiles
            )

    logger.info(f"Selected {len(selected_cmiles)} cmiles total")

    with open(output_count_file, "w") as f:
        json.dump(parameter_counts, f, indent=4)
    logger.info(f"Wrote CMILES counts per parameter to {output_count_file}")


    # now pick lowest energy versions of each
    seen_ids = set()
    optimization_results = []
    # now we also count how many records there are per parameter
    cmiles_to_parameter: dict[str, set[str]] = collections.defaultdict(set)
    for parameter_id, cmiles_set in parameter_to_cmiles.items():
        for cmiles in cmiles_set:
            cmiles_to_parameter[cmiles].add(parameter_id)
    parameter_record_counts = collections.Counter()

    for cmiles in tqdm.tqdm(
        selected_cmiles,
        total=len(selected_cmiles),
        desc="Selecting lowest energy"
    ):
        
        expression = pc.field("cmiles") == cmiles
        subset: ds.Dataset = table_dataset.filter(expression)
        subset_df: pd.DataFrame = subset.to_table(
            columns=["id", "energy"]
        ).to_pandas().sort_values("energy")
        lowest_energy_ids: list[int] = subset_df["id"].values[:n_conformers]

        # create entries
        for lowest_energy_id in lowest_energy_ids:
            if lowest_energy_id in seen_ids:
                continue
            result = OptimizationResult(
                type="optimization",
                record_id=lowest_energy_id,
                cmiles=cmiles,
                inchi_key=cmiles_to_inchi(cmiles),
            )
            optimization_results.append(result)
            seen_ids.add(lowest_energy_id)

        # update counts
        for parameter_id in cmiles_to_parameter[cmiles]:
            parameter_record_counts[parameter_id] += len(lowest_energy_ids)

    # print all parameters with low counts
    no_smiles: list[str] = [
        parameter_id
        for parameter_id, count in parameter_counts.items()
        if count == 0
    ]
    logger.info(f"Parameters with 0 CMILES: {', '.join(no_smiles)}")

    rare_smiles: list[str] = [
        parameter_id
        for parameter_id, count in parameter_counts.items()
        if count < 5
    ]
    logger.info(f"Parameters with < 5 CMILES: {', '.join(rare_smiles)}")

    no_records: list[str] = [
        parameter_id
        for parameter_id, count in parameter_record_counts.items()
        if count == 0
    ]
    logger.info(f"Parameters with 0 records: {', '.join(no_records)}")

    rare_records: list[str] = [
        parameter_id
        for parameter_id, count in parameter_record_counts.items()
        if count < 10
    ]
    logger.info(f"Parameters with < 10 records: {', '.join(rare_records)}")

    # convert to QCSubmit
    optimization_collection = OptimizationResultCollection(
        type="OptimizationResultCollection",
        entries={
            QCFRACTAL_URL: optimization_results
        }
    )

    with open(output_file, "w") as f:
        f.write(optimization_collection.json(indent=4))
    
    logger.info(f"Wrote {len(optimization_results)} records to {output_file}")



if __name__ == "__main__":
    main()
