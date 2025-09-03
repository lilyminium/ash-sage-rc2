"""
This script labels groups in a dataset and identifies groups
overrepresented in a high-error set, selected on a specified metric.
Groups are selected using Checkmol and a force field.

This saves two output files in the specified output directory.
The first is `labelled-properties.csv`.
This includes all columns of the input data, as well as
labels for the groups identified in each property.
\b
The second is `ratio_representation_{error_threshold}.csv`.
This includes the columns:
    - force field column (usually `FF`) (str): the force field
    - Group (str): the chemical group (either a Checkmol group or parameter ID)
    - n_dataset (int): the number of entries in the dataset
    - n_in_dataset (int): the number of entries in the dataset with this group present
    - n_high_error (int): the number of entries in the high-error set
    - n_in_high_error (int): the number of entries in the high-error set with this group present
    - fraction_in_dataset (float): n_in_dataset / n_dataset
    - fraction_in_high_error (float): n_in_high_error / n_high_error
    - ratio_representation (float): fraction_in_high_error / fraction_in_dataset

The `n_dataset` and `n_high_error` columns are expected to have
the same values for each group.
"""


from collections import defaultdict
import pathlib
import sys
import typing

import click
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from rdkit import Chem

from openff.toolkit import Molecule, ForceField

from loguru import logger

logger.remove()
logger.add(sys.stdout)

def sanitize_smiles(smiles: str) -> str:
    """Sanitize SMILES"""
    return Chem.MolToSmiles(
        Chem.MolFromSmiles(smiles)
    )


@click.command(help=__doc__)
@click.option(
    "--input",
    "-i",
    "input_paths",
    multiple=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=True, readable=True),
    help=(
        "Input paths. These can be files or directories. "
        "If a file, it is assumed to be a Pandas CSV. "
        "If a directory, it is assumed to contain a PyArrow dataset. "
        "These are loaded and concatenated into a single dataframe."
    )
)
@click.option(
    "--error-threshold",
    "-e",
    default=0.1,
    type=float,
    help=(
        "Error threshold for selecting high-error entries. "
        "If `error_metric='abs'`, this is in absolute units. "
        "If `error_metric='rel'`, this is a fraction."
    )
)
@click.option(
    "--output",
    "-o",
    "output_directory",
    default="output",
    type=click.Path(dir_okay=True, file_okay=False, writable=True),
    help=(
        "Output directory for the results. "
        "This directory will be created if it does not exist."
    )
)
@click.option(
    "--smiles-col",
    "-s",
    "smiles_columns",
    multiple=True,
    help=(
        "SMILES columns to use for grouping. "
        "These should be the names of the columns in the input dataframe."
    )
)
@click.option(
    "--error-metric",
    "-m",
    default="abs",
    type=click.Choice(["abs", "rel"], case_sensitive=False),
    help=(
        "Error metric to use for filtering. "
        "If 'abs', the error is in absolute units. "
        "If 'rel', the error is a fraction."
    )
)
@click.option(
    "--value-col",
    "-cv",
    type=str,
    default="value",
    help="Name of the column containing the values to analyze."
)
@click.option(
    "--method-col",
    "-cm",
    type=str,
    default="forcefield",
    help="Name of the column containing the method name (usually full name of FF)."
)
@click.option(
    "--force-field-col",
    "-cff",
    type=str,
    default="FF",
    help=(
        "Name of the column containing the human-friendly force field name, "
        "created with `name-and-forcefield`."
    )
)
@click.option(
    "--reference-value-col",
    "-crv",
    type=str,
    default="reference_value",
    help="Name of the column containing the reference values."
)
@click.option(
    "--id-col",
    "-cid",
    type=str,
    default="id",
    help="Name of the column containing the unique IDs to group results by."
)
@click.option(
    "--name-and-forcefield",
    "-nf",
    type=(str, str),
    multiple=True,
    help=(
        "Name and force field pairs to use for naming. "
        "These should be in the format -nf <name> <forcefield>"
    )
)
@click.option(
    "--restrict-to-shared-observables",
    "-r",
    is_flag=True,
    help=(
        "Restrict analysis to shared observables only. "
        "This means that only entries with all force fields present will be considered."
    )
)
@click.option(
    "--checkmol-labels-path",
    "-lcm",
    type=str,
    default="labels/checkmol",
    help="Path to the CheckMol labels dataset."
)
@click.option(
    "--forcefield-path",
    "-ff",
    type=str,
    default="../04_benchmark/forcefields/fb-fit-v1-single-mean-k100.offxml",
    help="Path to the force field file."
)
@click.option(
    "--forcefield-labels-path",
    "-lff",
    type=str,
    default="labels/forcefields/v1-k100",
    help="Path to the force field labels dataset."
)
def main(
    input_paths: list[str],
    error_threshold: float,
    output_directory: str = "output",
    smiles_columns: str = [],
    error_metric: typing.Literal["abs", "rel"] = "abs",
    value_col: str = "value",
    method_col: str = "forcefield",
    force_field_col: str = "FF",
    name_and_forcefield: list[tuple[str, str]] | None = None,
    reference_value_col: str | None = None,
    id_col: str | None = None,
    restrict_to_shared_observables: bool = False,
    checkmol_labels_path: str = "labels/checkmol",
    forcefield_path: str = "../04_benchmark/forcefields/fb-fit-v1-single-mean-k100.offxml",
    forcefield_labels_path: str = "labels/forcefields/v1-k100",
):
    dfs = []
    for input_path in input_paths:
        input_path = pathlib.Path(input_path)
        if input_path.is_dir():
            df_ = ds.dataset(input_path).to_table().to_pandas()
        else:
            df_ = pd.read_csv(input_path)
        logger.info(f"Read {len(df_)} entries from {input_path}")
        dfs.append(df_)
    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Read {len(df)} entries from {len(input_paths)} files")

    if not name_and_forcefield:
        name_and_forcefield = []

    # overwrite with human readable names
    STEM_TO_NAME = {
        v: k for k, v in name_and_forcefield
    }
    df[force_field_col] = [STEM_TO_NAME.get(v, v) for v in df[method_col]]
    if STEM_TO_NAME:
        df = df[df[method_col].isin(STEM_TO_NAME.keys())]
        logger.info(f"Filtered to {len(df)} entries with known forcefields")

    # compare like-to-like
    if restrict_to_shared_observables:
        assert id_col, "id_col must be specified if restrict_to_shared_observables is True"
        assert STEM_TO_NAME
        
        n_ff = len(STEM_TO_NAME)
        df_counts = df.groupby(id_col)[method_col].nunique()
        df_counts = df_counts[df_counts == n_ff]
        df = df[df[id_col].isin(df_counts.index)]
        logger.info(f"Filtered to {len(df)} entries with all {n_ff} forcefields")

    # get all smiles to check for
    all_smiles = set()
    for col in smiles_columns:
        all_smiles |= set(df[col].dropna().unique())
    logger.info(f"Found {len(all_smiles)} unique SMILES to label")

    # load labels
    expression = pc.field("smiles").isin(all_smiles)
    checkmol_labels = ds.dataset(checkmol_labels_path).filter(expression)
    checkmol_labels_df = checkmol_labels.to_table().to_pandas()
    logger.info(f"Loaded {len(checkmol_labels_df)} checkmol labels from {checkmol_labels_path}")
    assert all_smiles.issubset(set(checkmol_labels_df["smiles"].unique()))

    forcefield_labels = ds.dataset(forcefield_labels_path).filter(expression)
    forcefield_labels_df = forcefield_labels.to_table().to_pandas()
    logger.info(f"Loaded {len(forcefield_labels_df)} forcefield labels from {forcefield_labels_path}")
    assert all_smiles.issubset(set(forcefield_labels_df["smiles"].unique()))

    # assign to smiles
    CHECKMOL_GROUP_TO_SMILES = defaultdict(set)
    FORCEFIELD_GROUPS_TO_SMILES = defaultdict(set)
    ALL_CHECKMOL_GROUPS = set()
    ALL_FORCEFIELD_GROUPS = defaultdict(set)
    PARAMETER_IDS = {}
    PARAMETER_TYPES = ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW"]

    forcefield = ForceField(forcefield_path)

    for smiles, subdf in checkmol_labels_df.groupby("smiles"):
        groups = set(subdf["group"].unique())
        # SMILES_TO_GROUPS_CHECKMOL[smiles] = sorted(groups)
        ALL_CHECKMOL_GROUPS |= groups

        for group in groups:
            CHECKMOL_GROUP_TO_SMILES[group].add(smiles)

    for smiles, subdf in forcefield_labels_df.groupby("smiles"):
        for parameter_type, subsubdf in subdf.groupby("parameter_type"):
            parameter_handler = forcefield.get_parameter_handler(parameter_type)
            parameter_ids = [p.id for p in parameter_handler.parameters]
            PARAMETER_IDS[parameter_type] = parameter_ids

            groups = set(subsubdf["parameter_id"].unique())
            # SMILES_TO_GROUPS_FORCEFIELD[smiles][parameter_type] = sorted(groups)
            ALL_FORCEFIELD_GROUPS[parameter_type] |= groups

            for group in groups:
                FORCEFIELD_GROUPS_TO_SMILES[group].add(smiles)

    ALL_CHECKMOL_GROUPS -= {
        "Hydroxyl" # largely overlaps with Alcohol
    }
    ALL_CHECKMOL_GROUPS = sorted(ALL_CHECKMOL_GROUPS)
    logger.info(f"Found {len(ALL_CHECKMOL_GROUPS)} unique Checkmol groups")

    ALL_FORCEFIELD_GROUPS = {
        ptype: sorted(pids, key=lambda x: PARAMETER_IDS[ptype].index(x))
        for ptype, pids in ALL_FORCEFIELD_GROUPS.items()
    }
    for parameter_type, parameter_ids in ALL_FORCEFIELD_GROUPS.items():
        logger.info(f"Found {len(parameter_ids)} unique {parameter_type} forcefield parameters")

    # set column to presence of group in property
    for col in ALL_CHECKMOL_GROUPS:
        df[col] = df.apply(
            lambda row: any(
                row[smiles_col] in CHECKMOL_GROUP_TO_SMILES[col]
                for smiles_col in smiles_columns
            ),
            axis=1
        )
    df = df.copy() # stop fragmentation warnings
    for parameter_type in PARAMETER_TYPES:
        for col in ALL_FORCEFIELD_GROUPS[parameter_type]:
            df[col] = df.apply(
                lambda row: any(
                    row[smiles_col] in FORCEFIELD_GROUPS_TO_SMILES[col]
                    for smiles_col in smiles_columns
                ),
            axis=1
            )
        df = df.copy() # stop fragmentation warnings

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    output_file = output_directory / "labelled-properties.csv"
    df.to_csv(output_file, index=False)
    logger.info(f"Wrote {len(df)} entries to {output_file}")

    # calculate overrepresented groups

    # set up error calculation
    # calculate error
    if reference_value_col:
        df["error"] = df[value_col] - df[reference_value_col]
    else:
        # if not specified, could just be RMSD from QM metrics
        df["error"] = df[value_col]
    df["abs_error"] = df["error"].abs()

    N_DATASET = df[id_col].nunique()

    logger.info(f"Using {error_metric} error threshold {error_threshold}")

    ratio_rows = []
    for ffname, full_df in df.groupby(force_field_col):
        if not len(full_df) == N_DATASET:
            raise ValueError(
                f"Forcefield {ffname} does not have "
                f"all {N_DATASET} entries: {len(full_df)}"
            )

        # determine if we're selecting high error set by absolute error
        # or by top N entries
        error_metric = error_metric.lower()
        if error_metric == "abs":
            high_error_ids = full_df[full_df["abs_error"] > error_threshold][id_col].values
        elif error_metric == "rel":
            n_samples = int(error_threshold * N_DATASET)
            high_error_ids = full_df.sort_values(
                "abs_error", ascending=False
            )[id_col].values[:n_samples]

        high_error_df = full_df[full_df[id_col].isin(high_error_ids)]

        N_HIGH_ERROR = len(high_error_ids)
        logger.info(f"Found {N_HIGH_ERROR} high error entries for {ffname}")
        if not N_HIGH_ERROR:
            raise ValueError(f"No high error entries found for {ffname}")

        # calculate representation in high error set
        for col in ALL_CHECKMOL_GROUPS:
            N_IN_DATASET = int(full_df[col].sum())
            N_IN_HIGH_ERROR = int(high_error_df[col].sum())

            row = {
                force_field_col: ffname,
                "Group": col,
                "Type": "Checkmol",
                "n_dataset": N_DATASET,
                "n_in_dataset": N_IN_DATASET,
                "n_high_error": N_HIGH_ERROR,
                "n_in_high_error": N_IN_HIGH_ERROR,
            }
            ratio_rows.append(row)

        for parameter_type, parameter_ids in ALL_FORCEFIELD_GROUPS.items():
            for col in parameter_ids:
                N_IN_DATASET = int(full_df[col].sum())
                N_IN_HIGH_ERROR = int(high_error_df[col].sum())

                row = {
                    force_field_col: ffname,
                    "Group": col,
                    "Type": parameter_type,
                    "n_dataset": N_DATASET,
                    "n_in_dataset": N_IN_DATASET,
                    "n_high_error": N_HIGH_ERROR,
                    "n_in_high_error": N_IN_HIGH_ERROR,
                }
                ratio_rows.append(row)

    ratio_df = pd.DataFrame(ratio_rows)
    ratio_df["fraction_in_dataset"] = ratio_df.n_in_dataset / ratio_df.n_dataset
    ratio_df["fraction_in_high_error"] = ratio_df.n_in_high_error / ratio_df.n_high_error
    ratio_df["ratio_representation"] = ratio_df.fraction_in_high_error / ratio_df.fraction_in_dataset

    ratio_df = ratio_df.sort_values("ratio_representation", ascending=False)

    top_threshold = 1.5
    top_overrepresented = ratio_df[ratio_df["ratio_representation"] > top_threshold]
    logger.info(
        f"Found {len(top_overrepresented)} overrepresented groups "
        f"(ratio_representation > {top_threshold}) across {len(ratio_df[force_field_col].unique())} forcefields"
    )
    logger.info(f"Overrepresented groups: {', '.join(top_overrepresented.Group.values)}")

    ratio_csv = output_directory / f"ratio_representation_{error_threshold}.csv"
    ratio_df.to_csv(ratio_csv, index=False)
    logger.info(f"Wrote {len(ratio_df)} ratio representation entries to {ratio_csv}")


if __name__ == "__main__":
    main()
