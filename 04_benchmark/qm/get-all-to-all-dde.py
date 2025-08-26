"""
Get ddE values from all-to-all RMSD data.

This script reads all-to-all RMSD data from a specified input directory,
filters the data based on a given RMSD threshold, and calculates ddE values
for each force field. The results are saved in a specified output directory,
where each force field has its own Parquet file containing the ddE values.
Force fields to include and their user-friendly names are provided as input arguments.

Output data follows the following schema:
- inchi (str): The InChI string of the molecule.
- ff_qcarchive_id (int): The QCArchive ID of the force field conformer
- qm_qcarchive_id (int): The QCArchive ID of the QM conformer.
- ddE (float): The ddE value (kcal/mol).
- ff_de (float): The force field energy difference (kcal/mol).
- qm_de (float): The QM energy difference (kcal/mol).
- method (str): The name of the force field.
- Force field (str): The user-friendly name of the force field.
- n_conformers (int): The number of conformers for the molecule.
"""

import pathlib
import sys
import click

from loguru import logger

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

logger.remove()
logger.add(sys.stdout)


def get_ddEs(
    df: pd.DataFrame,
) -> list[dict]:
    """
    Calculate ddE values for each force field in the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing columns "Force field", "inchi", "ff_qcarchive_id", "qm_qcarchive_id",
        "ff_energy", "qm_energy".

    Returns
    -------
    list[dict]
        List of dictionaries containing ddE values for each force field.
        Each dictionary contains:
        - inchi (str): The InChI string of the molecule.
        - ff_qcarchive_id (int): The QCArchive ID of the force field conform
        - qm_qcarchive_id (int): The QCArchive ID of the QM conformer.
        - ddE (float): The ddE value (kcal/mol).
        - ff_de (float): The force field energy difference (kcal/mol).
        - qm_de (float): The QM energy difference (kcal/mol).
        - method (str): The name of the force field.
        - n_conformers (int): The number of conformers for the molecule.
    """
    assert len(df["Force field"].unique()) == 1, "DataFrame should only contain one force field"
    ff = df["Force field"].unique()[0]
    method = df["method"].unique()[0]

    entries: list[dict] = []

    for inchi, subdf in df.groupby("inchi"):
        # find lowest energy conformer ID for each inchi
        lowest_qm_energy_idx: int = subdf["qm_energy"].idxmin()
        lowest_qm_energy: float = subdf.loc[lowest_qm_energy_idx, "qm_energy"]

        lowest_energy_qm_id: int = subdf.loc[lowest_qm_energy_idx, "qm_qcarchive_id"]
        lowest_ff_energy: float = subdf.loc[
            subdf.qm_qcarchive_id == lowest_energy_qm_id,
            "ff_energy"
        ].values.min()
        assert len(subdf[subdf.ff_energy == lowest_ff_energy]) == 1, \
            "There should be only one conformer with the lowest energy"
        lowest_energy_ff_id: int = subdf[subdf.ff_energy == lowest_ff_energy]["ff_qcarchive_id"].values[0]
        
        logger.info(
            f"Lowest energy conformer ID for {inchi}: {lowest_energy_ff_id} ({ff}) -> "
            f"{lowest_energy_qm_id} (QM)"
            f"with QM energy {lowest_qm_energy:.3f} kcal/mol"
        )

        for _, row in subdf.iterrows():
            # skip self-to-self comparisons
            if row["ff_qcarchive_id"] == lowest_energy_ff_id:
                continue
            mm_de = row["ff_energy"] - lowest_ff_energy
            qm_de = row["qm_energy"] - lowest_qm_energy
            entry = {
                "inchi": inchi,
                "ff_qcarchive_id": row["ff_qcarchive_id"],
                "qm_qcarchive_id": row["qm_qcarchive_id"],
                "ddE": mm_de - qm_de,
                "ff_de": mm_de,
                "qm_de": qm_de,
                "method": method,
                "Force field": ff,
                "n_conformers": subdf.shape[0]
            }
            entries.append(entry)
    return entries


@click.command()
@click.option(
    "--ff-name-and-stem",
    "-ff",
    "ff_name_and_stems",
    type=(str, str),
    multiple=True,
    help=(
        "Force fields to include and their stem names. "
        "The first argument should be the name of the force field in the plot (e.g. 'Sage 2.2.1'), "
        "the second argument should be the stem name of the force field (e.g., 'openff_unconstrained-2.2.1'). "
        "This option can be specified multiple times to include multiple force fields."
    )
)
@click.option(
    "--input-directory",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="all-to-all-rmsd",
    help="Directory to read all-to-all RMSD data from.",
)
@click.option(
    "--rmsd-threshold",
    "-t",
    "rmsd_threshold",
    type=float,
    default=0.3,
    help="Heavy atom RMSD (A) threshold for filtering entries.",
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="ddEs",
    help="Directory to write output files to.",
)
def main(
    ff_name_and_stems: list[tuple[str, str]],
    input_directory: str = "all-to-all-rmsd",
    rmsd_threshold: float = 0.3,
    output_directory: str = "ddEs",
):
    FF_STEM_TO_NAME = {
        stem: name for name, stem in ff_name_and_stems
    }
    stem_to_name_str = ", ".join([
        f"{stem} -> {name}"
        for stem, name in FF_STEM_TO_NAME.items()
    ])
    logger.info(f"Mapping {stem_to_name_str}")

    # load input data
    dataset = ds.dataset(input_directory)
    logger.info(f"Loaded {dataset.count_rows()} rows from {input_directory}")

    # filter dataset based on rmsd threshold
    filtered_dataset = dataset.filter(
        pc.field("rmsd") < rmsd_threshold
    )
    logger.info(f"Filtered dataset to {filtered_dataset.count_rows()} rows with RMSD < {rmsd_threshold} A")

    df = filtered_dataset.to_table(
        columns=[
            "mapped_smiles",
            "ff_qcarchive_id", "ff_energy",
            "qm_qcarchive_id", "qm_energy",
            "inchi", "method"
        ]
    ).to_pandas()
    df["Force field"] = df["method"].map(FF_STEM_TO_NAME)
    unique_ffs = df["Force field"].unique()
    logger.info(f"Found unique FFs: {', '.join(unique_ffs)}")

    # log counts
    counts = df.groupby("Force field").size().reset_index(name='count')
    count_str = ", ".join(
        f"{row['Force field']}: {row['count']}" for _, row in counts.iterrows()
    )
    logger.info(f"Counts of entries per force field -- {count_str}")

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    # get ddEs
    for ff_name, subdf in df.groupby("Force field"):
        entries = get_ddEs(subdf)
        table = pa.Table.from_pylist(entries)
        table_file = output_directory / f"{ff_name}.parquet"
        pq.write_table(table, table_file)
        logger.info(f"Wrote {len(entries)} ddE entries for {ff_name} to {table_file}")


if __name__ == "__main__":
    main()
