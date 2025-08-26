"""
Extract optimization data from a QCSubmit OptimizationResultCollection
and save it as a Parquet dataset for benchmarking. The output file is saved
with the same name as the input file, but with a .parquet extension,
in the specified output directory.

The script also writes an intermediate pickle file to speed up re-runs.

\b
Data is saved with the following schema:
- qcarchive_id (int): The QCArchive optimization ID.
- cmiles (str): The CMILES string for the molecule.
- mapped_smiles (str): The mapped SMILES string for the molecule.
- coordinates (List[float]): The flattened list of coordinates for the conformation in angstrom.
- energy (float): The energy of the conformation in kcal/mol.
- method (str): The method or forcefield used to generate the data, i.e., "qm".
- dataset (str): The name of the dataset, derived from the input file name
"""

import pathlib
import pickle
import click
import tqdm
from loguru import logger

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds

import qcelemental
import qcportal as ptl
from openff.units import unit
from openff.qcsubmit.utils.utils import portal_client_manager
from openff.qcsubmit.results import OptimizationResultCollection

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

hartree2kcalmol = qcelemental.constants.hartree2kcalmol

@click.command()
@click.option(
    "--input-file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="Path to the input file.",
)
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    default="qm-data",
    help="Path to the output directory.",
)
def main(
    input_file: str,
    output_directory: str,
    force: bool = False,
):
    output_directory = pathlib.Path(output_directory) / "qm"
    output_directory.mkdir(parents=True, exist_ok=True)

    # look for existing data to avoid duplicates
    existing_dataset = ds.dataset(output_directory)
    existing_qcarchive_ids = []
    n_files = 0
    if existing_dataset.count_rows():
        existing_qcarchive_ids = existing_dataset.to_table(
            columns=["qcarchive_id"]
        ).to_pydict()["qcarchive_id"]
        n_files = len(existing_dataset.files)

    input_file_name = pathlib.Path(input_file).stem
    pickle_file = f"{input_file_name}.pkl"

    collection = OptimizationResultCollection.parse_file(input_file)
    qcarchive_id_to_cmiles = {
        record.record_id: record.cmiles
        for record in collection.entries[QCFRACTAL_URL]
    }

    # check if pickle file already exists and reload
    if pathlib.Path(pickle_file).exists() and not force:
        logger.info(f"Loading records and molecules from {pickle_file}")
        with open(pickle_file, "rb") as f:
            records_and_molecules = pickle.load(f)
    else:
        with portal_client_manager(
            lambda x: ptl.PortalClient(x, cache_dir="../../03_fit-valence/02_curate-data/")
        ):
            records_and_molecules = collection.to_records()
        
        with open(pickle_file, "wb") as f:
            pickle.dump(records_and_molecules, f)
    
    # create new dataset
    entries = []
    for record, molecule in tqdm.tqdm(records_and_molecules):
        if record.id in existing_qcarchive_ids:
            continue
        assert len(molecule.conformers) == 1
        geometry = molecule.conformers[0].m_as(unit.angstrom)
        entry = {
            "qcarchive_id": record.id,
            "cmiles": qcarchive_id_to_cmiles[record.id],
            "mapped_smiles": molecule.to_smiles(mapped=True),
            "coordinates": geometry.flatten().tolist(),
            "energy": record.energies[-1] * hartree2kcalmol,
            "method": "qm",
            "dataset": input_file_name,
        }
        entries.append(entry)
    
    table = pa.Table.from_pylist(entries)
    table_file = output_directory / f"{input_file_name}.parquet"
    pq.write_table(table, table_file)
    logger.info(f"Wrote {len(entries)} entries to {table_file}")


if __name__ == "__main__":
    main()
