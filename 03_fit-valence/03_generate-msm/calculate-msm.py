from collections import defaultdict
import logging
import os
import json
import pathlib
import pickle
import multiprocessing
import typing

import click
import numpy as np
from openff.toolkit import Molecule
from openff.units import unit

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds

import qcportal as ptl
from openff.qcsubmit.utils.utils import portal_client_manager

import tqdm
# suppress stereochemistry warnings
logging.getLogger("openff").setLevel(logging.ERROR)

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"


if typing.TYPE_CHECKING:
    from openff.toolkit import Molecule, ForceField
    from qcportal.models import OptimizationRecord, ResultRecord

from openff.qcsubmit.results import BasicResultCollection

def calculate_parameters_single(
    basic_result,
    qc_record: "ResultRecord",
    molecule: "Molecule",
) -> typing.Dict[str, typing.Dict[str, typing.List[unit.Quantity]]]:
    """
    Calculate the modified seminario parameters for the given input molecule
    and store them by OFF SMIRKS.
    """
    from qubekit.molecules import Ligand
    from qubekit.bonded.mod_seminario import ModSeminario

    mod_sem = ModSeminario()

    # create the qube molecule, this should be in the same order as the off_mol
    qube_mol = Ligand.from_rdkit(molecule.to_rdkit(), name="offmol")
    qube_mol.hessian = qc_record.return_result

    # calculate the modified seminario parameters and store in the molecule
    try:
        qube_mol = mod_sem.run(qube_mol)
    except Exception as e:
        print(f"Failed to calculate parameters for {basic_result.record_id}: {e}")
        return []

    entries = []
    for bond_parameter in qube_mol.BondForce:
        indices: tuple[int, ...] = bond_parameter.atoms
        if indices[-1] < indices[0]:
            indices = tuple(reversed(indices))

        entry = {
            "id": basic_result.record_id,
            "cmiles": basic_result.cmiles,
            "inchi_key": basic_result.inchi_key,
            "parameter_type": "Bonds",
            "indices": list(indices),
            "eq": bond_parameter.length, # nm
            "force_constant": bond_parameter.k, # kJ/nm2
        }
        entries.append(entry)
    
    for angle_parameter in qube_mol.AngleForce:
        indices: tuple[int, ...] = angle_parameter.atoms
        if indices[-1] < indices[0]:
            indices = tuple(reversed(indices))

        entry = {
            "id": basic_result.record_id,
            "cmiles": basic_result.cmiles,
            "inchi_key": basic_result.inchi_key,
            "parameter_type": "Angles",
            "indices": list(indices),
            "eq": angle_parameter.angle, # rad
            "force_constant": angle_parameter.k, # kJ/rad2
        }
        entries.append(entry)

    return entries



@click.command()
@click.option(
    "--input-dataset",
    "-i",
    default="qm-data/hessian-data.json",
    type=click.Path(exists=True, dir_okay=False),
    help=(
        "Path to the input dataset file. "
        "This is a JSON file containing the hessian data."
    ),
)
@click.option(
    "--output-directory",
    "-o",
    default="msm-data",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    help=(
        "Directory to write the output data. "
        "This will be a directory containing parquet files."
    ),
)
def main(
    input_dataset: str,
    output_directory: str,
):
    hessian_set = BasicResultCollection.parse_file(input_dataset)
    print(f"Found {hessian_set.n_results} hessian calculations")
    print(f"Found {hessian_set.n_molecules} hessian molecules")

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    # load existing data
    existing_data = ds.dataset(output_directory)
    n_files = 0
    if existing_data.count_rows():
        n_files = len(existing_data.files)
        found = set(existing_data.to_table(columns=["id"]).to_pydict()["id"])
        print(f"Found {len(found)} existing records")

        filtered_entries = [
            entry
            for entry in hessian_set.entries[QCFRACTAL_URL]
            if entry.record_id not in found
        ]
        hessian_set.entries[QCFRACTAL_URL] = filtered_entries
        print(f"Filtered to {len(filtered_entries)} records")

    if pathlib.Path("records_and_molecules.pkl").is_file():
        # reload to save on time

        with open("records_and_molecules.pkl", "rb") as f:
            records_and_molecules = pickle.load(f)
    else:

        with portal_client_manager(
            lambda x: ptl.PortalClient(x, cache_dir="../02_curate-data")
        ):
            records_and_molecules = list(hessian_set.to_records())

        with open("records_and_molecules.pkl", "wb") as f:
            pickle.dump(records_and_molecules, f)

    print(f"Loaded {len(records_and_molecules)} records and molecules")

    # filter for complete
    records_and_molecules = [
        (record, molecule)
        for record, molecule in records_and_molecules
        if record.status.upper() == "COMPLETE"
    ]
    print(f"Filtered to {len(records_and_molecules)} complete records and molecules")

    # match cmiles to records and molecules
    id_to_entry = {
        entry.record_id: entry
        for entry in hessian_set.entries[QCFRACTAL_URL]
    }
    entries_records_molecules = []
    for record, molecule in records_and_molecules:
        entry = id_to_entry[record.id]
        # try:
        #     molecule = Molecule.from_rdkit(molecule.to_rdkit(), allow_undefined_stereo=True)
        # except Exception as e:
        #     print(f"Failed to load molecule {record.id}: {e}")
        #     continue
        entries_records_molecules.append((entry, record, molecule))

    print(f"Working with {len(entries_records_molecules)} records and molecules")

    # save it in batches
    batch_size = 1000
    n_total = len(entries_records_molecules)
    n_batches = n_total // batch_size + 1
    file_number = n_files
    for i in tqdm.tqdm(
        range(0, n_total, batch_size),
        desc="Processing batches",
        total=n_batches
    ):
        batch = entries_records_molecules[i:i + batch_size]
        # calculate the parameters
        results = []
        for entry, record, molecule in tqdm.tqdm(batch):
            result = calculate_parameters_single(entry, record, molecule)
            results.extend(result)

        table = pa.Table.from_pylist(results)

        batch_file = output_directory / f"batch-{file_number:04d}.parquet"
        pq.write_table(table, batch_file)
        print(f"Wrote {len(results)} results to {batch_file}")

        file_number += 1


if __name__ == "__main__":
    main()
