import functools
import pickle
import multiprocessing
import pathlib
import typing

import click
import tqdm

import pyarrow.dataset as ds
import pyarrow.compute as pc
from openff.units import unit
from openff.toolkit import Molecule, ForceField
import MDAnalysis as mda

from openff.qcsubmit.results import (
    BasicResultCollection,
    OptimizationResultCollection,
    TorsionDriveResultCollection,
    OptimizationResult,
    TorsionDriveResult,
)
import qcportal as ptl
from openff.qcsubmit.utils.utils import portal_client_manager


QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

def compute_angle_difference(mm_value, qm_value):
    difference = mm_value - qm_value
    if difference > 180:
        difference -= 360
    if difference < -180:
        difference += 360
    return difference


def has_seesaw_sx4(
    item
):
    record, molecule = item
    matches = molecule.chemical_environment_matches(
        "[*:1]~[#16X4:2]~[*:3]"
    )
    if not matches:
        return True, record.id
    
    matches = [list(x) for x in matches]
    
    u = mda.Universe(molecule.to_rdkit())

    threshold = 10
    target = 180
    
    for conf in molecule.conformers:
        u.atoms.positions = conf.m_as(unit.angstrom)
        for match in matches:
            angle = u.atoms[match].angle.value()
            abs_diff = abs(
                compute_angle_difference(angle, target)
            )
            if abs_diff < threshold:
                return False, record.id
    return True, record.id


@click.command()
@click.option(
    "--input-dataset",
    "-i",
    "input_datasets",
    type=click.Path(exists=True, dir_okay=True),
    multiple=True,
    help="Input dataset(s) to filter.",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(exists=False, dir_okay=False),
    default="failed_smiles.dat",
    help="Output file to save the failed SMILES.",
)
@click.option(
    "--n-processes",
    "-np",
    type=int,
    default=1,
    help="Number of processes to use for parallel processing.",
)
@click.option(
    "--data-type",
    "-dt",
    type=click.Choice(["optimization", "torsiondrive"]),
    default="optimization",
    help="Type of data to process.",
)
def main(
    input_datasets: list[str],
    output_file: str,
    n_processes: int = 1,
    data_type: typing.Literal["optimization", "torsiondrive"] = "optimization",
):
    pickle_file = pathlib.Path(f"s-records-{data_type}.pkl")
    if pickle_file.is_file():
        with open(pickle_file, "rb") as f:
            records_and_molecules = pickle.load(f)
    else:
        all_cmiles = set()
        for input_dataset in input_datasets:
            dataset = ds.dataset(input_dataset)
            cmiles = dataset.to_table(columns=["cmiles"]).to_pydict()["cmiles"]
            all_cmiles |= set(cmiles)

        print(f"Loaded {len(all_cmiles)} unique cmiles")

        # filter for including S
        all_cmiles = [
            cmiles for cmiles in all_cmiles
            if "S" in cmiles
        ]
        print(f"Filtered to {len(all_cmiles)} cmiles with S")

        if data_type == "optimization":
            collection_klass = OptimizationResultCollection
            result_klass = OptimizationResult
        elif data_type == "torsiondrive":
            collection_klass = TorsionDriveResultCollection
            result_klass = TorsionDriveResult


        all_entries = []
        seen_ids = set()
        for input_dataset in input_datasets:
            dataset = ds.dataset(input_dataset)
            subset = dataset.filter(pc.field("cmiles").isin(all_cmiles))
            rows = subset.to_table().to_pylist()
            for row in tqdm.tqdm(rows, desc="Processing rows"):
                if row["id"] in seen_ids:
                    continue
                seen_ids.add(row["id"])
                result = result_klass(
                    cmiles=row["cmiles"],
                    inchi_key=row["inchi_key"],
                    record_id=row["id"]
                )
                all_entries.append(result)
        print(f"Loaded {len(all_entries)} entries")
        collection = collection_klass(
            entries={QCFRACTAL_URL: all_entries}
        )

        with portal_client_manager(lambda x: ptl.PortalClient(x, cache_dir=".")):
            records_and_molecules = collection.to_records()
        with open(pickle_file, "wb") as f:
            pickle.dump(records_and_molecules, f)
    
    failed = set()

    with multiprocessing.Pool(n_processes) as pool:
        results = list(
            tqdm.tqdm(
                pool.imap(has_seesaw_sx4, records_and_molecules),
                total=len(records_and_molecules),
            )
        )
        batch_failed = [
            qcarchive_id
            for success, qcarchive_id in results
            if not success
        ]
        failed |= set(batch_failed)
        with open(output_file, "a") as f:
            f.write("\n".join(list(map(str, batch_failed))))


    print(f"Total number of failed IDs: {len(failed)}")



if __name__ == "__main__":
    main()
