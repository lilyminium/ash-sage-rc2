import functools
import pathlib
import click
import tqdm
import multiprocessing

import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
import pandas as pd

from rdkit import Chem
import qcportal as ptl
from openff.toolkit import Molecule

from openff.qcsubmit.results import (
    BasicResultCollection,
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.qcsubmit.utils.utils import portal_client_manager
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    RecordStatusFilter,
    UnperceivableStereoFilter,
    HydrogenBondFilter,
    ElementFilter,
)


QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

IGNORE_IODINE = [
    "OpenFF Discrepancy Benchmark 1",
    "OpenFF Gen 2 Opt Set 2 Coverage",
    "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy",
    "SMIRNOFF Coverage Set 1",
    "OpenFF Ehrman Informative Optimization v0.2",
    "FDA optimization dataset 1",
    "Kinase Inhibitors: WBO Distributions",

    "OpenFF Gen 2 Torsion Set 2 Coverage 2",
    "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2",
]

OPTIMIZATION_WHITELISTS = [
    "OpenFF Optimization Set 1",
    "SMIRNOFF Coverage Set 1",
    "OpenFF VEHICLe Set 1",
    "OpenFF Discrepancy Benchmark 1",
    "OpenFF Ehrman Informative Optimization v0.2",
    "Pfizer discrepancy optimization dataset 1",
    "FDA optimization dataset 1",
    "Kinase Inhibitors: WBO Distributions",
    "OpenFF Gen 2 Opt Set 1 Roche",
    "OpenFF Gen 2 Opt Set 2 Coverage",
    "OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy",
    "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy",
    "OpenFF Gen 2 Opt Set 5 Bayer",
    "OpenFF Sandbox CHO PhAlkEthOH v1.0",
    "OpenFF Aniline Para Opt v1.0",
    "OpenFF Industry Benchmark Season 1 v1.2",
    "OpenFF Gen2 Optimization Dataset Protomers v1.0",
    "OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0",
    "OpenFF Iodine Chemistry Optimization Dataset v1.0",
    "XtalPi Shared Fragments OptimizationDataset v1.0",
    "XtalPi 20-percent Fragments OptimizationDataset v1.0",
    "OpenFF Torsion Benchmark Supplement v1.0",
    "OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0",
    "OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0",
    "OpenFF Iodine Fragment Opt v1.0",
    "OpenFF Sulfur Optimization Training Coverage Supplement v1.0",
    "OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0",
    "OpenFF Lipid Optimization Training Supplement v1.0",
    "OpenFF Lipid Optimization Benchmark Supplement v1.0",
    "OpenFF Cresset Additional Coverage Optimizations v4.0",
    "OpenFF Protein PDB 4-mers v4.0",
    "OpenFF Additional Generated ChEMBL Optimizations v4.0"
]

TORSIONDRIVE_WHITELISTS = [
    "OpenFF Group1 Torsions",
    "OpenFF Group1 Torsions 2",
    "OpenFF Group1 Torsions 3",
    "SMIRNOFF Coverage Torsion Set 1",
    "OpenFF Substituted Phenyl Set 1",
    "Pfizer discrepancy torsion dataset 1",
    "OpenFF Primary Benchmark 1 Torsion Set",
    "OpenFF Gen 2 Torsion Set 1 Roche 2",
    "OpenFF Gen 2 Torsion Set 2 Coverage 2",
    "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2",
    "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2",
    "OpenFF Gen 2 Torsion Set 5 Bayer 2",
    "OpenFF Gen 2 Torsion Set 6 supplemental 2",
    "OpenFF Fragmenter Validation 1.0",
    "OpenFF DANCE 1 eMolecules t142 v1.0",
    "OpenFF Rowley Biaryl v1.0",
    "OpenFF-benchmark-ligand-fragments-v1.0",
    "OpenFF Protein Fragments TorsionDrives v1.0",
    "OpenFF WBO Conjugated Series v1.0",
    "OpenFF Amide Torsion Set v1.0",
    "OpenFF-benchmark-ligand-fragments-v2.0",
    "OpenFF multiplicity correction torsion drive data v1.1",
    "OpenFF Protein Capped 3-mer Omega v1.0",
    "XtalPi Shared Fragments TorsiondriveDataset v1.0",
    "OpenFF Torsion Coverage Supplement v1.0",
    "OpenFF RNA Dinucleoside Monophosphate TorsionDrives v1.0",
    "XtalPi 20-percent Fragments TorsiondriveDataset v1.0",
    "OpenFF Torsion Drive Supplement v1.0",
    "OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0",
    "OpenFF Phosphate Torsion Drives v1.0",
    "OpenFF Alkane Torsion Drives v1.0",
    "OpenFF Cresset Additional Coverage TorsionDrives v4.0",
    "OpenFF Additional Generated ChEMBL TorsionDrives 4.0",
    "OpenFF Additional Generated Guanidine Amidine Derivative and O-Linker TorsionDrives 4.0",
    "OpenFF Gen3 Torsion Set v1.0"
]

@functools.cache
def sanitize_smiles(smiles: str) -> str:
    """Sanitize a SMILES string by removing atom map numbers and Hs."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def get_canonical_smiles(entry) -> str:
    """Get the canonical SMILES from the entry."""
    KEY = "canonical_isomeric_explicit_hydrogen_mapped_smiles"
    try:
        return entry.attributes[KEY]
    except KeyError:
        return getattr(entry.initial_molecule.identifiers, KEY)

@functools.cache
def cmiles_to_inchi(cmiles: str) -> str:
    return Molecule.from_mapped_smiles(
        cmiles,
        allow_undefined_stereo=True
    ).to_inchikey(fixed_hydrogens=True)


def single_process_optimization(item):
    record, molecule, cmiles, dataset_name = item
    try:
        if not record.status.upper() == "COMPLETE":
            return
        if not ConnectivityFilter()._filter_function(None, None, molecule):
            return
        if not UnperceivableStereoFilter()._filter_function(None, None, molecule):
            return
    # can have some rdkit connection errors
    except:
        return

    try:
        smiles = sanitize_smiles(cmiles)
    except ValueError as e:
        print(e)
        return []
    try:
        inchi_key = cmiles_to_inchi(cmiles)
    except:
        return []
    energy = record.energies[-1]
    data_entry = {
        "id": record.id,
        "inchi_key": inchi_key,
        "cmiles": cmiles,
        "smiles": smiles,
        "dataset_name": dataset_name,
        "energy": energy
    }
    return data_entry



def download_optimization(client, dataset_name: str, n_processes: int = 4) -> pa.Table:
    """
    Download an Optimization dataset from QCArchive and return a pyarrow table.
    """
    data_entries = []
    try:

        optimization_result_collection = OptimizationResultCollection.from_server(
            client=client,
            datasets=[dataset_name],
            spec_name="default"
        )
    except KeyError:
        optimization_result_collection = OptimizationResultCollection.from_server(
            client=client,
            datasets=[dataset_name],
            spec_name="spec_1"
        )

    ids_to_cmiles = {
        entry.record_id: entry.cmiles
        for entries in optimization_result_collection.entries.values()
        for entry in entries
    }
    with portal_client_manager(lambda x: ptl.PortalClient(x, cache_dir=".")):
        records_and_molecules = optimization_result_collection.to_records()

    items = [
        (record, molecule, ids_to_cmiles[record.id], dataset_name)
        for record, molecule in records_and_molecules
    ]
    with multiprocessing.Pool(n_processes) as pool:
        results = list(
            tqdm.tqdm(
                pool.imap(single_process_optimization, items),
                total=len(items),
                desc="Processing records and molecules"
            )
        )
    data_entries = [result for result in results if result is not None]

    table = pa.Table.from_pylist(data_entries)
    return table

def single_process_torsiondrive(item):
    record, molecule, cmiles, dataset_name = item
    try:
        if not record.status.upper() == "COMPLETE":
            return
        if not ConnectivityFilter()._filter_function(None, None, molecule):
            return
        if not UnperceivableStereoFilter()._filter_function(None, None, molecule):
            return
        if not HydrogenBondFilter()._filter_function(None, None, molecule):
            return
    except:
        return

    try:
        smiles = sanitize_smiles(cmiles)
    except ValueError as e:
        print(e)
        return []
    try:
        inchi_key = cmiles_to_inchi(cmiles)
    except:
        return []
    dihedrals = record.specification.keywords.dihedrals
    dihedral = []
    for dih in dihedrals:
        dihedral.extend(dih)
    data_entry = {
        "id": record.id,
        "inchi_key": inchi_key,
        "cmiles": cmiles,
        "smiles": smiles,
        "dataset_name": dataset_name,
        "dihedral": list(dihedral),
        "n_dihedrals": len(dihedrals),
    }
    return data_entry


def download_torsiondrive(client, dataset_name: str, n_processes: int = 4) -> pa.Table:
    """
    Download a TorsionDrive dataset from QCArchive and return a pyarrow table.
    """
    data_entries = []

    try:
        torsion_result_collection = TorsionDriveResultCollection.from_server(
            client=client,
            datasets=[dataset_name],
            spec_name="default"
        )
    except KeyError:
        torsion_result_collection = TorsionDriveResultCollection.from_server(
            client=client,
            datasets=[dataset_name],
            spec_name="spec_1"
        )

    ids_to_cmiles = {
        entry.record_id: entry.cmiles
        for entries in torsion_result_collection.entries.values()
        for entry in entries
    }
    with portal_client_manager(lambda x: ptl.PortalClient(x, cache_dir=".")):
        records_and_molecules = torsion_result_collection.to_records()

    items = [
        (record, molecule, ids_to_cmiles[record.id], dataset_name)
        for record, molecule in records_and_molecules
    ]
    with multiprocessing.Pool(n_processes) as pool:
        results = list(
            tqdm.tqdm(
                pool.imap(single_process_torsiondrive, items),
                total=len(items),
                desc="Processing records and molecules"
            )
        )
    data_entries = [result for result in results if result is not None]

    table = pa.Table.from_pylist(data_entries)
    return table


@click.command()
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="data/tables",
    help="Directory to save the downloaded tables."
)
def main(
    output_directory: str = "data/tables"
):
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    client = ptl.PortalClient(address=QCFRACTAL_URL, cache_dir=".")
    
    # download optimizations
    opt_directory = output_directory / "optimization"
    opt_directory.mkdir(parents=True, exist_ok=True)
    for dsname in tqdm.tqdm(OPTIMIZATION_WHITELISTS, desc="Downloading Optimizations"):
        table = download_optimization(client, dsname)
        # filter out iodines
        if dsname in IGNORE_IODINE:
            df = table.to_pandas()
            mask = np.array(["I" in smi for smi in df.smiles.values])
            df = pd.DataFrame(df[~mask])
            table = pa.Table.from_pandas(df)

        table_file = opt_directory / f"{dsname}.parquet"
        pq.write_table(table, table_file)
        print(f"Saved {table_file}")

    # # download torsiondrives
    # td_directory = output_directory / "torsiondrive"
    # td_directory.mkdir(parents=True, exist_ok=True)
    # for dsname in tqdm.tqdm(TORSIONDRIVE_WHITELISTS, desc="Downloading TorsionDrives"):
    #     table = download_torsiondrive(client, dsname)
    #     if dsname in IGNORE_IODINE:
    #         df = table.to_pandas()
    #         mask = np.array(["I" in smi for smi in df.smiles.values])
    #         df = pd.DataFrame(df[~mask])
    #         table = pa.Table.from_pandas(df)

    #     table_file = td_directory / f"{dsname}.parquet"
    #     pq.write_table(table, table_file)
    #     print(f"Saved {table_file}")

    
        
if __name__ == "__main__":
    main()
