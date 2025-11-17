import multiprocessing
import pathlib
import tqdm
import time
import sys
from loguru import logger

from click_option_group import optgroup

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

import MDAnalysis as mda
from openff.toolkit import Molecule

logger.remove()
logger.add(sys.stdout)

def compute_single(
    row,
):
    from openff.units import unit
    mol = Molecule.from_mapped_smiles(
        row["mapped_smiles"],
        allow_undefined_stereo=True
    )
    positions = np.array(row["coordinates"]).reshape((-1, 3))
    mol._conformers = [positions * unit.angstrom]

    base_entry = {
        "qcarchive_id": row["qcarchive_id"],
        "mapped_smiles": row["mapped_smiles"],
        "method": row["method"]
    }

    u = mda.Universe(mol.to_rdkit(), to_guess=["angles", "dihedrals", "impropers"])

    entries = []

    # compute impropers from OpenFF
    for improper in mol.impropers:
        atom_indices = [atom.molecule_atom_index for atom in improper]
        off_indices = list(atom_indices)
        middle = atom_indices.pop(1)
        mda_order = [middle, *atom_indices]
        angle_deg = u.atoms[mda_order].improper.value()

        entry = dict(base_entry)
        entry.update(
            {
                "topology": "ImproperTorsions",
                "atom_indices": off_indices,
                "value": angle_deg,
            }
        )
        entries.append(entry)
    return entries




def main(
    input_directory: str = "data/optimization",
    output_directory: str = "improper-values",
    n_processes: int = 8,
):

    logger.info(f"{time.ctime()} - Starting batch optimization")
    start_time = time.time()

    input_directory = pathlib.Path(input_directory)
    input_dataset = ds.dataset(input_directory)
    logger.info(f"Loaded {input_dataset.count_rows()} rows from {input_directory}")

    unique_methods = set(input_dataset.to_table(
            columns=["method"]
        ).to_pydict()["method"]
    )

    logger.info(f"Found {len(unique_methods)} unique methods")
    

    for method in unique_methods:
        output_directory = pathlib.Path(output_directory) / method
        output_directory.mkdir(parents=True, exist_ok=True)
        output_dataset = ds.dataset(output_directory)
        n_files = 0
        expression = pc.field("method") == method

        if output_dataset.count_rows() > 0:
            seen_ids = set(output_dataset.to_table(
                columns=["qcarchive_id"]
            ).to_pydict()["qcarchive_id"]
            )
            logger.info(f"Method {method} already has {len(seen_ids)} processed entries")
            expression = (expression) & (~pc.field("qcarchive_id").isin(seen_ids))
            n_files = len(output_dataset.files)

        all_rows = input_dataset.to_table(
            filter=expression,
            columns=["qcarchive_id", "mapped_smiles", "coordinates", "method"]
        ).to_pylist()
        logger.info(f"Processing {len(all_rows)} rows")

        batch_size = 5000
        for i in range(0, len(all_rows), batch_size):
            batch_rows = all_rows[i : i + batch_size]
            logger.info(f"Processing batch {i // batch_size} with {len(batch_rows)} rows")
            entries = []
            with multiprocessing.Pool(n_processes) as pool:
                for entry in tqdm.tqdm(pool.imap(compute_single, batch_rows), total=len(batch_rows)):
                    entries.extend(entry)
            if entries:
                table = pa.Table.from_pylist(entries)
                output_path = output_directory / f"batch-{n_files:03d}.parquet"
                pq.write_table(table, output_path)
                logger.info(f"Wrote {len(entries)} entries to {output_path}")
                n_files += 1

    
    logger.info(f"{time.ctime()} - Finished batch optimization")
    elapsed_time = time.time() - start_time
    logger.info(f"Elapsed time: {elapsed_time / 60:.2f} min")
    logger.info("Done!")


if __name__ == "__main__":
    main()
