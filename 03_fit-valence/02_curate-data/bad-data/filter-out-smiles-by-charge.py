import functools
import multiprocessing
import pathlib

import click
import tqdm

import pyarrow.dataset as ds
from openff.toolkit import Molecule, ForceField


@functools.cache
def can_charge_smiles(smiles: str) -> bool:
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    try:
        mol.assign_partial_charges("am1bccelf10")
    except:
        return False
    return True

@functools.cache
def can_charge_cmiles(cmiles: str) -> bool:
    mol = Molecule.from_mapped_smiles(cmiles, allow_undefined_stereo=True)
    return (can_charge_smiles(mol.to_smiles()), cmiles)


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
def main(
    input_datasets: list[str],
    output_file: str,
    n_processes: int = 1,
):
    all_cmiles = set()
    for input_dataset in input_datasets:
        dataset = ds.dataset(input_dataset)
        cmiles = dataset.to_table(columns=["cmiles"]).to_pydict()["cmiles"]
        all_cmiles |= set(cmiles)

    print(f"Loaded {len(all_cmiles)} unique cmiles")

    failed = set()

    output_file = pathlib.Path(output_file)
    if output_file.exists():
        with open(output_file, "r") as f:
            failed = set([
                x.strip() for x in f.readlines()
            ])
        print(f"Loaded {len(failed)} failed cmiles from {output_file}")
        all_cmiles -= failed
        print(f"Remaining {len(all_cmiles)} cmiles to process")

    all_cmiles = sorted(all_cmiles)

    batch_size = 100
    n_batches = len(all_cmiles) // batch_size + 1
    print(f"Splitting into {n_batches} batches of size {batch_size}")
    all_cmiles = [
        all_cmiles[i * batch_size:(i + 1) * batch_size]
        for i in range(n_batches)
    ]
    print(f"Total number of batches: {len(all_cmiles)}")

    with multiprocessing.Pool(n_processes) as pool:
        for batch in tqdm.tqdm(all_cmiles, desc="Processing batches"):
            results = list(
                tqdm.tqdm(
                    pool.imap(can_charge_cmiles, batch),
                    total=len(batch),
                )
            )
            batch_failed = [
                cmiles
                for success, cmiles in results
                if not success
            ]
            failed |= set(batch_failed)
            with open(output_file, "a") as f:
                f.write("\n".join(batch_failed))

    # with multiprocessing.Pool(n_processes) as pool:
    #     results = list(
    #         tqdm.tqdm(
    #             pool.imap(can_charge_cmiles, all_cmiles),
    #             total=len(all_cmiles),
    #         )
    #     )
    print(f"Total number of failed cmiles: {len(failed)}")
    with open(output_file, "w") as f:
        f.write("\n".join(failed))


if __name__ == "__main__":
    main()
