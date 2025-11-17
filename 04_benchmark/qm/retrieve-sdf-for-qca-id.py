import click

from rdkit import Chem
from openff.toolkit import Molecule, unit

import numpy as np
import pyarrow.dataset as ds
import pyarrow.compute as pc

# hardcode the methods we're interested in
METHODS = {
    "fb-fit-v3-single-mean-k100_unconstrained": "Sage 2.3.0",
    "openff_unconstrained-2.1.0": "Sage 2.1.0",
    "openff_unconstrained-2.2.1": "Sage 2.2.1",
    "qm": "QM",
}

@click.command()
@click.option(
    "--qca-id",
    "-q",
    type=int,
    required=True,
    help="The QCArchive ID to write SDF for.",
)
@click.option(
    "--output-path",
    "-o",
    type=click.Path(dir_okay=False, file_okay=True, exists=False),
    required=True,
    help="The output SDF path.",
)
@click.option(
    "--input-path",
    "-i",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data/optimization",
    help="The input data directory.",
)
def main(
    qca_id: int,
    output_path: str,
    input_path: str = "data/optimization",
):
    dataset = ds.dataset(input_path)
    expression = (
        (pc.field("qcarchive_id") == qca_id)
        & (pc.field("method").isin(list(METHODS.keys())))
    )

    rows = dataset.filter(expression).to_table(
        columns=[
            "mapped_smiles",
            "cmiles",
            "coordinates",
            "energy",
            "dataset",
            "method"
        ]
    ).to_pylist()
    rows = sorted(rows, key=lambda r: METHODS[r["method"]])

    writer = Chem.SDWriter(output_path)
    for row in rows:
        mol = Molecule.from_mapped_smiles(row["mapped_smiles"], allow_undefined_stereo=True)
        mol._conformers = [np.array(row["coordinates"]).reshape((-1, 3)) * unit.angstrom]
        rdmol = mol.to_rdkit()
        rdmol.SetProp("cmiles", row["cmiles"])
        rdmol.SetDoubleProp("energy", row["energy"])
        rdmol.SetProp("dataset", row["dataset"])
        rdmol.SetProp("method", row["method"])
        rdmol.SetProp("FF", METHODS[row["method"]])
        writer.write(rdmol)

    writer.close()


if __name__ == "__main__":
    main()
