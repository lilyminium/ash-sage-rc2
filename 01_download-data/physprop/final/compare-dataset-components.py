from collections import defaultdict, Counter
import pathlib
import json
import logging
import typing

import click
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit.Chem import Draw
from openff.toolkit import Molecule
from openff.evaluator.datasets import PhysicalPropertyDataSet
from yammbs.checkmol import analyze_functional_groups

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")


def save_highres_image(smiles: list[str], filename: str):
    rdmols = [Chem.MolFromSmiles(smi) for smi in smiles if smi]
    img = Draw.MolsToGridImage(
        rdmols,
        molsPerRow=5,
        subImgSize=(300, 300),
        useSVG=True,
    )
    with open(filename, "w") as f:
        f.write(img)#.data)
    print(f"Saved high-resolution image to {filename}")

@click.command()
@click.option(
    "--new-dataset-file",
    "-n",
    default="new-dataset.json",
    help="The JSON file containing the new dataset",
)
@click.option(
    "--old-dataset-file",
    "-o",
    default="sage-training-set.json",
    help="The JSON file containing the old dataset",
)
@click.option(
    "--output-file",
    "-f",
    default="functional-group-profiles/unique-components.json",
    help="The output file to save the unique components to",
)
@click.option(
    "--image-directory",
    "-i",
    default="images",
    help="Directory to save the high-resolution images of added/removed components",
)
def main(
    new_dataset_file: str = "new-dataset.json",
    old_dataset_file: str = "sage-training-set.json",
    output_file: str = "functional-group-profiles/unique-components.json",
    image_directory: str = "images",
):
    new_dataset = PhysicalPropertyDataSet.from_json(new_dataset_file)
    old_dataset = PhysicalPropertyDataSet.from_json(old_dataset_file)

    new_unique_components = set()
    old_unique_components = set()

    for physprop in old_dataset.properties:
        for component in physprop.substance.components:
            old_unique_components.add(component.smiles)

    for physprop in new_dataset.properties:
        for component in physprop.substance.components:
            new_unique_components.add(component.smiles)

    print(f"Old unique components: {len(old_unique_components)}")
    print(f"New unique components: {len(new_unique_components)}")


    added_components = new_unique_components - old_unique_components
    removed_components = old_unique_components - new_unique_components

    print(f"Added components: {len(added_components)}")
    print(f"Removed components: {len(removed_components)}")


    data = {
        "old_unique_components": sorted(old_unique_components),
        "new_unique_components": sorted(new_unique_components),
        "added_components": sorted(added_components),
        "removed_components": sorted(removed_components),
    }
    with open(output_file, "w") as f:
        json.dump(data, f, indent=4)

    print(f"Saved unique components to {output_file}")

    # plot high-resolution RDKit image
    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)
    save_highres_image(added_components, str(image_directory / "highres_added.svg"))
    save_highres_image(removed_components, str(image_directory / "highres_removed.svg"))


if __name__ == "__main__":
    main()
