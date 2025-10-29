"""
This script compares the training set of SMILES strings against the
maximum common substructures (MCS) of ligands in a given RBFE directory. It identifies
which training set molecules match the MCS of ligands in each system and
outputs the results to JSON files. Optionally, it can also generate images
of the MCS cores.
\b
It saves the following outputs:
- ligand_cores.json: A JSON file mapping each group and system to its MCS SMARTS string.
- training_set_overlaps.json: A JSON file mapping each training set SMILES to the
  groups and systems whose MCS it matches.
- unmatched_systems.json: A JSON file listing systems whose MCS did not match
  any training set SMILES.
"""

from collections import defaultdict
import pathlib
import json

import click
import tqdm
from loguru import logger
from openff.toolkit import Molecule

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS

def get_all_mcs(
    input_directory: pathlib.Path,
    image_directory: pathlib.Path | None = None,
    plot_cores: bool = True,
) -> dict[str, dict[str, str]]:
    """Get the maximum common substructure for all ligands in each system.
    
    Parameters
    ----------
    input_directory : pathlib.Path
        The input directory containing ligand files.
    image_directory : pathlib.Path | None, optional
        The directory to save images of the cores, by default None.
    plot_cores : bool, optional
        Whether to plot the cores, by default True.
        If False or image_directory is None, no images will be saved.
    
    Returns
    -------
    dict[str, dict[str, str]]
        A nested dictionary mapping group names, to system names, to SMARTS strings of the MCS.
    """
    ligand_files = sorted(input_directory.glob("*/*/ligands.sdf"))
    if not ligand_files:
        raise ValueError(
            f"No ligand files found in {input_directory}. "
            "Expected structure: <input_directory>/<group_name>/<system_name>/ligands.sdf"
        )

    logger.info(f"Found {len(ligand_files)} ligand files.")

    ligand_patterns: dict[str, dict[str, str]] = defaultdict(dict)
    for ligand_file in tqdm.tqdm(ligand_files, desc="Finding ligand cores"):
        system = ligand_file.parent.stem
        group = ligand_file.parent.parent.stem

        molecules = Molecule.from_file(ligand_file, "SDF", allow_undefined_stereo=True)
        rdmols = [mol.to_rdkit() for mol in molecules]

        res = rdFMCS.FindMCS(rdmols, maximizeBonds=True, matchValences=True)
        ligand_patterns[group][system] = res.smartsString
        
        if plot_cores and image_directory is not None:
            draw_rdmols = [Chem.MolFromSmarts(res.smartsString)]
            legends = ["core"]
            for i, rdmol in enumerate(rdmols):
                draw_rdmols.append(Chem.MolFromSmiles(Chem.MolToSmiles(rdmol)))
                legends.append(f"{group}-{system}-{i}")
            img = Draw.MolsToGridImage(draw_rdmols, subImgSize=(500,500), legends=legends)
            groupdir = image_directory / group
            groupdir.mkdir(exist_ok=True, parents=True)
            systemfile = groupdir / f"{system}.png"
            with systemfile.open("wb+") as f:
                f.write(img.data)
            logger.info(f"Wrote core image to {systemfile}")

    return ligand_patterns


@click.command()
@click.option(
    "--rbfe-directory",
    "-r",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=pathlib.Path),
    help=(
        "The RBFE directory containing ligand files. "
        "Expected structure: <rbfe_directory>/<group_name>/<system_name>/ligands.sdf"
    )
)
@click.option(
    "--training-smiles-file",
    "-t",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
    help="A file containing SMILES strings of the training set, one per line.",
)
@click.option(
    "--output-directory",
    "-o",
    default=".",
    type=click.Path(file_okay=False, path_type=pathlib.Path),
    help="The output directory to write results to.",
)
@click.option(
    "--image-directory",
    "-im",
    default=None,
    type=click.Path(file_okay=False, path_type=pathlib.Path),
    help="The directory to save images of the cores. If not provided, no images will be saved.",
)
@click.option(
    "--no-plot-cores/--plot-cores",
    "plot_cores",
    default=True,
    is_flag=True,
    help="Whether to plot the cores. If false, no images will be saved.",
)
def main(
    rbfe_directory: pathlib.Path,
    training_smiles_file: pathlib.Path,
    output_directory: pathlib.Path = ".",
    image_directory: pathlib.Path | None = None,
    plot_cores: bool = True,
):
    rbfe_directory = pathlib.Path(rbfe_directory)
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)
    if image_directory is not None:
        image_directory = pathlib.Path(image_directory)
        image_directory.mkdir(exist_ok=True, parents=True)

    logger.info(f"Analyzing RBFE directory: {rbfe_directory}")
    ligand_cores = get_all_mcs(
        input_directory=rbfe_directory,
        image_directory=image_directory,
        plot_cores=plot_cores,
    )

    core_patterns_file = output_directory / "ligand_cores.json"
    with core_patterns_file.open("w") as f:
        json.dump(ligand_cores, f, indent=4)
    logger.info(f"Wrote ligand cores to {core_patterns_file}")

    pattern_to_grouping: dict[str, tuple[str, str]] = {}
    for group, subdict in ligand_cores.items():
        for system, pattern in subdict.items():
            pattern_to_grouping[pattern] = (group, system)

    matches = defaultdict(list)

    with open(training_smiles_file, "r") as f:
        training_smiles = [line.strip() for line in f.readlines() if line.strip()]
    for smi in tqdm.tqdm(training_smiles, desc="Checking training set for matches"):
        mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        for pattern, (group, system) in pattern_to_grouping.items():
            if mol.chemical_environment_matches(pattern):
                matches[smi].append((group, system))

    output_file = output_directory / "training_set_overlaps.json"
    with output_file.open("w") as f:
        json.dump(matches, f, indent=4)
    logger.info(f"Wrote training set overlaps to {output_file}")

    # get cores that had no matches
    cores_to_matches = defaultdict(lambda: defaultdict(list))
    for smi, smi_matches in matches.items():
        for group, system in smi_matches:
            cores_to_matches[group][system].append(smi)

    unmatched_systems = {}

    for group, subgroup in ligand_cores.items():
        matchgroup = cores_to_matches[group]
        missing = set(subgroup) - set(matchgroup)
        unmatched_systems[group] = sorted(list(missing))
    
    unmatched_file = output_directory / "unmatched_systems.json"
    with unmatched_file.open("w") as f:
        json.dump(unmatched_systems, f, indent=4)
    logger.info(f"Wrote unmatched systems to {unmatched_file}")

    for group, systems in unmatched_systems.items():
        logger.info(f"Group '{group}' has {len(systems)} unmatched systems: {systems}")

if __name__ == "__main__":
    main()