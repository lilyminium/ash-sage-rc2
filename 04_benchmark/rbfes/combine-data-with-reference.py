"""
This script combines experimental and computed free energy data from RBFEs.
It writes two files to the output directory.
\b
The first is dG.csv with the columns:
    - System (str): name of the system (e.g. "tyk2")
    - Ligand (str): name of the ligand (e.g. "ejm_55")
    - SMILES (str): SMILES representation of the ligand
    - Mapped SMILES (str): mapped SMILES representation of the ligand
    - Value (kcal / mol) (float): the computed dG in kcal/mol
    - Uncertainty (kcal / mol) (float): the computed dG uncertainty in kcal/mol
    - Reference Value (kcal / mol) (float): the experimental dG in kcal/mol
    - Reference Uncertainty (kcal / mol) (float): the experimental dG uncertainty in kcal/mol
    - Force Field (str): the force field used for the calculation
\b
The second is ddG.csv with the columns:
    - System (str): name of the system (e.g. "tyk2")
    - Transformation (str): name of the transformation (e.g. "ejm_55 -> ejm_54")
    - Ligand 1 (str): The name of ligand_i
    - SMILES 1 (str): SMILES representation of ligand_i
    - Mapped SMILES 1 (str): mapped SMILES representation of ligand_i
    - Ligand 2 (str): The name of ligand_j
    - SMILES 2 (str): SMILES representation of ligand_j
    - Mapped SMILES 2 (str): mapped SMILES representation of ligand_j
    - Value (kcal / mol) (float): the computed ddG in kcal/mol
    - Uncertainty (kcal / mol) (float): the computed ddG uncertainty in kcal/mol
    - Reference Value (kcal / mol) (float): the experimental ddG in kcal/mol
    - Reference Uncertainty (kcal / mol) (float): the experimental ddG uncertainty in kcal/mol
    - Force Field (str): the force field used for the calculation
"""

from collections import defaultdict
import pathlib
import sys

import click
import tqdm
import networkx as nx
from rdkit import Chem
import pandas as pd
import cinnabar
from openff.units import unit
from openff.toolkit import Molecule

from loguru import logger

logger.remove()
logger.add(sys.stdout)

FORCEFIELD_NAMES = {
    "gaff-1.81": "GAFF 1.81",
    "gaff-2.11": "GAFF 2.11",
    "openff-1.2.1": "Sage 1.2.1",
    "openff-2.0.0": "Sage 2.0.0",
    "openff-2.2.1_am1bcc": "Sage 2.2.1+AM1-BCC",
    "openff-2.2.1": "Sage 2.2.1 + ELF10",
    "openff-2.2.1_nagl": "Sage 2.2.1 + AshGC",
    "openff-2.3.0rc0_elf10": "Sage 2.3.0rc0 + ELF10",
    "openff-2.3.0rc0_nagl": "Sage 2.3.0rc0 + AshGC",
    "openff-2.3.0rc1_elf10": "Sage 2.3.0rc1 + ELF10",
    "openff-2.3.0rc1_nagl": "Sage 2.3.0rc1 + AshGC",
}

def sanitize_smiles(smiles: str) -> str:
    """Sanitize SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    return ""

def generate_femap(
    experimental_data: dict[str, tuple[float, float]],
    computed_data: dict[tuple[str, str], tuple[float, float]]
):
    """
    Create a cinnabar FEMap
    
    Parameters
    ----------
    experimental_data: dict[str, tuple[float, float]]
      The experimental dG data, keyed by ligand name.
      The first element is the value in kcal/mol,
      the second element is the uncertainty in kcal/mol.
    computed_data: dict[tuple[str, str], tuple[float, float]]
      The calculated data, keyed by (ligand_i, ligand_j) pairs.
      The first element of the value is the value in kcal/mol,
      the second element is the uncertainty in kcal/mol.

    Returns
    -------
    femap : cinnabar.FEMap
      A cinnabar FEMap object with generated absolute values.
    """
    from cinnabar import Measurement, ReferenceState, FEMap

    # define a ground data
    ground = ReferenceState()

    fe_results = {'Experimental': {}, 'Calculated': []}

    # Lead the Measurements
    for ligand_name, experimental_values in experimental_data.items():
        # add experimental data
        m = Measurement(
            labelA=ground,
            labelB=ligand_name,
            DG=experimental_values[0] * unit.kilocalorie_per_mole,
            uncertainty=experimental_values[1] * unit.kilocalorie_per_mole,
            computational=False
        )
        fe_results['Experimental'][m.labelB] = m

    for (ligand_i, ligand_j), (ddG, ddG_err) in computed_data.items():
        # add computed data
        if not ligand_i in fe_results['Experimental']:
            logger.info(f"ligand not found: {ligand_i}")
            continue
        if not ligand_j in fe_results['Experimental']:
            logger.info(f"ligand not found: {ligand_j}")
            continue

        m = Measurement(
            labelA=ligand_i,
            labelB=ligand_j,
            DG=ddG * unit.kilocalorie_per_mole,
            uncertainty=ddG_err * unit.kilocalorie_per_mole,
            computational=True
        )
        fe_results['Calculated'].append(m)
    
    # Feed into the FEMap object
    femap = FEMap()

    for entry in fe_results['Experimental'].values():
        femap.add_measurement(entry)

    for entry in fe_results['Calculated']:
        femap.add_measurement(entry)

    # generate dGs
    femap.generate_absolute_values()
    return femap


@click.command(help=__doc__)
@click.option(
    "--input-structures",
    "-is",
    "input_structures_path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help=(
        "Path to the directory containing prepared structures for each system. "
        "The expected structure is that systems have the same name "
        "and ligands are saved in a file called `ligands.sdf` with names."
    ),
    default="/Users/lily/pydev/IndustryBenchmarks2024/industry_benchmarks/input_structures/prepared_structures/jacs_set"
)
@click.option(
    "--input-data",
    "-i",
    "input_data_path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data/jacs",
    help=(
        "Path to the directory containing computed data for each system. "
        "The expected structure is that subdirectories contain the system name "
        "and data files are saved in a sub-subdirectory called `data`. "
        "Computed data is saved in files with the `.ddG` extension. "
        "Experimental data is saved in files ending with `out.csv`."
    )
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="output",
    help=(
        "Path to the directory where output files will be saved. "
        "Two file are saved: dG.csv and ddG.csv"
    )
)
def main(
    input_structures_path: str = "/Users/lily/pydev/IndustryBenchmarks2024/industry_benchmarks/input_structures/prepared_structures/jacs_set",
    input_data_path: str = "data/jacs",
    output_directory: str = "output",
):

    input_structures_path = pathlib.Path(input_structures_path)
    input_data_path = pathlib.Path(input_data_path)
    reference_csvs = sorted(input_data_path.glob("*/data/*out.csv"))

    # keys in experimental_dGs_by_name are system names
    # values are dicts of {ligand_name: (dG, dG error)} in kcal/mol
    experimental_dGs_by_name: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    name_to_system: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)

    # load experimental data
    for csv in tqdm.tqdm(reference_csvs, desc="read reference"):
        df = pd.read_csv(csv, header=0, dtype={
            "Ligand name": str,
            "Exp. dG (kcal/mol)": float,
            "Exp. dG error (kcal/mol)": float,
        })
        system_name = csv.parent.parent.name

        sdf_file = input_structures_path / system_name / "ligands.sdf"
        assert sdf_file.exists(), f"SDF file {sdf_file} does not exist"

        molecules_by_name: dict[str, Molecule] = {
            mol.name: mol
            for mol in Molecule.from_file(sdf_file, "SDF")
        }
        logger.info(f"Loaded {len(molecules_by_name)} molecules from {sdf_file}")

        # placeholder value
        if not "Exp. dG error (kcal/mol)" in df.columns:
            df["Exp. dG error (kcal/mol)"] = 0

        for _, row in df.iterrows():
            ligand_name = str(row["Ligand name"])
            molecule = molecules_by_name[ligand_name]
            mapped_smiles = molecule.to_smiles(mapped=True)
            smiles = sanitize_smiles(molecule.to_smiles())

            experimental_dGs_by_name[system_name][ligand_name] = (
                row["Exp. dG (kcal/mol)"],
                row["Exp. dG error (kcal/mol)"]
            )
            name_to_system[system_name][ligand_name] = (mapped_smiles, smiles)

    ddG_rows = []
    dG_rows = []

    # load in computed data
    ddG_files = sorted(input_data_path.glob("*/data/*.ddG"))
    for ddG_file in tqdm.tqdm(ddG_files, desc="Parsing ddG files and generating dGs"):
        system_name = ddG_file.parent.parent.name
        ddg_df = pd.read_csv(
            ddG_file, header=0, sep="\t",
            dtype={
                "ligand_i": str,
                "ligand_j": str,
                "DDG(i->j) (kcal/mol)": float,
                "uncertainty (kcal/mol)": float
            }
        )
        forcefield_name = FORCEFIELD_NAMES.get(ddG_file.stem, ddG_file.stem)

        # skip invalid data
        # opls-4 MCL system has some ligands I can't find
        # and not sure what tmp is
        if forcefield_name in ("tmp", "opls-4"):
            continue

        computed_data = {}
        for _, row in ddg_df.iterrows():
            ligand_i = str(row["ligand_i"])
            ligand_j = str(row["ligand_j"])
            ddG = row["DDG(i->j) (kcal/mol)"]
            ddG_err = row["uncertainty (kcal/mol)"]
            computed_data[(ligand_i, ligand_j)] = (ddG, ddG_err)

            # create ddG row
            mapped_1, smiles_1 = name_to_system[system_name][ligand_i]
            mapped_2, smiles_2 = name_to_system[system_name][ligand_j]
            exp_value_j, exp_err_j = experimental_dGs_by_name[system_name][ligand_j]
            exp_value_i, exp_err_i = experimental_dGs_by_name[system_name][ligand_i]

            ddg_row = {
                "System": system_name,
                "Transformation": f"{smiles_1} -> {smiles_2}",
                "Ligand 1": ligand_i,
                "SMILES 1": smiles_1,
                "Mapped SMILES 1": mapped_1,
                "Ligand 2": ligand_j,
                "SMILES 2": smiles_2,
                "Mapped SMILES 2": mapped_2,
                "Value (kcal / mol)": ddG,
                "Uncertainty (kcal / mol)": ddG_err,
                "Reference Value (kcal / mol)": exp_value_j - exp_value_i,
                "Reference Uncertainty (kcal / mol)": exp_err_j + exp_err_i,
                "Force Field": forcefield_name
            }
            ddG_rows.append(ddg_row)

        experimental_data = experimental_dGs_by_name[system_name]

        # generate femap for absolute dGs
        femap = generate_femap(
            experimental_data=experimental_data,
            computed_data=computed_data
        )
        graph: nx.DiGraph = femap.to_legacy_graph()
        for ligand_name, data in graph.nodes(data=True):
            mapped_smiles, smiles = name_to_system[system_name][ligand_name]
            dG_row = {
                "System": system_name,
                "Ligand": ligand_name,
                "SMILES": smiles,
                "Mapped SMILES": mapped_smiles,
                "Value (kcal / mol)": data["calc_DG"],
                "Uncertainty (kcal / mol)": data["calc_dDG"],
                "Reference Value (kcal / mol)": data["exp_DG"],
                "Reference Uncertainty (kcal / mol)": data["exp_dDG"],
                "Force Field": forcefield_name
            }
            dG_rows.append(dG_row)

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    dG_df = pd.DataFrame(dG_rows)

    unique_systems = dG_df["System"].unique()
    logger.info(
        f"Found {len(unique_systems)} unique systems: "
        f"{', '.join(sorted(unique_systems))}"
    )

    unique_force_fields = dG_df["Force Field"].unique()
    logger.info(
        f"Found {len(unique_force_fields)} unique force fields: "
        f"{', '.join(sorted(unique_force_fields))}"
    )

    dG_file = output_directory / "dG.csv"
    dG_df.to_csv(dG_file, index=False)
    logger.info(f"Wrote {len(dG_df)} dG values to {dG_file}")

    ddG_df = pd.DataFrame(ddG_rows)
    ddG_file = output_directory / "ddG.csv"
    ddG_df.to_csv(ddG_file, index=False)
    logger.info(f"Wrote {len(ddG_df)} ddG values to {ddG_file}")


if __name__ == "__main__":
    main()
