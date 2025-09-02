"""
Generate an alchemical network for solvation free energy calculations for a set of ligands.
"""


import pathlib
import sys

import zstandard
import click
from gufe import (
    AlchemicalNetwork,
    ChemicalSystem,
    SmallMoleculeComponent,
    Transformation,
)
from openff.units import unit
from rdkit import Chem
import pandas as pd
from loguru import logger

from pontibus.components import ExtendedSolventComponent
from pontibus.protocols.solvation import ASFEProtocol, ASFESettings
from pontibus.protocols.solvation.settings import PackmolSolvationSettings

logger.remove()
logger.add(sys.stdout)

def get_nonwater_settings(forcefield: str) -> ASFESettings:
    """
    Get the settings for solvation free energy calculations.
    These settings are tuned to give reasonable results on the FreeSolv dataset
    with a balance of accuracy and speed.
    """
    # The settings here are effectively the "fast" settings
    # shown in the validation.
    settings: ASFEProtocol = ASFEProtocol.default_settings()
    # Because it's Alchemiscale, you set protocol_repeats to 1 and then
    # run the Transformation task multiple times to get repeats.
    # Locally, the recommendation would be to set this to 3 so that you can
    # get a standard deviation uncertainty. It's not super necessary since
    # SFEs converge well, but hey with Alchemiscale why not?!
    settings.protocol_repeats = 1
    settings.solvent_forcefield_settings.forcefields = [
        # To use a custom force field, just pass an OFFXML string
        # just like you would to openff.toolkit.ForceField
        forcefield,
    ]
    settings.vacuum_forcefield_settings.forcefields = [
        forcefield,  # as above
    ]
    settings.vacuum_engine_settings.compute_platform = "CUDA"
    settings.solvent_engine_settings.compute_platform = "CUDA"
    settings.solvation_settings = PackmolSolvationSettings(
        # In our tests 750 gave quasi equivalent results to the 1999 used in the Sage
        # benchmarks
        number_of_solvent_molecules=750,
        box_shape="cube",
        # We set assign_solvent_charges to True because we don't have LibraryCharges.
        # If False it will only attempt to use LibraryCharges.
        # Note that if True and you don't have any charges on the SmallMoleculeComponent
        # passed to ExtendedSolventComponent, the Protocol will attempt to automatically
        # assign partial charges (default is AmberTools am1bcc, but it's controllable
        # using `partial_charge_settings`.
        assign_solvent_charges=True,
        solvent_padding=None,
    )
    # Below are the default lambda & replica settings, so you don't have to
    # actually set them, but if you want to change things, you can alter them
    # by defining them this way. Note: you have to update n_replica to match
    # the number of lambda windows (and all lambda window lists must be of the same length).
    settings.lambda_settings.lambda_elec = [
        0.0,
        0.25,
        0.5,
        0.75,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    ]
    settings.lambda_settings.lambda_vdw = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.12,
        0.24,
        0.36,
        0.48,
        0.6,
        0.7,
        0.77,
        0.85,
        1.0,
    ]
    settings.lambda_settings.lambda_restraints = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    settings.vacuum_simulation_settings.n_replicas = 14
    settings.solvent_simulation_settings.n_replicas = 14
    # This set the time per replica exchange, the default is 1 ps but
    # hurts performance. I would recommend 2.5 ps
    settings.solvent_simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    settings.vacuum_simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    # Below are the default simulation lengths we use in Pontibus,
    # so you don't need to set them. However, you can do so manually
    # This is the pre-alchemical equilibration lengths
    # NVT equilibration -> NPT equilibration -> NPT "production" (more equilibration)
    # In vacuum, we set the NVT equilibration to None since it's all gas phase
    settings.solvent_equil_simulation_settings.equilibration_length_nvt = (
        0.5 * unit.nanosecond
    )
    settings.solvent_equil_simulation_settings.equilibration_length = (
        0.5 * unit.nanosecond
    )
    settings.solvent_equil_simulation_settings.production_length = 9.5 * unit.nanosecond
    settings.vacuum_equil_simulation_settings.equilibration_length_nvt = None
    settings.vacuum_equil_simulation_settings.equilibration_length = (
        0.2 * unit.nanosecond
    )
    settings.vacuum_equil_simulation_settings.production_length = 0.5 * unit.nanosecond
    # This is the alchemical equilibration length
    settings.solvent_simulation_settings.equilibration_length = 1.0 * unit.nanosecond
    settings.vacuum_simulation_settings.equilibration_length = 0.5 * unit.nanosecond
    # This is the alchemical production length
    settings.solvent_simulation_settings.production_length = 10.0 * unit.nanosecond
    settings.vacuum_simulation_settings.production_length = 2.0 * unit.nanosecond
    return settings


def get_transformation(system, forcefield: str) -> Transformation:
    """Get single transformation for solvation free energy calculation."""
    settings = get_nonwater_settings(forcefield)

    # An SFE transformation in GUFE formalism is defined as
    # going from a solute + solvent (stateA) to just solvent (stateB)
    # The vacuum states are created automatically, and this transformation
    # includes both the vacuum and solvent legs
    stateA = system
    stateB = ChemicalSystem({"solvent": system.components["solvent"]})
    protocol = ASFEProtocol(settings=settings)
    return Transformation(
        stateA=stateA, stateB=stateB, mapping=None, protocol=protocol, name=stateA.name
    )


def load_ligands_to_dict(ligands: pathlib.Path) -> dict[str, SmallMoleculeComponent]:
    """
    Load ligands from an SDF file into a dictionary of SmallMoleculeComponents,
    keyed by their SMILES strings.
    This allows easy lookup of molecules by their SMILES.
    
    Parameters
    ----------
    ligands : pathlib.Path
        Path to the SDF file containing the ligands.
    """
    molecules = {}

    rdmols = Chem.SDMolSupplier(ligands, removeHs=False)

    for mol in rdmols:
        if "SMILES" not in mol.GetPropNames():
            continue
        smiles = mol.GetProp("SMILES")

        molecules[smiles] = SmallMoleculeComponent(mol)

    return molecules


def get_chemical_systems(
    smcs: dict[str, SmallMoleculeComponent],
    csv_benchmark_data: pathlib.Path,
) -> list[ChemicalSystem]:
    """
    Using the benchmark data file, create a set of ChemicalSystems
    that contain the solute and solvent Components.

    Parameters
    ----------
    smcs : dict[str, SmallMoleculeComponent]
        Dictionary of SmallMoleculeComponents keyed by their SMILES strings.
    csv_benchmark_data : pathlib.Path
        Path to the CSV file containing the benchmark data.

    Returns
    -------
    list[ChemicalSystem]
        List of ChemicalSystems for the benchmark.
    """

    benchmark_data = pd.read_csv(csv_benchmark_data)

    systems = []
    for _index, entry in benchmark_data.iterrows():
        # ExtendedSolventComponent is a special case of SolventComponent
        # it takes a SmallMoleculeComponent on construction and retains
        # its properties. Technically you could also just use a standard
        # SolventComponent, but this allows you to define the solvent's
        # conformer before packing, and also pass through user charges.
        assert entry["Role 1"] == "Solvent"
        solvent = ExtendedSolventComponent(solvent_molecule=smcs[entry["Component 1"]])
        solute = smcs[entry["Component 2"]]

        csystem = ChemicalSystem(
            {
                "solute": solute,
                "solvent": solvent,
            },
            name=f"{entry["Component 2"]} | {entry["Component 1"]}",
        )
        systems.append(csystem)

    return systems


@click.command(help=__doc__)
@click.option(
    "--input-ligands",
    "-i",
    type=click.Path(exists=True, path_type=pathlib.Path, dir_okay=False),
    required=True,
    help="Path to an SDF file containing the ligands to include in the network.",
)
@click.option(
    "--csv-file",
    "-c",
    type=click.Path(exists=True, path_type=pathlib.Path, dir_okay=False),
    required=True,
    default=pathlib.Path("data/sage-fsolv-test-v1.csv"),
    show_default=True,
    help=(
        "Path to a CSV file containing the solvent information. "
        "Only the solutes in this file will be included in the network."
    ),
)
@click.option(
    "--forcefield-file",
    "-ff",
    type=str,
    required=True,
    help="Path to the forcefield file to use.",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(path_type=pathlib.Path, dir_okay=False),
    default=pathlib.Path("alchemical_network.json"),
    show_default=True,
    help="Path to write the alchemical network JSON file.",
)
def main(
    input_ligands: pathlib.Path,
    csv_file: pathlib.Path,
    forcefield_file: str,
    output_file: pathlib.Path = pathlib.Path("alchemical_network.json"),
):
    """
    Create an alchemical network.
    """
    smcs = load_ligands_to_dict(input_ligands)
    logger.info(f"{len(smcs)=}")

    chemical_systems = get_chemical_systems(
        smcs,
        csv_file,
    )
    logger.info(f"{len(chemical_systems)=}")

    # load force field
    if pathlib.Path(forcefield_file).exists():
        forcefield = pathlib.Path(forcefield_file).read_text()
        logger.info(f"Read forcefield from {forcefield_file}")
    else:
        forcefield = forcefield_file
        logger.info(f"Using built-in forcefield {forcefield_file}")

    transformations = []
    for chemical_system in chemical_systems:
        transformations.append(get_transformation(chemical_system, forcefield))
    logger.info(f"{len(transformations)=}")

    alchemical_network = AlchemicalNetwork(transformations)
    alchemical_network.to_json(output_file)
    logger.info(f"Wrote alchemical network with {len(transformations)} transformations to {output_file}")


if __name__ == "__main__":
    main()