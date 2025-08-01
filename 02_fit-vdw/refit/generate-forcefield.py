import collections
import logging
import pathlib

import click
import tqdm

from openff.toolkit import ForceField, Molecule
from openff.evaluator.datasets import PhysicalPropertyDataSet

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

@click.command()
@click.option(
    "--input-forcefield",
    "-i",
    default="openff-2.2.1.offxml",
    help="The input force field file to modify",
)
@click.option(
    "--output-forcefield",
    "-o",
    default="forcefield/force-field.offxml",
    help="The output force field file to save",
)
@click.option(
    "--training-set",
    "-t",
    default="../../01_download-data/physprop/final/output/training-set.json",
    help="The training set file to use for parameterization",
)
@click.option(
    "--n-properties",
    "-n",
    default=5,
    help="The minimum number of properties a vdw parameter must be associated with to be trained"
    "in the output force field",
)
def main(
    input_forcefield: str = "openff-2.2.1.offxml",
    output_forcefield: str = "forcefield/force-field.offxml",
    training_set: str = "../../01_download-data/physprop/final/output/training-set.json",
    n_properties: int = 5,
):
    forcefield = ForceField(input_forcefield)

    training_dataset = PhysicalPropertyDataSet.from_json(training_set)
    
    label_counter = collections.Counter()

    for prop in tqdm.tqdm(training_dataset.properties):
        all_vdw_label_ids = set()
        for component in prop.substance.components:
            mol = Molecule.from_smiles(component.smiles, allow_undefined_stereo=True)
            labels = forcefield.label_molecules(mol.to_topology())[0]["vdW"]
            for parameter in labels.values():
                all_vdw_label_ids.add(parameter.id)

        for label_id in all_vdw_label_ids:
            label_counter[label_id] += 1

    vdw_handler = forcefield.get_parameter_handler("vdW")
    for parameter in vdw_handler.parameters:
        property_count = label_counter.get(parameter.id, 0)
        logger.info(
            f"Parameter {parameter.id} {parameter.smirks} has {property_count} properties associated with it."
        )
        if property_count >= n_properties:
            # parameter.add_cosmetic_attribute("parameterize", "epsilon, rmin_half")
            
            # we want to keep epsilon non-zero
            # there are more elegant and continuous ways to do this, but we also want
            # something that renders as nonzero in a reasonable number of decimal places
            # so we arbitrarily set the minimum value to 1e-5
            parameter.add_cosmetic_attribute("parameterize", "rmin_half, constrained_epsilon")
            parameter.add_cosmetic_attribute("constrained_epsilon", parameter.epsilon)
            # I believe we just ignore units here -- see PR below
            # https://github.com/leeping/forcebalance/pull/281
            parameter.add_cosmetic_attribute(
                "parameter_eval",
                f"epsilon=min(1e-5, PRM['vdW/Atom/constrained_epsilon/{parameter.smirks}'])"
            )

            logger.info(f"Training {parameter.id} {parameter.smirks}")


    # replace AM1-BCC with NAGL
    forcefield.deregister_parameter_handler("ToolkitAM1BCC")
    forcefield.get_parameter_handler(
        "ChargeIncrementModel",
        {
            "version": 0.3,
            "number_of_conformers": 1,
            "partial_charge_method": "openff-gnn-am1bcc-0.1.0-rc.3.pt"
        }
    )

    output_forcefield = pathlib.Path(output_forcefield)
    output_forcefield.parent.mkdir(parents=True, exist_ok=True)
    forcefield.to_file(output_forcefield)
    logger.info(f"Saved force field to {output_forcefield}")


if __name__ == "__main__":
    main()
