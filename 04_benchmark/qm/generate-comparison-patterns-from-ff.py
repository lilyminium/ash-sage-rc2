"""
Generate comparison patterns from a force field.
These can be used with `compare-topology-pattern-from-file.py`.
\b
A file is generated in JSON format with a list of dictionaries with the fields:
    - pattern: The SMIRKS pattern
    - output_stem: The output stem for the files generated with this pattern ("<topology_group>-<parameter_id>")
    - topology_group: The topology group (Bonds, Angles, ProperTorsions, ImproperTorsions)
"""

import click
import json
import pathlib
import sys

from openff.toolkit import ForceField
from loguru import logger

logger.remove()
logger.add(sys.stdout)


@click.command()
@click.option(
    "--input-force-field",
    "-i",
    "input_force_field",
    type=str,
    default="openff-2.2.1.offxml",
    help="The path of the force field to generate patterns from.",
    show_default=True,
)
@click.option(
    "--output-file",
    "-o",
    "output_file",
    type=str,
    default="comparison-patterns/patterns-openff-2.2.1.json",
    help="The path to the output JSON file.",
    show_default=True,
)
def main(
    input_force_field: str = "openff-2.2.1.offxml",
    output_file: str = "comparison-patterns/patterns-openff-2.2.1.json",
):
    logger.info(f"Generating comparison patterns from {input_force_field} to {output_file}")
    
    patterns = []

    ff = ForceField(input_force_field)
    parameter_types = ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"]
    for parameter_type in parameter_types:
        handler = ff.get_parameter_handler(parameter_type)
        for param in handler.parameters:
            pattern = {
                "pattern": param.smirks,
                "output_stem": f"{parameter_type}-{param.id}",
                "topology_group": parameter_type
            }
            patterns.append(pattern)
    
    logger.info(f"Found {len(patterns)} patterns")

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(patterns, f, indent=4)
    logger.info(f"Wrote patterns to {output_file}")


if __name__ == "__main__":
    main()
