import click
import sys

from loguru import logger
import pandas as pd

logger.remove()
logger.add(sys.stdout)

@click.command()
@click.option(
    "--input_file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default="../../initial/input/thermoml.csv",
    help=(
        "Input file to identify baseline dhmixes. "
        "This should be a valid PhysicalPropertyDataSet CSV file containing ThermoML data."
    ),
)
@click.option(
    "--output_file",
    "-o",
    type=click.Path(dir_okay=False, file_okay=True),
    default="baseline-dhmixes.csv",
    help=(
        "Output file to save identified baseline dhmixes. "
        "This will be a valid PhysicalPropertyDataSet CSV file."
    ),
)
def main(
    input_file: str = "../../initial/output/thermoml.csv",
    output_file: str = "baseline-dhmixes.csv",
):
    df = pd.read_csv(input_file)
    logger.info(f"Loaded {input_file} with {len(df)} rows")

    # identify dhmixes
    dhmixes = df[pd.notna(df["EnthalpyOfMixing Value (kJ / mol)"])]
    logger.info(f"Identified {len(dhmixes)} rows with EnthalpyOfMixing data")

    is_zero = dhmixes[dhmixes["EnthalpyOfMixing Value (kJ / mol)"] == 0]
    logger.info(f"Identified {len(is_zero)} rows with 0 kJ/mol EnthalpyOfMixing data, indicating baseline values")

    dois = is_zero["Source"].unique()
    logger.info(f"These baseline dhmixes come from {len(dois)} unique DOIs: {dois}")
    cols = [
        "Id","Temperature (K)","Pressure (kPa)","Phase",
        "N Components","Component 1","Role 1","Mole Fraction 1","Exact Amount 1",
        "Component 2","Role 2","Mole Fraction 2","Exact Amount 2",
        "EnthalpyOfMixing Value (kJ / mol)","EnthalpyOfMixing Uncertainty (kJ / mol)",
        "Source"
    ]
    is_zero[cols].to_csv(output_file, index=False)
    logger.info(f"Saved baseline dhmixes to {output_file}")

if __name__ == "__main__":
    main()