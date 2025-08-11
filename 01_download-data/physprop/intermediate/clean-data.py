"""
This script cleans the ThermoML data by cleaning or removing known bad points:
1. Swapping mole fractions for certain compounds.
2. Correcting ester names.
3. Filtering out enthalpy of mixing values of 0.
4. Removing data points from specific sources.

It is intended to be run after the initial download of ThermoML data.
"""
from collections import defaultdict
import logging

import pandas as pd

import click

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "--input_file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default="../initial/input/thermoml.csv",
    help=(
        "Input file to clean. "
        "This should be a valid PhysicalPropertyDataSet CSV file containing ThermoML data."
    ),
)
@click.option(
    "--output_file",
    "-o",
    type=click.Path(dir_okay=False, file_okay=True),
    default="intermediate/thermoml-cleaned.csv",
    help=(
        "Output file to save cleaned data. "
        "This will be a valid PhysicalPropertyDataSet CSV file."
    ),
)
def main(
    input_file: str = "../initial/input/thermoml.csv",
    output_file: str = "output/thermoml-cleaned.csv",
):
    input_df = pd.read_csv(input_file)
    logger.info(f"Loaded {input_file} with {len(input_df)} rows")

    # filter out all data points from 10.1016/j.jct.2007.12.002
    # under 1 MPa
    input_df = input_df[input_df["Source"] != "10.1016/j.jct.2007.12.002"]
    logger.info(f"Filtered out 10.1016/j.jct.2007.12.002, now {len(input_df)} rows")


    # swap the mole fractions of all data points in the below
    bad_br_data_df = pd.read_csv("clean-data/bad-br-data.csv")
    bad_p_data_df = pd.read_csv("clean-data/bad-p-data.csv")
    bad_ring_data_df = pd.read_csv("clean-data/bad-ring-data.csv")
    bad_dimethyl_data_df = pd.read_csv("clean-data/bad-dimethyl-carbonate-data.csv")
    mixed_mole_fractions = pd.concat(
        [
            bad_br_data_df,
            bad_p_data_df,
            bad_ring_data_df,
            bad_dimethyl_data_df
        ]
    )
    mixed_mole_fraction_by_source = defaultdict(list)
    for _, row in mixed_mole_fractions.iterrows():
        doi = row["doi"]
        fraction = row["Mole fraction value"]
        mixed_mole_fraction_by_source[doi].append(fraction)

    # for each of the above dataframes, swap the mole fractions
    # for the first and second compound
    seen_mixed_mole_fractions = defaultdict(list)
    new_rows = []
    for _, row in input_df.iterrows():
        doi = row["Source"]
        frac_1 = row["Mole Fraction 1"]
        frac_2 = row["Mole Fraction 2"]
        new_row = dict(row)
        seen = seen_mixed_mole_fractions[doi]
        vals = mixed_mole_fraction_by_source[doi]
        # if frac_1 in seen or frac_2 in seen:
        #     raise ValueError(
        #         f"Already seen this mole fraction pair {(frac_1, frac_2)} for {doi}"
        #     )
        if frac_1 in vals or frac_2 in vals:
            seen_mixed_mole_fractions[doi].extend([frac_1, frac_2])
            new_row["Mole Fraction 1"] = frac_2
            new_row["Mole Fraction 2"] = frac_1
        
        new_rows.append(new_row)

    for doi, fractions in seen_mixed_mole_fractions.items():
        original_fractions = set(mixed_mole_fraction_by_source[doi])
        seen_fractions = set(fractions)
        if not original_fractions.issubset(seen_fractions):
            raise ValueError(
                f"Seen {fractions} fractions for {doi} but expected {original_fractions}"
            )
        # if len(fractions) != len(original_fractions):
        #     raise ValueError(
        #         f"Seen {len(fractions)} fractions for {doi} but expected {len(original_fractions)}:\n"
        #         f"Seen: {fractions}\n"
        #         f"Expected: {original_fractions}"
        #     )

    input_df = pd.DataFrame(new_rows)
    
    # the following have mixed up the the ester namings
    bad_ester_data = pd.read_csv("clean-data/bad-ester-data.csv")
    ESTER_SWITCHES = {
        "CCOC=O": "COC(C)=O",
        "CCCOC=O": "CCC(=O)OC",
        "CCCCOC=O": "CCCC(=O)OC",
        "CCCCCOC=O": "CCCCC(=O)OC",
    }
    ester_sources = bad_ester_data["doi"].unique()
    new_rows = []
    for _, row in input_df.iterrows():
        doi = row["Source"]
        row = dict(row)
        if doi in ester_sources:
            if row["Component 1"] in ESTER_SWITCHES:
                row["Component 1"] = ESTER_SWITCHES[row["Component 1"]]
            if row["Component 2"] in ESTER_SWITCHES:
                row["Component 2"] = ESTER_SWITCHES[row["Component 2"]]
        new_rows.append(row)
    input_df = pd.DataFrame(new_rows)

    # ignore methanol data -- can't find it in paper
    # dubious_methanol_data = pd.read_csv("clean-data/dubious-methanol-data.csv")
    # input_df = input_df[input_df["Source"] != "10.1016/j.tca.2006.02.028"]
    # logger.info(f"Filtered out 10.1016/j.tca.2006.02.028, now {len(input_df)} rows")

    # ignore all dhmix enthalpies of 0, these are likely to be ternary or quarternary systems
    # and this is just a baseline
    dhmix_zeros = input_df[input_df["EnthalpyOfMixing Value (kJ / mol)"] == 0]
    logger.info(f"Found {len(dhmix_zeros)} enthalpy of mixing values of 0")
    input_df = input_df[input_df["EnthalpyOfMixing Value (kJ / mol)"] != 0]
    logger.info(f"Filtered out enthalpy of mixing values of 0, now {len(input_df)} rows")

    input_df.to_csv(output_file, index=False)
    logger.info(f"Saved cleaned data to {output_file}")


if __name__ == "__main__":
    main()