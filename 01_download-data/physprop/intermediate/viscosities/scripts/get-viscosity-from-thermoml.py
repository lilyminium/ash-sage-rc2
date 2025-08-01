"""
Obtain experimental viscosity values from ThermoML.

This script reads ThermoML XML files and extracts viscosity data.
It saves a number of files:

* `thermoml-with-viscosities.csv`: The full dataset with all properties supported by Evaluator, and the viscosity data.
* `viscosities.csv`: The subset of the dataset with only viscosity data.
* `viscosities-subset.csv`: A subset of *single-component* viscosity data in a specific temperature and pressure range (295-305 K, 100-102 kPa).
* `viscosity-stats.csv`: A summary of the viscosity data, including max, min, mean, standard deviation, and count values of each SMILES
* `viscosity-stats.json`: The same summary as `viscosity-stats.csv`, but in JSON format.
"""

from collections import defaultdict
import json
import pathlib
import tqdm
import multiprocessing

import click
import numpy as np
import pandas as pd

from openff.evaluator.plugins import register_default_plugins

from viscosity_stub import _process_archive

register_default_plugins()


@click.command
@click.option(
    "--input-directory",
    "-i",
    default="ThermoML.v2020-09-30",
    help="The directory containing the ThermoML XML files",
)
@click.option(
    "--output-directory",
    "-o",
    default="viscosities",
    help="The directory to save the output CSV files",
)
@click.option(
    "--n-processes",
    "-np",
    default=1,
    help="The number of processes to use for loading the data",
)
def main(
    input_directory: str = "ThermoML.v2020-09-30",
    output_directory: str = "viscosities",
    n_processes: int = 1,
):
    input_directory = pathlib.Path(input_directory)
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    xmls = sorted(input_directory.glob("*/*.xml"))
    print(f"Found {len(xmls)} XML files")

    # do an initial filter for viscosity
    viscosity_xmls = []
    for xml in tqdm.tqdm(xmls, desc="Initial filter for viscosity data"):
        with open(xml, "r") as file:
            try:
                if "Viscosity" in file.read():
                    viscosity_xmls.append(xml)
            except:
                # occasionally there's a non utf-8 friendly character
                # so just add it to the list anyway
                viscosity_xmls.append(xml)

    print(f"Found {len(viscosity_xmls)} XML files with viscosity data")
    
    # load the data
    with multiprocessing.Pool(n_processes) as pool:
      dfs = list(
          tqdm.tqdm(
              pool.imap(_process_archive, viscosity_xmls),
              total=len(viscosity_xmls)
            )
        )
    
    df = pd.concat(dfs)
    thermoml_with_viscosity_path = output_directory / "thermoml-with-viscosities.csv"
    df.to_csv(thermoml_with_viscosity_path)
    print(
        f"Loaded {len(df)} data points "
        f"and saved to {thermoml_with_viscosity_path}"
    )

    viscosity = df[df['Viscosity Value (Pa * s)'].notna()]

    cols = [
        'Id', 'Temperature (K)', 'Pressure (kPa)', 'Phase', 'N Components',
        'Component 1', 'Role 1', 'Mole Fraction 1', 'Exact Amount 1',
        'Component 2', 'Role 2', 'Mole Fraction 2', 'Exact Amount 2',
        'Component 3', 'Role 3', 'Mole Fraction 3', 'Exact Amount 3',
        'Viscosity Value (Pa * s)', 'Viscosity Uncertainty (Pa * s)', 'Source',
    ]
    viscosity_path = output_directory /"viscosities.csv"
    viscosity[cols].to_csv(viscosity_path)
    print(
        f"Found {len(viscosity)} data points with viscosity data "
        f"and saved to {viscosity_path}"
    )


    # filter for general interest
    viscosity_subdf = viscosity[
        (viscosity['Temperature (K)'] < 305)
        & (viscosity['Temperature (K)'] > 295)
        & (viscosity['Pressure (kPa)'] > 100)
        & (viscosity['Pressure (kPa)'] < 102)
        & (viscosity["N Components"] == 1)
    ]
    viscosity_subset_path = output_directory / "viscosities-subset.csv"

    viscosity_subdf[cols].to_csv(viscosity_subset_path)
    print(
        f"Found {len(viscosity_subdf)} data points "
        "with viscosity data in the range of interest "
        f"and saved to {viscosity_subset_path}"
    )

    # get a JSON file with max/min/mean stats
    viscosities_by_component = defaultdict(list)
    for _, row in viscosity.iterrows():
        # special case water -- water is ~1 centipoise, i.e. 1e-3 Pa.s
        # and apparently typo-ed as ~1 Pa.s
        if row["Component 1"] == "O" and row["Viscosity Value (Pa * s)"] >= 1.0:
            continue
        component = row[f"Component 1"]
        viscosities_by_component[component].append(row["Viscosity Value (Pa * s)"])

    entries = []
    for smi, values in viscosities_by_component.items():
        # special case water, can get some weird values
        if smi == "O":
            max_val = 0.001
        else:
            max_val = max(values)
        entries.append(
            {
                "smiles": smi,
                "max": max_val,
                "min": min(values),
                "mean": np.mean(values),
                "std": np.std(values),
                "n": len(values),
            }
        )
    
    viscosity_stats = pd.DataFrame(entries)
    viscosity_stats_csv_path = output_directory / "viscosity-stats.csv"
    viscosity_stats.to_csv(viscosity_stats_csv_path)
    print(
        f"Saved stats to {viscosity_stats_csv_path}"
    )

    viscosity_stats_by_component = {}
    for _, row in viscosity_stats.iterrows():
        data = dict(row)
        data.pop("smiles")
        viscosity_stats_by_component[row["smiles"]] = data

    viscosity_stats_json_path = output_directory / "viscosity-stats.json"
    with viscosity_stats_json_path.open("w") as file:
        json.dump(viscosity_stats_by_component, file, indent=2)
    print(
        f"Saved stats to {viscosity_stats_json_path}"
    )


if __name__ == "__main__":
    main()
