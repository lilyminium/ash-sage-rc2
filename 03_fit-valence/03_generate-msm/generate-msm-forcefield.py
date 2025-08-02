import collections
import json
import pathlib
import typing

import click
import tqdm

import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds

from openff.toolkit import Molecule, ForceField
from openff.units import unit



@click.command()
@click.option(
    "--input-forcefield",
    "-i",
    type=str,
    help="Input forcefield file to modify.",
)
@click.option(
    "--output-forcefield",
    "-o",
    type=str,
    help="Output forcefield file to write.",
)
@click.option(
    "--output-msm",
    "-om",
    type=str,
    help="Output MSM file to write.",
)
@click.option(
    "--msm-data-directory",
    "-im",
    default="msm-data",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help=(
        "Directory containing the input MSM data. "
        "This is a directory containing parquet files."
    ),
)
@click.option(
    "--aggregator",
    "-a",
    default="mean",
    type=click.Choice(["mean", "median"]),
    help=(
        "Aggregator to use for the MSM data. "
        "This is either 'mean' or 'median'."
    ),
)
def main(
    input_forcefield: str,
    output_msm: str,
    output_forcefield: str,
    msm_data_directory: str = "msm-data",
    aggregator: typing.Literal["mean", "median"] = "mean",
):
    msm_data = ds.dataset(msm_data_directory)
    print(f"Loaded {msm_data.count_rows()} MSM records")

    df = msm_data.to_table().to_pandas()
    n_cmiles = len(df.cmiles.unique())

    ff = ForceField(input_forcefield)

    output_msm = pathlib.Path(output_msm)
    output_msm.parent.mkdir(parents=True, exist_ok=True)
    if False: # output_msm.is_file():
        with output_msm.open("r") as f:
            all_msm_values = json.load(f)
    else:
        # collect all MSM values by parameter
        all_msm_values = collections.defaultdict(
            lambda: collections.defaultdict(
                lambda: collections.defaultdict(list)
            )
        )

        for cmiles, subdf in tqdm.tqdm(
            df.groupby("cmiles"),
            total=n_cmiles,
            desc="Processing molecules to FF parameters"
        ):
            mol = Molecule.from_mapped_smiles(
                cmiles,
                allow_undefined_stereo=True,
            )
            labels = ff.label_molecules(mol.to_topology())[0]
            
            for parameter_type, parameter_df in subdf.groupby("parameter_type"):
                parameter_labels = labels[parameter_type]
                all_parameters = set()
                parameter_type_dict = all_msm_values[parameter_type]

                for _, row in parameter_df.iterrows():
                    indices = tuple(row["indices"])
                    try:
                        parameter = parameter_labels[indices]
                    except KeyError:
                        continue

                    parameter_dict = parameter_type_dict[parameter.id]
                    parameter_dict["eq"].append(row["eq"])
                    parameter_dict["k"].append(row["force_constant"])
                    parameter_dict["id"].append(row["id"])
                    parameter_type_dict[parameter.id]["cmiles"].append(cmiles)
                    all_parameters.add(parameter)

                for parameter in all_parameters:
                    parameter_type_dict[parameter.id]["smirks"] = parameter.smirks
        
        # save output json for debugging
        with output_msm.open("w") as f:
            json.dump(
                all_msm_values,
                f,
                indent=4,
            )
        print(f"Wrote MSM partitions to {output_msm}")

    # aggregate and update FF
    kj_per_mol_per_nm2 = unit.kilojoule_per_mole / (unit.nanometer ** 2)
    kcal_per_mol_per_a2 = unit.kilocalorie_per_mole / (unit.angstrom ** 2)
    kj_per_mol_per_rad2 = unit.kilojoule_per_mole / (unit.radian ** 2)
    kcal_per_mol_per_rad2 = unit.kilocalorie_per_mole / (unit.radian ** 2)

    if aggregator == "mean":
        agg_func = np.mean
    elif aggregator == "median":
        agg_func = np.median
    else:
        raise ValueError(
            f"Aggregator must be 'mean' or 'median', not {aggregator}"
        )
    print(f"Aggregating with {aggregator}")

    for parameter_type, parameter_type_dict in all_msm_values.items():
        handler = ff.get_parameter_handler(parameter_type)
        for parameter_id, parameter_dict in parameter_type_dict.items():
            parameter = handler[parameter_dict["smirks"]]
            k = agg_func(parameter_dict["k"])
            eq = agg_func(parameter_dict["eq"])

            if parameter_type == "Bonds":
                k = k * kj_per_mol_per_nm2
                eq = eq * unit.nanometer

                parameter.length = eq.to(unit.angstrom)
                parameter.k = k.to(kcal_per_mol_per_a2)

            elif parameter_type == "Angles":
                k = k * kj_per_mol_per_rad2
                eq = eq * unit.radian
                
                parameter.k = k.to(kcal_per_mol_per_rad2)

                if np.isclose(parameter.angle.m, 180.0, atol=0, rtol=0):
                    # this is linear, ignore
                    continue
                parameter.angle = eq.to(unit.degree)


    ff.to_file(output_forcefield)
    print(f"Wrote forcefield to {output_forcefield}")



if __name__ == "__main__":
    main()
