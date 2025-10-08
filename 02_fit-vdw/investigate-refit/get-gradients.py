"""
Get gradients from a force field fit
\b
This script extracts the gradient contributions from each physical property
in a ForceBalance fitting directory over the course of the optimization.
It saves a CSV file with the following columns:
    - Index: The index of the parameter in the ForceBalance optimization.
    - Id: The ID of the parameter in the force field.
    - smirks: The SMIRKS pattern of the parameter.
    - SMIRKS: The SMIRKS pattern of the parameter, with common atom types replaced.
    - Attribute: The attribute of the parameter (e.g., epsilon, rmin_half).
    - pval: The physical gradient contribution.
    - mval: The mathematical gradient contribution.
    - Unit: The unit of the gradient contribution.
    - smiles_1: The SMILES of the first component in the physical property.
    - smiles_2: The SMILES of the second component in the physical property (if applicable).
    - property_type: The type of the physical property (e.g., Density, EnthalpyOfMixing).
    - property_id: The ID of the physical property.
    - diff: The difference between the estimated and reference value of the property.
    - denominator: The denominator used in the objective function for this property type.
    - grad_contribution: The contribution of this gradient to the overall objective function.
    - obj_contribution: The contribution of this property to the overall objective function.
    - Iteration: The iteration number of the optimization.
"""
import os
import pathlib
import sys
import tempfile

import click
import tqdm

import numpy as np
import pandas as pd
from openff.evaluator.forcefield import ParameterGradientKey
from openff.evaluator.datasets import PhysicalPropertyDataSet
from openff.toolkit import ForceField

from forcebalance.nifty import lp_load
from forcebalance.evaluator_io import Evaluator_SMIRNOFF
from forcebalance.parser import parse_inputs
from forcebalance.forcefield import FF
from forcebalance.objective import Objective
from forcebalance.optimizer import Optimizer

from loguru import logger

logger.remove()
logger.add(sys.stdout)


class HidePrinter:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def setup_gradient_keys(evaluator_target):
    """Instead of copying from ForceBalance source code, we copy from the logfile"""

    LOGTEXT = """   0 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#6X4]
   1 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#6X4]
   2 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]
   3 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]
   4 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]
   5 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]
   6 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])(-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]
   7 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])(-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]
   8 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#6X3]
   9 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#6X3]
  10 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]
  11 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]
  12 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#6X3](~[#7,#8,#9,#16,#17,#35])~[#7,#8,#9,#16,#17,#35]
  13 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#6X3](~[#7,#8,#9,#16,#17,#35])~[#7,#8,#9,#16,#17,#35]
  14 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#7]
  15 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#7]
  16 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#1:1]-[#8]
  17 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#1:1]-[#8]
  18 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#6:1]
  19 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#6:1]
  20 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#6X2:1]
  21 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#6X2:1]
  22 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#6X4:1]
  23 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#6X4:1]
  24 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#8:1]
  25 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#8:1]
  26 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#8X2H0+0:1]
  27 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#8X2H0+0:1]
  28 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#8X2H1+0:1]
  29 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#8X2H1+0:1]
  30 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#7:1]
  31 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#7:1]
  32 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#17:1]
  33 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#17:1]
  34 [    vdW/Atom/rmin_half            : 1.00000e+00 ] : vdW/Atom/rmin_half/[#35:1]
  35 [    vdW/Atom/constrained_epsilon  : 1.29714e+01 ] : vdW/Atom/constrained_epsilon/[#35:1]"""
    lines = LOGTEXT.split("\n")
    for line in lines:
        line = line.strip().split()
        if not line:
            continue
        index = int(line[0])
        last_fields = line[-1].strip().split("/")
        attr = last_fields[-2]
        smirks = last_fields[-1]
        if attr == "constrained_epsilon":
            attr = "epsilon"
        parameter_gradient_key = ParameterGradientKey(tag="vdW", smirks=smirks, attribute=attr)
        parameter_value, is_cosmetic = evaluator_target._parameter_value_from_gradient_key(
            parameter_gradient_key
        )
        evaluator_target._gradient_key_mappings[parameter_gradient_key] = index
        evaluator_target._parameter_units[parameter_gradient_key] = parameter_value.units

    # index_counter = 0
    # parameter_gradient_keys = []
    # for field_list in evaluator_target.FF.pfields:
    #     string_key = field_list[0]
    #     key_split = string_key.split("/")

    #     if len(key_split) == 3 and key_split[0] == "":
    #         parameter_tag = key_split[1].strip()
    #         parameter_smirks = None
    #         parameter_attribute = key_split[2].strip()
    #     elif len(key_split) == 4:
    #         parameter_tag = key_split[0].strip()
    #         parameter_smirks = key_split[3].strip()
    #         parameter_attribute = key_split[2].strip()
    #     else:
    #         raise NotImplementedError()

    #     # Use the full attribute name (e.g. k1) for the gradient key.
    #     parameter_gradient_key = ParameterGradientKey(
    #         tag=parameter_tag,
    #         smirks=parameter_smirks,
    #         attribute=parameter_attribute,
    #     )

    #     # Find the unit of the gradient parameter.
    #     parameter_value, is_cosmetic = evaluator_target._parameter_value_from_gradient_key(
    #         parameter_gradient_key
    #     )

    #     if parameter_value is None or is_cosmetic:
    #         # We don't wan't gradients w.r.t. cosmetic parameters.
    #         continue

    #     parameter_unit = parameter_value.units
    #     parameter_gradient_keys.append(parameter_gradient_key)

    #     evaluator_target._gradient_key_mappings[parameter_gradient_key] = index_counter
    #     evaluator_target._parameter_units[parameter_gradient_key] = parameter_unit

    #     index_counter += 1



@click.command(help=__doc__)
@click.option(
    "--input-directory",
    "-i",
    default="../refit",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="Path to the ForceBalance fitting directory.",
)
@click.option(
    "--output-file",
    "-o",
    default="output/gradient-contributions.csv",
    type=click.Path(dir_okay=False, file_okay=True, writable=True),
    help="Path to save the output CSV file with gradient contributions.",
)
def main(
    input_directory: str = "../refit",
    output_file: str = "output/gradient-contributions.csv"
):
    # ForceBalance messes with files.
    # copy the fitting directory to a temp directory and work there.
    fitting_directory = pathlib.Path(input_directory).resolve()
    current_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        tmp_fitting_directory = tmpdir / fitting_directory.name
        os.system(f"cp -r {fitting_directory} {tmp_fitting_directory}")
        os.chdir(tmp_fitting_directory)

        # ForceBalance prints too much, so hide the output
        with HidePrinter():
            options, tgt_opts = parse_inputs("optimize.in")
            forcefield  = FF(options)
            ## The objective function
            objective   = Objective(options, tgt_opts, forcefield)
            evaluator_target = objective.Targets[0]
            setup_gradient_keys(evaluator_target)

        os.chdir(current_dir)

    fitting_directory = pathlib.Path(fitting_directory)
    iter_dir = fitting_directory / "optimize.tmp" / "phys-prop"

    # Get the files
    objective_files = sorted(iter_dir.glob("*/objective.p"))

    # load the training set
    training_set = PhysicalPropertyDataSet.from_json(fitting_directory / "targets/phys-prop/training-set.json")
    TRAINING_SET_BY_ID = {
        prop.id: prop
        for prop in training_set.properties
    }
    n_density = sum(1 for prop in training_set.properties if "Dens" in type(prop).__name__)
    n_dhmix = len(TRAINING_SET_BY_ID) - n_density

    # set up denominators / weights
    DENOMINATORS = {
        "Density": 0.05,
        "EnthalpyOfMixing": 1.6,
    }
    WEIGHTS = {
        "Density": 1 / n_density,
        "EnthalpyOfMixing": 1 / n_dhmix,
    }

    smirnoff_ff = ForceField(
        str(fitting_directory / "forcefield" / "force-field.offxml"),
        allow_cosmetic_attributes=True,
    )
    vdw_handler = smirnoff_ff.get_parameter_handler("vdW")

    all_dfs = []
    for obj_file in tqdm.tqdm(objective_files, desc="Loading gradients"):
        contents = lp_load(obj_file)
        result_file = obj_file.parent / "results.json"
        mval_file = obj_file.parent / "mvals.txt"

        mvals = np.loadtxt(mval_file)
        jacobian = evaluator_target._build_pvals_jacobian(mvals)
        request_result = PhysicalPropertyDataSet.from_json(result_file)

        gradient_rows = []
        grad = np.zeros(len(evaluator_target._gradient_key_mappings))
        for prop in request_result.estimated_properties.properties:
            refprop = TRAINING_SET_BY_ID[prop.id]
            diff = prop.value.m - refprop.value.m
            property_type = type(prop).__name__
            denominator = DENOMINATORS[property_type]
            weight = WEIGHTS[property_type]
            estimated_pvals = np.zeros(len(evaluator_target._gradient_key_mappings))
            # get physical parameters
            for gradkey in prop.gradients:
                index = evaluator_target._gradient_key_mappings[gradkey.key]
                estimated_pvals[index] = gradkey.value.m

            # get mathematical parameters
            estimated_mvals = np.matmul(jacobian, estimated_pvals)
            assert estimated_mvals.shape == mvals.shape, f"Shape mismatch: {estimated_mvals.shape} vs {mvals.shape}"
            gradient_contribution = 2 * weight * diff * estimated_mvals / denominator ** 2
            assert gradient_contribution.shape == grad.shape, f"Shape mismatch: {gradient_contribution.shape} vs {grad.shape}"
            grad += gradient_contribution

            hessian_contribution = (
                2.0
                * weight
                * (np.outer(estimated_mvals, estimated_mvals))
                / denominator ** 2
            )

            # objective function contribution
            obj_contrib = weight * (diff / denominator) ** 2

            for gradkey, index in evaluator_target._gradient_key_mappings.items():
                parameter = vdw_handler[gradkey.smirks]
                
                row = {
                    "Index": index,
                    "Id": parameter.id,
                    "smirks": gradkey.smirks,
                    "SMIRKS": gradkey.smirks.replace("#7,#8,#9,#16,#17,#35", "ENA"),
                    "Attribute": gradkey.attribute,
                    "pval": estimated_pvals[index],
                    "mval": estimated_mvals[index],
                    "smiles_1": prop.substance.components[0].smiles,
                    "smiles_2": prop.substance.components[1].smiles if len(prop.substance.components) > 1 else "",
                    "property_type": property_type,
                    "property_id": prop.id,
                    "diff": diff,
                    "denominator": denominator,
                    "grad_contribution": gradient_contribution[index],
                    "obj_contribution": obj_contrib,
                }
                hess_row = {
                    "hess_contribution_" + str(i): hessian_contribution[index, i]
                    for i in range(len(evaluator_target._gradient_key_mappings))
                }
                row.update(hess_row)
                gradient_rows.append(row)
        iteration_df = pd.DataFrame(gradient_rows)
        iteration_df["Iteration"] = int(obj_file.parent.name.split("_")[-1])

        # sanity check
        obj = iteration_df.groupby("property_id").first().sum()["obj_contribution"]
        assert np.isclose(obj, contents["X"], atol=1e-6), f"Objective function mismatch: {obj} vs {contents['X']}"

        summed = iteration_df.groupby("Index").sum().reset_index()
        grad_contrib_sum = summed.sort_values("Index")["grad_contribution"].values
        assert np.allclose(grad_contrib_sum, contents["G"]), f"Gradient mismatch: {grad_contrib_sum} vs {contents['G']}, {grad}"

        hessian_summed = iteration_df.groupby("Index").sum().reset_index()
        hess_cols = [col for col in hessian_summed.columns if col.startswith("hess_contribution_")]
        hess_contrib_sum = hessian_summed.sort_values("Index")[hess_cols].values
        assert np.allclose(hess_contrib_sum, contents["H"]), f"Hessian mismatch: {hess_contrib_sum} vs {contents['H']}"

        all_dfs.append(iteration_df)
    
    df = pd.concat(all_dfs, ignore_index=True)
    df["Method"] = "vdw-refit"
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)
    logger.info(f"Wrote gradient contributions for {len(objective_files)} iterations to {output_file}")


if __name__ == "__main__":
    main()
