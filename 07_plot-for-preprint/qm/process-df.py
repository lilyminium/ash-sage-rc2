import pathlib

import pandas as pd
import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds

import seaborn as sns
import tqdm
from matplotlib import pyplot as plt

from openff.toolkit import Molecule


FFS = {
    "Sage 2.1.0": "Sage 2.1.0",
    "Sage 2.2.1": "Sage 2.2.1",
    "v3-k100": "Sage 2.3.0",
}
METHODS = {
    "fb-fit-v3-single-mean-k100_unconstrained": "Sage 2.3.0",
    "openff_unconstrained-2.1.0": "Sage 2.1.0",
    "openff_unconstrained-2.2.1": "Sage 2.2.1",
}

def compute_angle_difference(mm_value, qm_value):
    difference = mm_value - qm_value
    if difference > 180:
        difference -= 360
    if difference < -180:
        difference += 360
    return difference


def main():
    QM_DIR = pathlib.Path("../../04_benchmark/qm/")
    dde_dataset = ds.dataset(QM_DIR / "ddes").filter(
        pc.field("Force field").isin(list(FFS.keys()))
    )
    # subset of the data with just FFs of interest
    rmsds_dataset = ds.dataset(QM_DIR / "rmsd-tfd")

    ic_dataset = ds.dataset(QM_DIR / "topology-values")

    dde_df = dde_dataset.to_table().to_pandas()
    dde_df["FF"] = [FFS[x] for x in dde_df["Force field"].values]

    rmsds_df = rmsds_dataset.to_table().to_pandas()
    rmsds_df["FF"] = [METHODS[x] for x in rmsds_df.method.values]

    # filter to only molecules with all FFS present
    counts = rmsds_df.groupby("qcarchive_id").count()
    matching_indices = counts[counts.FF == len(FFS)].index
    dde_df = dde_df[dde_df.ff_qcarchive_id.isin(matching_indices)]
    rmsds_df = rmsds_df[rmsds_df.qcarchive_id.isin(matching_indices)]
    ic_dataset = ic_dataset.filter(pc.field("qcarchive_id").isin(matching_indices))

    # rename and combine bulk qm properties
    dde_df = dde_df.melt(
        id_vars=["ff_qcarchive_id", "FF"],
        value_vars=["ddE"],
        var_name="Property",
        value_name="Value",
    ).reset_index().rename(columns={"ff_qcarchive_id": "qcarchive_id"})
    abs_dde_df = pd.DataFrame(dde_df)
    abs_dde_df["Value"] = abs_dde_df["Value"].abs()
    abs_dde_df["Property"] = "|ddE|"
    bulk_qm_df = pd.concat([
        rmsds_df.melt(
            id_vars=["qcarchive_id", "FF"],
            value_vars=["rmsd", "tfd"],
            value_name="Value",
            var_name="Property",
        ),
        abs_dde_df,
    ]).reset_index()

    # process internal coordinates
    topologies = ["Bonds", "Angles", "ProperTorsions"]
    topology_rows = []
    QCARCHIVE_TO_MAPPED_SMI = {}
    # do it group by group to avoid memory issues
    for topology in topologies:
        subset = ic_dataset.filter(pc.field("topology") == topology)
        topologydf = subset.to_table().to_pandas()
        topologydf["atom_indices_str"] = ["-".join(map(str, x)) for x in tqdm.tqdm(topologydf.atom_indices.values)]
        
        qm = topologydf[topologydf.method == "qm"].sort_values(by=["qcarchive_id", "atom_indices_str"])
        
        for method, subdf in topologydf[topologydf.method != "qm"].groupby("method"):
            subdf = pd.DataFrame(subdf.sort_values(by=["qcarchive_id", "atom_indices_str"]))
            assert len(subdf) == len(qm), f"Length mismatch for {topology} -- {method}, {len(subdf)} vs {len(qm)}"

            subdf["reference_value"] = qm["value"].values
            for qca_id, smidf in tqdm.tqdm(subdf.groupby("qcarchive_id"), desc=f"{topology} -- {method}"):
                mapped_smi = smidf.mapped_smiles.values[0]
                QCARCHIVE_TO_MAPPED_SMI[qca_id] = mapped_smi
                if topology == "Bonds":
                    diff = smidf.value - smidf.reference_value
                else:
                    diff = []
                    if topology == "ProperTorsions":
                        mol = Molecule.from_mapped_smiles(
                            mapped_smi,
                            allow_undefined_stereo=True
                        )
                        # look for linear torsions
                        linear = []
                        for pattern in [
                            "[*:1]-[*:2]#[*:3]-[*:4]",
                            "[*:1]~[*:2]-[*:3]#[*:4]",
                            "[*:1]~[*:2]=[#6,#7,#16,#15;X2:3]=[*:4]"
                        ]:
                            for match in mol.chemical_environment_matches(pattern):
                                linear.append("-".join(list(map(str, match))))
                                linear.append("-".join(list(map(str, match[::-1]))))
                        smidf = smidf[~smidf.atom_indices_str.isin(linear)]
                    for x, y in zip(smidf.value, smidf.reference_value):
                        # ignore angles too close to linear
                        if abs(compute_angle_difference(y, 180)) < 5:
                            continue
                        else:
                            diff.append(compute_angle_difference(x, y))

                if not len(diff):
                    continue
                diff = np.array(diff)
                rmse = (diff ** 2).mean() ** 0.5
                row = {
                    "qcarchive_id": qca_id,
                    "mapped_smiles": mapped_smi,
                    "FF": METHODS[method],
                    "Property": f"{topology} ICRMSD",
                    "Value": rmse,
                }
                topology_rows.append(row)
                
    topology_df = pd.DataFrame(topology_rows)
    bulk_qm_df["mapped_smiles"] = [
        QCARCHIVE_TO_MAPPED_SMI.get(x, x) for x in bulk_qm_df.qcarchive_id.values
    ]
    columns = ["qcarchive_id", "mapped_smiles", "FF", "Property", "Value"]
    combined_df = pd.concat([bulk_qm_df[columns], topology_df]).reset_index()
    combined_df.to_csv("processed_bulk_qm_properties.csv", index=False)
    print("Wrote processed_bulk_qm_properties.csv")


if __name__ == "__main__":
    main()
