import tqdm

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.compute as pc

from loguru import logger

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
    dataset = ds.dataset("../../04_benchmark/qm/improper-values/")

    # filter
    counts = dataset.to_table(
        columns=["qcarchive_id", "method"]
    ).to_pandas().groupby(
        by=["qcarchive_id", "method"]
    ).first().reset_index().groupby("qcarchive_id").count().reset_index()

    filtered_df = dataset.filter(
        pc.field("qcarchive_id").isin(
            counts[counts.method == 4].qcarchive_id.values
        )
    ).to_table().to_pandas()

    filtered_df["atom_indices_str"] = [
        "-".join(map(str, x))
        for x in tqdm.tqdm(filtered_df.atom_indices.values)
    ]

    others = filtered_df[filtered_df.method != "qm"]

    topology_rows = []
    qm = filtered_df[filtered_df.method == "qm"].sort_values(by=["qcarchive_id", "atom_indices_str"])
    for method, subdf in others.groupby("method"):
        subdf = pd.DataFrame(subdf.sort_values(by=["qcarchive_id", "atom_indices_str"]))
        assert len(subdf) == len(qm), f"Length mismatch for impropers -- {method}, {len(subdf)} vs {len(qm)}"

        subdf["reference_value"] = qm["value"].values

        for qca_id, smidf in tqdm.tqdm(subdf.groupby("qcarchive_id"), desc=f"Impropers -- {method}"):
            diff = []
            for x, y in zip(smidf.value, smidf.reference_value):
                diff.append(compute_angle_difference(x, y))
            if not len(diff):
                continue
            diff = np.array(diff)
            rmse = (diff ** 2).mean() ** 0.5
            row = {
                "qcarchive_id": qca_id,
                "mapped_smiles": smidf.mapped_smiles.values[0],
                "FF": METHODS[method],
                "Property": f"ImproperTorsions ICRMSD",
                "Value": rmse,
            }
            topology_rows.append(row)
    topology_df = pd.DataFrame(topology_rows)
    topology_df.to_csv("processed_improper_torsions.csv", index=False)


if __name__ == "__main__":
    main()
