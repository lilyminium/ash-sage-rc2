"""
Label RMSD-TFD results with checkmol functional group information.
"""

import pyarrow.dataset as ds

import pandas as pd

import tqdm
from yammbs.checkmol import analyze_functional_groups, ChemicalEnvironment

FFS = {
    "fb-fit-v3-single-mean-k100_unconstrained": "Sage 2.3.0",
    "openff_unconstrained-2.2.1": "Sage 2.2.1",
    "openff_unconstrained-2.1.0": "Sage 2.1.0",
}

def main(
    threshold: float = 1
):
    dataset = ds.dataset("../../04_benchmark/qm/data/optimization/qm/")
    QCARCHIVE_ID_TO_SMILES = {}
    for row in dataset.to_table(columns=["qcarchive_id", "cmiles"]).to_pylist():
        QCARCHIVE_ID_TO_SMILES[row["qcarchive_id"]] = row["cmiles"]

    rmsd = ds.dataset("../../04_benchmark/qm/rmsd-tfd/")
    df = rmsd.to_table().to_pandas()
    df["FF"] = [FFS[x] for x in df.method.values]

    # label each row
    labelled_rows = []
    base = {x.value: False for x in ChemicalEnvironment}
    for _, row in tqdm.tqdm(df.iterrows()):
        smi = QCARCHIVE_ID_TO_SMILES[row["qcarchive_id"]]
        row2 = {
            "qcarchive_id": row["qcarchive_id"],
            "rmsd": row["rmsd"],
            "FF": row["FF"],
            "mapped_smiles": smi,
        }
        row2.update(base)
        groups = [gp.value for gp in analyze_functional_groups(smi)]
        for gp in groups:
            row2[gp] = True
        labelled_rows.append(row2)

    labelled_df = pd.DataFrame(labelled_rows)

    # save
    labelled_df.to_csv("checkmol-labelled-rmsd.csv")
    print("Saved to checkmol-labelled-rmsd.csv")

    counts = labelled_df.groupby("qcarchive_id").count()
    counts = counts.reset_index()
    labelled_df = labelled_df[labelled_df.qcarchive_id.isin(counts[counts.FF == 3].qcarchive_id.values)]
    labelled_df.to_csv("checkmol-labelled-rmsd-filtered.csv")
    print("Saved to checkmol-labelled-rmsd-filtered.csv")

    # get high-error set
    subdf = labelled_df[
        labelled_df.FF == "Sage 2.2.1"
    ]
    n_total = len(subdf)

    high_error_set = labelled_df[
        labelled_df.rmsd > threshold
    ]
    group_cols = labelled_df.columns[4:]
    counts_in_base = dict(subdf[group_cols].sum())

    representation_rows = []
    for ff, subdf in high_error_set.groupby("FF"):
        n_high = len(subdf)
        for col in tqdm.tqdm(group_cols, desc=ff):
            row = {
                "FF": ff,
                "Group": col,
                "n_base": n_total,
                "n_in_base": counts_in_base[col],
                "n_high": n_high,
                "n_in_high": sum(subdf[col])
            }
            representation_rows.append(row)

    representation_df = pd.DataFrame(representation_rows)
    representation_df["base_ratio"] = representation_df["n_in_base"] / representation_df["n_base"]
    representation_df["high_ratio"] = representation_df["n_in_high"] / representation_df["n_high"]
    representation_df["high_over_base"] = representation_df["high_ratio"] / representation_df["base_ratio"]
    representation_df.to_csv("checkmol-representation.csv")
    print("Saved to checkmol-representation.csv")

if __name__ == "__main__":
    main()
