import tqdm

import pyarrow.dataset as ds
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt

from rdkit import Chem
from openff.toolkit import Molecule, ForceField

sns.set_context("paper")
sns.set_style("ticks")
mpl.rcParams['font.sans-serif'] = ["muli"]

FFS = {
    "fb-fit-v3-single-mean-k100_unconstrained": "Sage 2.3.0",
    "openff_unconstrained-2.2.1": "Sage 2.2.1",
    "openff_unconstrained-2.1.0": "Sage 2.1.0",
}


def plot_with_handler(
    ax,
    df,
    parameter_handler,
    unit: str = "$\AA"
):
    
    bond_unique = set(df.parameter_id.unique())
    ax = sns.boxplot(
        data=df,
        y="parameter_id",
        order=[x.id for x in parameter_handler.parameters if x.id in bond_unique],
        x="difference",
        hue="FF",
        fliersize=2,
        whis=(5, 95),
        hue_order=["Sage 2.1.0", "Sage 2.2.1", "Sage 2.3.0"]
    )
    ax.axvline(0, ls="--", lw=1, color="gray")
    ax.grid(axis="x", color="gray", alpha=0.3, ls="--", lw=1)
    ax.grid(axis="y", color="gray", alpha=0.3, ls="--", lw=1)
    ax.set_ylabel("Parameter ID")
    ax.set_xlabel(f"MM - QM [{unit}]")
    plt.tight_layout()
    return ax
    
    

def main():
    # load force field to map parameter ids to smirks
    ff = ForceField("../../04_benchmark/forcefields/fb-fit-v3-single-mean-k100_unconstrained.offxml")
    bond_handler = ff.get_parameter_handler("Bonds")
    angle_handler = ff.get_parameter_handler("Angles")
    torsion_handler = ff.get_parameter_handler("ProperTorsions")

    # load dataset and label FF
    dataset = ds.dataset("../../04_benchmark/qm/topology-labeled-by-parameter/")
    df = dataset.to_table().to_pandas()
    df["FF"] = [FFS[x] for x in tqdm.tqdm(df.method.values)]

    # filter to shared observations
    qcarchive_ids = set(df.qcarchive_id.values)
    print(f"Starting with {len(qcarchive_ids)}")
    for ff, subdf in df.groupby("FF"):
        qcarchive_ids &= set(subdf.qcarchive_id.values)
        print(f"Filtered to {len(qcarchive_ids)} with {ff}")

    df = df[df.qcarchive_id.isin(qcarchive_ids)]
    angles = df[df.topology == "Angles"]
    bonds = df[df.topology == "Bonds"]
    torsions = df[df.topology == "ProperTorsions"]


    # plot bonds
    fig, ax = plt.subplots(figsize=(6, 20))
    ax = plot_with_handler(
        ax,
        bonds,
        bond_handler,
        unit="$\AA$"
    )
    plt.savefig("images/bond_mm-vs-qm.png", dpi=300)
    plt.close()

    # plot angles
    fig, ax = plt.subplots(figsize=(6, 16))
    ax = plot_with_handler(
        ax,
        angles,
        angle_handler,
        unit="$\degree$"
    )
    plt.savefig("images/angle_mm-vs-qm.png", dpi=300)
    plt.close()

    # plot torsions
    fig, ax = plt.subplots(figsize=(6, 40))
    ax = plot_with_handler(
        ax,
        torsions,
        torsion_handler,
        unit="$\degree$"
    )
    plt.savefig("images/torsion_mm-vs-qm.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
