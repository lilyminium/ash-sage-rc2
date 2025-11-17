import pathlib
from helpers import plot_facetgrid, read_into_df, print_latex_stats, FORCEFIELD_ORDER

KWARGS = dict(
    force_field_order=FORCEFIELD_ORDER,
    force_field_col="FF",
    dodge_spacing=0.15,
    height=2.5,
    aspect=0.8,
)


def main():
    image_directory = pathlib.Path("images")
    image_directory.mkdir(exist_ok=True, parents=True)

    base_dir = pathlib.Path("../../05_analysis/output/sfes/")

    # plot freesolv
    df = read_into_df(base_dir / "freesolv" / "stats")
    data_all = df[df.group == "All"]
    n_props = data_all["n"].values[0]
    g = plot_facetgrid(
        data_all,
        unit_str="[kcal/mol]",
        title=f"∆G$_{{aq}}$ (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "sfes-freesolv-validation.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    data_all.to_csv("processed_freesolv_validation_stats.csv", index=False)
    print_latex_stats(data_all)

    # plot mnsol
    df = read_into_df(base_dir / "mnsol" / "stats")
    data_all = df[df.group == "All"]
    n_props = data_all["n"].values[0]
    g = plot_facetgrid(
        data_all,
        unit_str="[kcal/mol]",
        title=f"∆G$_{{non-aq}}$ (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "sfes-mnsol-validation.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    data_all.to_csv("processed_mnsol_validation_stats.csv", index=False)
    print_latex_stats(data_all)

    # plot ketone mnsol solutes
    df = read_into_df(pathlib.Path("../../05_analysis/output/sfes-solutes/") / "mnsol" / "stats")
    data_all = df[df.group == "Ketone"]
    n_props = data_all["n"].values[0]
    g = plot_facetgrid(
        data_all,
        unit_str="[kcal/mol]",
        title=f"∆G$_{{non-aq}}$, Ketone Solutes (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "sfes-mnsol-ketone-solutes-validation.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    data_all.to_csv("processed_mnsol_ketone_solutes_validation_stats.csv", index=False)
    print_latex_stats(data_all)


    # plot tfes
    df = read_into_df(base_dir / "tfes" / "stats")
    data_all = df[df.group == "All"]
    n_props = data_all["n"].values[0]
    g = plot_facetgrid(
        data_all,
        unit_str="[kcal/mol]",
        title=f"∆G$_{{trans}}$ (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "sfes-tfes-validation.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    data_all.to_csv("processed_tfes_validation_stats.csv", index=False)
    print_latex_stats(data_all)

    


if __name__ == "__main__":
    main()
