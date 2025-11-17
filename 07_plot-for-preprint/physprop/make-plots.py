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

    base_dir = pathlib.Path("../../05_analysis/output/phys-prop/")

    # plot density validation
    df = read_into_df(base_dir / "density-validation" / "stats")
    density_all = df[df.group == "All"]
    n_props = density_all["n"].values[0]
    g = plot_facetgrid(
        density_all,
        unit_str="[g/mL]",
        title=f"Density (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "physprop-density-validation.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    density_all.to_csv("processed_density_validation_stats.csv", index=False)
    print_latex_stats(density_all)

    # plot density all
    df = read_into_df(base_dir / "density-all" / "stats")
    density_all = df[df.group == "All"]
    n_props = density_all["n"].values[0]
    g = plot_facetgrid(
        density_all,
        unit_str="[g/mL]",
        title=f"Density (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "physprop-density-all.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    density_all.to_csv("processed_density_all_stats.csv", index=False)
    print_latex_stats(density_all)

    # plot dhmix validation
    df = read_into_df(base_dir / "dhmix-validation" / "stats")
    dhmix_all = df[df.group == "All"]
    n_props = dhmix_all["n"].values[0]
    g = plot_facetgrid(
        dhmix_all,
        unit_str="[kJ/mol]",
        title=f"∆Hmix (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "physprop-dhmix-validation.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    dhmix_all.to_csv("processed_dhmix_validation_stats.csv", index=False)
    print_latex_stats(dhmix_all)
    

    # plot dhmix all
    df = read_into_df(base_dir / "dhmix-all" / "stats")
    dhmix_all = df[df.group == "All"]
    n_props = dhmix_all["n"].values[0]
    g = plot_facetgrid(
        dhmix_all,
        unit_str="[kJ/mol]",
        title=f"∆Hmix (n={n_props})",
        **KWARGS
    )
    filename = image_directory / "physprop-dhmix-all.png"
    g.figure.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    dhmix_all.to_csv("processed_dhmix_all_stats.csv", index=False)
    print_latex_stats(dhmix_all)


if __name__ == "__main__":
    main()
