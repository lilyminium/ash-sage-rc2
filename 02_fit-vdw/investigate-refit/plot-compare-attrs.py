"""
Plot comparison of attribute changes (e.g., epsilon or rmin_half) for vdw parameters over iterations.
\b
Generates a grid of subplots with one row per SMIRKS pattern, showing:
- Left subplot: bar chart comparing initial, rc1 final, and v3-k100 final values.
- Right subplot: line chart showing the attribute value over iterations for rc1 and v3-k100, with a dashed line for the initial value.
Saves the resulting figure as a PNG file.
"""

import click
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from openff.toolkit import ForceField


@click.command(help=__doc__)
@click.option(
    "--output-dir",
    "-o",
    default="output",
    help="The directory containing the output CSV files from the refit",
)
@click.option(
    "--attrname",
    "-a",
    default="epsilon",
    type=click.Choice(["epsilon", "rmin_half"], case_sensitive=False),
    help="The attribute name to plot (e.g., 'epsilon' or 'rmin_half')",
)
@click.option(
    "--image-directory",
    "-im",
    type=str,
    default="images",
    help="Path to the output directory for images.",
)
def main(
    output_dir: str = "output",
    attrname: str = "epsilon",
    image_directory: str = "images"
):
    parsley = ForceField("openff-1.3.0.offxml")
    parsley_vdw_handler = parsley.get_parameter_handler("vdW")

    output_dir = pathlib.Path(output_dir)

    # Parameter change snapshots
    p_changes = pd.read_csv(output_dir / "parameter-changes.csv")
    rc1_p_changes = pd.read_csv(output_dir / "rc1-parameter-changes.csv")

    # Per-iteration time series
    p_iters = pd.read_csv(output_dir / "parameters-over-iterations.csv")
    rc1_p_iters = pd.read_csv(output_dir / "rc1-parameters-over-iterations.csv")

    # Filter iteration tables to epsilon only (for the right-hand plots)
    def filter_attr(df, attrname):
        if "Attribute" not in df.columns:
            raise ValueError("Iteration CSVs must contain an 'Attribute' column.")
        # Robust compare: some files may capitalize differently
        return df[df["Attribute"].astype(str).str.lower() == attrname.lower()].copy()

    p_eps = filter_attr(p_iters, attrname)
    rc1_eps = filter_attr(rc1_p_iters, attrname)

    # Make sure Iteration is numeric and Value is float
    for frame in (p_eps, rc1_eps):
        frame["Iteration"] = pd.to_numeric(frame["Iteration"], errors="coerce")
        frame["Value"] = pd.to_numeric(frame["Value"], errors="coerce")
        frame.dropna(subset=["Iteration", "Value"], inplace=True)

    # ----------------------------
    # Determine SMIRKS order and n
    # ----------------------------
    # Use the order appearing in parameter-changes.csv as the canonical list
    smirks_list = p_changes["SMIRKS"].astype(str).tolist()
    SMIRKS_to_smirks = {
        row["SMIRKS"]: row["smirks"]
        for _, row in p_changes.iterrows()
    }
    n = len(smirks_list)
    if n == 0:
        raise ValueError("No rows found in parameter-changes.csv")

    # Build fast lookups for initial and finals
    p_changes_idx = p_changes.set_index("SMIRKS")
    rc1_p_changes_idx = rc1_p_changes.set_index("SMIRKS")

    # ----------------------------
    # Plot
    # ----------------------------
    # Create n x 2 subplot grid
    fig, axes = plt.subplots(
        nrows=n, ncols=2,
        gridspec_kw={"width_ratios": [1, 2]},  # make left column thinner for bar charts
        figsize=(5, max(2.5*n, 3.5)),  # scale vertically with n
        constrained_layout=True,
        sharey='row'  # share y-axis per row
    )

    # If n == 1, axes would be 1D arrays; normalize to 2D indexing
    if n == 1:
        axes = np.array([axes])

    bar_colors = {
        "1.3.0": "tab:gray",
        "2.2.1": "tab:blue",
        "2.3.0rc1": "tab:green",
        "v3-k100": "tab:orange",
    }


    for i, smk in enumerate(smirks_list):
        # Right subplot: time series
        ax_line = axes[i, 1]
        # slice series for this SMIRKS
        rc1_series = rc1_eps[rc1_eps["SMIRKS"] == smk].sort_values("Iteration")
        v3_series = p_eps[p_eps["SMIRKS"] == smk].sort_values("Iteration")

        parsley_eps = getattr(parsley_vdw_handler[SMIRKS_to_smirks[smk]], attrname).m
        min_y = min([rc1_series["Value"].min(), v3_series["Value"].min(), parsley_eps])
        max_y = max([rc1_series["Value"].max(), v3_series["Value"].max(), parsley_eps])
        y_range = max_y - min_y
        min_y = min_y - 0.1 * y_range
        max_y = max_y + 0.1 * y_range

        # Plot lines if data exists
        line_handles = []
        if not rc1_series.empty:
            h1, = ax_line.plot(
                rc1_series["Iteration"], rc1_series["Value"],
                label="2.3.0rc1", color=bar_colors["2.3.0rc1"]
            )
            line_handles.append(h1)
        if not v3_series.empty:
            h2, = ax_line.plot(
                v3_series["Iteration"], v3_series["Value"],
                label="v3-k100", color=bar_colors["v3-k100"]
            )
            line_handles.append(h2)

        # Left subplot: bars
        ax_bar = axes[i, 0]
        # Pull values with graceful fallback to NaN if missing
        init_eps = p_changes_idx.loc[smk, f"initial_{attrname}"] if smk in p_changes_idx.index else np.nan
        rc1_final = rc1_p_changes_idx.loc[smk, f"final_{attrname}"] if smk in rc1_p_changes_idx.index else np.nan
        v3_final = p_changes_idx.loc[smk, f"final_{attrname}"] if smk in p_changes_idx.index else np.nan

        x = np.arange(4)
        heights = [parsley_eps, init_eps, rc1_final, v3_final]
        labels = ["1.3.0", "2.2.1", "2.3.0rc1", "v3-k100"]
        colors = [bar_colors[lbl] for lbl in labels]

        # Plot bars; if a value is NaN, Matplotlib will draw an empty bar (height 0)
        ax_bar.bar(x, heights, color=colors)
        ax_bar.set_xticks(x, labels, rotation=90)
        ax_bar.set_ylabel(attrname)
        ax_bar.set_title(smk, loc="left", fontsize=10)

        # Dashed initial-epsilon reference
        if pd.notna(init_eps):
            h3 = ax_line.axhline(
                init_eps, linestyle="--", linewidth=1.0,
                color=bar_colors["2.2.1"], label="2.2.1"
            )
            # Include in legend only once per row
            line_handles.append(h3)

        ax_line.set_xlabel("Iteration")
        ax_line.set_ylabel(attrname)
        # Put legend only if we plotted something
        if line_handles:
            # Avoid duplicate labels on repeated artists
            handles, lbls = ax_line.get_legend_handles_labels()
            uniq = dict(zip(lbls, handles))
            ax_line.legend(uniq.values(), uniq.keys(), frameon=False, fontsize=8, loc="best")

        # Sync y-limits
        ylim = (min_y, max_y)
        ax_bar.set_ylim(ylim)
        ax_line.set_ylim(ylim)

    # Nice overall figure title (optional)
    fig.suptitle(f"{attrname} changes by SMIRKS", y=1.02)
    plt.tight_layout()

    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{image_directory}/{attrname}_changes_by_SMIRKS.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
