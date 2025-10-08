"""
Calculate statistics for labelled data.
If a force field is specified, filter to that force field only.
\b
This script saves a summary CSV file with the computed statistics.
The CSV contains the columns:
    - force field column (usually 'FF') (str): force field
    - stat (str): statistic in question (RMSE, MUE, MSE, rho, or R2)
    - n (int): the number of observations
    - mle (float): the MLE value
    - mean (float): the mean value
    - low (float): the lower bound of the confidence interval
    - high (float): the upper bound of the confidence interval
    - group (str): the functional group
"""

import click
import pathlib
import sys

import tqdm
import numpy as np
import pandas as pd

from loguru import logger
import sklearn.metrics
import scipy.stats

logger.remove()
logger.add(sys.stdout)

def bootstrap_statistic(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    dy_true: np.ndarray | None = None,
    dy_pred: np.ndarray | None = None,
    ci: float = 0.95,
    statistic: str = "RMSE",
    nbootstrap: int = 1000,
    include_true_uncertainty: bool = False,
    include_pred_uncertainty: bool = False,
) -> dict:
    """Compute mean and confidence intervals of specified statistic.

    Copied from cinnabar.stats.bootstrap_statistic, with additions.

    Parameters
    ----------
    y_true : ndarray with shape (N,)
        True values
    y_pred : ndarray with shape (N,)
        Predicted values
    dy_true : ndarray with shape (N,) or None
        Errors of true values. If None, the values are assumed to have no errors
    dy_pred : ndarray with shape (N,) or None
        Errors of predicted values. If None, the values are assumed to have no errors
    ci : float, optional, default=0.95
        Interval for confidence interval (CI)
    statistic : str
        Statistic, one of ['RMSE', 'MUE', 'R2', 'rho','KTAU','RAE', 'MSE']
    nbootstrap : int, optional, default=1000
        Number of bootstrap samples
    include_true_uncertainty : bool, default False
        whether to account for the uncertainty in y_true when bootstrapping
    include_pred_uncertainty : bool, default False
        whether to account for the uncertainty in y_pred when bootstrapping

    Returns
    -------
    rmse_stats : dict of float
        'mean' : mean RMSE
        'stderr' : standard error
        'low' : low end of CI
        'high' : high end of CI
    """

    def compute_statistic(
        y_true_sample: np.ndarray,
        y_pred_sample: np.ndarray,
        statistic: str
    ) -> float:
        """Compute requested statistic.

        Parameters
        ----------
        y_true : ndarray with shape (N,)
            True values
        y_pred : ndarray with shape (N,)
            Predicted values
        statistic : str
            Statistic, one of ['RMSE', 'MUE', 'R2', 'rho','RAE','KTAU']

        """

        def calc_RAE(y_true_sample: np.ndarray, y_pred_sample: np.ndarray):
            MAE = sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
            mean = np.mean(y_true_sample)
            MAD = np.sum([np.abs(mean - i) for i in y_true_sample]) / float(len(y_true_sample))
            return MAE / MAD

        def calc_RRMSE(y_true_sample: np.ndarray, y_pred_sample: np.ndarray):
            rmse = np.sqrt(sklearn.metrics.mean_squared_error(y_true_sample, y_pred_sample))
            mean_exp = np.mean(y_true_sample)
            mds = np.sum([(mean_exp - i) ** 2 for i in y_true_sample]) / float(len(y_true_sample))
            rrmse = np.sqrt(rmse**2 / mds)
            return rrmse

        if statistic == "RMSE":
            return sklearn.metrics.root_mean_squared_error(y_true_sample, y_pred_sample)
        elif statistic == "MUE":
            return sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
        # mean signed error
        elif statistic == "MSE":
            return np.mean(y_pred_sample - y_true_sample)
        elif statistic == "R2":
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y_true_sample, y_pred_sample)
            return r_value**2
        elif statistic == "rho":
            return scipy.stats.pearsonr(y_true_sample, y_pred_sample)[0]
        elif statistic == "RAE":
            return calc_RAE(y_true_sample, y_pred_sample)
        elif statistic == "KTAU":
            return scipy.stats.kendalltau(y_true_sample, y_pred_sample)[0]
        else:
            raise Exception("unknown statistic '{}'".format(statistic))

    if dy_true is None:
        dy_true = np.zeros_like(y_true)
    if dy_pred is None:
        dy_pred = np.zeros_like(y_pred)

    assert len(y_true) == len(y_pred)
    assert len(y_true) == len(dy_true)
    assert len(y_true) == len(dy_pred)

    sample_size = len(y_true)
    s_n = np.zeros([nbootstrap], np.float64)  # s_n[n] is the statistic computed for bootstrap sample n
    for replicate in range(nbootstrap):
        y_true_sample = np.zeros_like(y_true)
        y_pred_sample = np.zeros_like(y_pred)
        for i, j in enumerate(np.random.choice(np.arange(sample_size), size=[sample_size], replace=True)):
            stddev_true = np.fabs(dy_true[j]) if include_true_uncertainty else 0
            stddev_pred = np.fabs(dy_pred[j]) if include_pred_uncertainty else 0
            y_true_sample[i] = np.random.normal(loc=y_true[j], scale=stddev_true)
            y_pred_sample[i] = np.random.normal(loc=y_pred[j], scale=stddev_pred)
        try:
            s_n[replicate] = compute_statistic(y_true_sample, y_pred_sample, statistic)
        except Exception as e:
            print(y_true_sample, y_pred_sample)
            raise e

    rmse_stats = dict()
    rmse_stats["mle"] = compute_statistic(y_true, y_pred, statistic)
    rmse_stats["stderr"] = np.std(s_n)
    rmse_stats["mean"] = np.mean(s_n)
    # TODO: Is there a canned method to do this?
    s_n = np.sort(s_n)
    low_frac = (1.0 - ci) / 2.0
    high_frac = 1.0 - low_frac
    rmse_stats["low"] = s_n[int(np.floor(nbootstrap * low_frac))]
    rmse_stats["high"] = s_n[int(np.ceil(nbootstrap * high_frac))]

    return rmse_stats



def compute_stats(
    df: pd.DataFrame,
    forcefield_col: str = "FF",
    x_col: str = None,
    y_col: str = None,
    xerr_col: str = None,
    yerr_col: str = None,
) -> pd.DataFrame:
    stat_rows = []

    for ff, subdf in df.groupby(by=forcefield_col):
        y = subdf[y_col].values
        yerr = subdf[yerr_col].values
        x = subdf[x_col].values
        xerr = subdf[xerr_col].values

        for i, stat in enumerate(["RMSE", "MUE", "MSE", "R2", "rho"], 1):
            s = bootstrap_statistic(
                x,
                y,
                xerr,
                yerr,
                statistic=stat,
                include_true_uncertainty=True,
                include_pred_uncertainty=True,
            )
            stat_row = {forcefield_col: ff, "stat": stat, "n": len(subdf)}
            stat_row.update(s)
            stat_rows.append(stat_row)
    return pd.DataFrame(stat_rows)


@click.command(help=__doc__)
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="Input CSV file containing labelled data to calculate stats for.",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    help="Output CSV file for saving computed statistics.",
)
@click.option("--x-col", "-cx", "x_col", type=str, help="Column name for x-axis.")
@click.option("--y-col", "-cy", "y_col", type=str, help="Column name for y-axis.")
@click.option(
    "--xerr-col",
    "-cxe",
    "xerr_col",
    type=str,
    required=False,
    default=None,
    help="Column name for x-error.",
)
@click.option(
    "--yerr-col",
    "-cye",
    "yerr_col",
    type=str,
    required=False,
    default=None,
    help="Column name for y-error.",
)
@click.option(
    "--forcefield-col",
    "-cff",
    "forcefield_col",
    type=str,
    default="FF",
    help="Column name for forcefield.",
)
@click.option(
    "--forcefield",
    "-ff",
    type=str,
    help="Force field to compute for.",
)
def main(
    input_file: str,
    output_file: str,
    x_col: str = None,
    y_col: str = None,
    xerr_col: str = None,
    yerr_col: str = None,
    forcefield_col: str = "FF",
    forcefield: str | None = None,
):
    df = pd.read_csv(input_file)
    df["All"] = True

    logger.info(f"Read {len(df)} rows from {input_file}")

    if forcefield:
        df = df[df[forcefield_col] == forcefield]
        logger.info(f"Filtered to {len(df)} rows with specified force field {forcefield}")

    # restrict only to data points with same number of observations per FF

    # get groups to compute
    index = list(df.columns).index(forcefield_col)
    groups_to_compute = df.columns[index + 1 :].tolist()
    # move 'All' to the start
    try:
        groups_to_compute.remove("All")
    except ValueError:
        pass
    else:
        groups_to_compute.insert(0, "All")
    logger.info(f"Found {len(groups_to_compute)} groups to compute")

    stat_dfs = []
    for group in tqdm.tqdm(groups_to_compute, desc="Computing stats"):
        subdf = df[df[group]]
        if len(subdf) <= 5:
            logger.warning(f"Skipping group {group} with only {len(subdf)} entries")
            continue
        stat_df = compute_stats(
            subdf, x_col=x_col, y_col=y_col, xerr_col=xerr_col, yerr_col=yerr_col
        )
        stat_df["group"] = group
        stat_dfs.append(stat_df)

    # Concatenate all stats DataFrames
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if not (output_file.is_file() and output_file.exists()):
        all_stats_df = pd.concat(stat_dfs, ignore_index=True)
        all_stats_df.to_csv(output_file, index=False)
        logger.info(f"Saved combined stats to {output_file}")


if __name__ == "__main__":
    main()
