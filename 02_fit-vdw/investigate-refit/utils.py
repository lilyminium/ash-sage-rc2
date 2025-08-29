"""This module contains common functions I use often for plotting and analysis"""

def get_limits(subdf, y) -> tuple[float, float]:
    """
    Get sensible axis limits for a square scatter plot based on the distribution
    of the data. This assumes that the reference column is called "Reference".

    Parameters
    ----------
    subdf : pd.DataFrame
        The subset of the DataFrame containing the data to plot.
    y : str
        The name of the target column to use for the y-axis.

    Returns
    -------
    tuple[float, float]
        Returns the lower and upper values of the axis limit.
    """
    min_, max_ = min(subdf[y].values), max(subdf[y].values)
    min2_, max2_ = min(subdf["Reference"].values), max(subdf["Reference"].values)
    min_ = min([min_, min2_])
    max_ = max([max_, max2_])
    range_ = max_ - min_
    inc = range_ / 20
    return min_ - inc, max_ + inc