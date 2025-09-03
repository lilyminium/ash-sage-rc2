import pandas as pd

def get_limits(df: pd.DataFrame, x: str, y: str) -> tuple[float, float]:
    """Get sensible scatterplot limits"""
    min_, max_ = min(df[y].values), max(df[y].values)
    min2_, max2_ = min(df[x].values), max(df[x].values)
    min_ = min([min_, min2_])
    max_ = max([max_, max2_])
    range_ = max_ - min_
    inc = range_ / 20
    return min_ - inc, max_ + inc