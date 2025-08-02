

def get_limits(subdf, y):
    min_, max_ = min(subdf[y].values), max(subdf[y].values)
    min2_, max2_ = min(subdf["Reference"].values), max(subdf["Reference"].values)
    min_ = min([min_, min2_])
    max_ = max([max_, max2_])
    range_ = max_ - min_
    inc = range_ / 20
    return min_ - inc, max_ + inc