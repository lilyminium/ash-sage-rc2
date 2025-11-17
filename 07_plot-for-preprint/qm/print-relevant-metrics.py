import pandas as pd


def main():
    df = pd.read_csv("processed_bulk_qm_properties.csv")

    ddes = df[df.Property == "|ddE|"]
    for ff, ffdf in ddes.groupby("FF"):
        mean_ddE = ffdf["Value"].mean()
        median_ddE = ffdf["Value"].median()
        print(f"{ff} |ddE|: n={len(ffdf)}, mean={mean_ddE:.2f} kcal/mol, median={median_ddE:.2f} kcal/mol")


if __name__ == "__main__":
    main()