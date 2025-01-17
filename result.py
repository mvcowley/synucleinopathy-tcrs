import pathlib

import polars as pl

if __name__ == "__main__":
    path = pathlib.Path("out/query/result")
    files = path.glob("*.csv")
    data = {
        file.resolve().__str__().split("/")[-1].split(".")[0]: pl.read_csv(file)
        for file in files
    }
    dropped = {name: df.drop_nulls() for name, df in data.items()}
    print(dropped)
