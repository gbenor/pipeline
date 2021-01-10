import pandas as pd
from numpy import int64
from pandas import DataFrame
import click

from consts.global_consts import ROOT_PATH
from consts.pipeline_steps import READ_PATH
from utils.logger import logger
from utils.utils import get_wrapper, to_csv



def read() -> DataFrame:
    # Consts
    file_name = ROOT_PATH / "papers" / \
                "Pairing_beyond_the_Seed_Supports_MicroRNA_Targeting_Specificity" / \
                "1-s2.0-S1097276516305214-mmc3.xlsx"
    # nrows = 50
    nrows = None

    skiprows = 3
    validation_string = 'chr'

    usecols = ['chr', 'start', 'stop', 'name', 'strand']
    dtype = {'chr': str,
             'start': int,
             'stop': int,
             'strand': str,
             'name': str}

    # Logic
    logger.info(f"Reading file {file_name}")
    df: DataFrame = pd.read_excel(file_name,  nrows=nrows, usecols=usecols,
                            dtype=dtype, skiprows=skiprows)
    assert df.columns[0] == validation_string, f"reader validation error: {df.columns[0]}"
    df["name"] = df["name"].apply(lambda x: x[:x.find("_")])
    return df


def change_columns_names(df: DataFrame) -> DataFrame:
    return df.rename(columns={"name": "miRNA ID",
                              "stop": "end"})



def add_meta_data(df: DataFrame) -> DataFrame:
    paper_name = "Pairing_beyond_the_Seed_Supports_MicroRNA_Targeting_Specificity"
    df.insert(0, "paper name", paper_name)
    df.insert(0, "organism", "celegans")
    df.insert(0, "paper region", "None")

    df["key"] = df.reset_index().index

    return df


def save(df: DataFrame, file_name: str):
    full_path = READ_PATH / file_name
    to_csv(df, full_path)


@click.command()
@click.argument('out_filename', type=str)
def run(out_filename: str):
    df = read()
    df = change_columns_names(df)
    df = add_meta_data(df)
    save(df, out_filename)


if __name__ == '__main__':
    run()
