import pandas as pd
from numpy import int64
from pandas import DataFrame
import click

from consts.global_consts import ROOT_PATH
from consts.pipeline_steps import READ_PATH
from utils.logger import logger
from utils.utils import to_csv


def read() -> DataFrame:
    # Consts
    file_name = ROOT_PATH / "papers" / \
        "Dynamic_miRNA-mRNA_interactions_coordinate_gene_expression_in_adult_Anopheles_gambiae" / \
        "journal.pgen.1008765.s018.xlsx"
    skiprows = 2
    nrows = 11483
    validation_string = "miRNA ID"
    usecols = ['miRNA ID', 'Transcript ID', 'chimera_start', 'chimera_end']
    dtype = {'miRNA ID': str,
             'Transcript ID': str,
             'chimera_start': int64,
             'chimera_end': int64}

    # Logic
    logger.info(f"Reading file {file_name}")

    df: DataFrame = pd.read_excel(file_name, skiprows=skiprows, nrows=nrows, usecols=usecols, dtype=dtype)
    assert df.columns[0] == validation_string, "reader validation error"

    return df


def change_columns_names(df: DataFrame) -> DataFrame:
    return df.rename(columns={"Transcript ID": "mRNA ID"})


def add_meta_data(df: DataFrame) -> DataFrame:
    paper_name = "Dynamic_miRNA-mRNA_interactions_coordinate_gene_expression_in_adult_Anopheles_gambiae"
    df.insert(0, "paper name", paper_name)
    df.insert(0, "organism", "Anopheles_gambiae")
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
