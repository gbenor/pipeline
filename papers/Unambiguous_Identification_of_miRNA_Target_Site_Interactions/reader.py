from pathlib import Path
from typing import Dict

import pandas as pd
from numpy import int64
from pandas import DataFrame
import click
from consts.global_consts import ROOT_PATH
from consts.pipeline_steps import READ_PATH
from utils.logger import logger
from utils.utils import get_wrapper, to_csv


def read(organism: str) -> DataFrame:
    # Consts
    file_name = ROOT_PATH / "papers" / \
                "Unambiguous_Identification_of_miRNA_Target_Site_Interactions" / \
                "1-s2.0-S1097276514003566-mmc3.xls"

    sheets = {"celegans": "Sheet2",
              "human": "Sheet3",
              "viral": "Sheet4",
              "mouse": "Sheet5"}

    skiprows = 3
    # nrows = 20
    nrows = None

    usecols = ['miRNA ID', 'miRNA sequence', 'target sequence']
    dtype = {'miRNA ID': str,
             'miRNA sequence': str,
             'target sequence': str}

    sheet = sheets[organism]
    logger.info(f"Reading file {file_name}, organism={organism} sheet={sheet}")

    inter_df = pd.read_excel(file_name, sheet_name=sheet, header=None)
    assert_org = "elegans" if "elegans" in organism else organism
    assert assert_org in inter_df.iloc[0, 0].lower(), "Read the wrong sheet. no {} in the first cell".format(organism)

    df: DataFrame = pd.read_excel(file_name, sheet_name=sheet, nrows=nrows, usecols=usecols,
                                  dtype=dtype, skiprows=skiprows)

    seq_cols = ['miRNA sequence', 'target sequence']
    df[seq_cols] = df[seq_cols].replace(to_replace='T', value='U', regex=True)

    return df


def change_columns_names(df: DataFrame) -> DataFrame:
    return df.rename(columns={"miRNA": "miRNA ID",
                              "target sequence": "site",
                              "region": "paper region"})



def add_meta_data(df: DataFrame, organism: str) -> DataFrame:
    paper_name = "Unambiguous Identification of miRNA:Target Site Interactions by Different Types of Ligation Reactions"
    df.insert(0, "paper name", paper_name)
    df.insert(0, "organism", organism)
    df.insert(0, "paper region", "None")

    df["key"] = df.reset_index().index

    return df


def save(df: DataFrame, file_name: str):
    full_path = READ_PATH / file_name
    to_csv(df, full_path)


@click.command()
@click.argument('organism', type=str)
@click.argument('out_filename', type=str)
def run(organism: str, out_filename: str):
    df: DataFrame = read(organism)
    df = change_columns_names(df)
    df = add_meta_data(df, organism)
    save(df, out_filename)


if __name__ == '__main__':
    run()
