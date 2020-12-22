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


def read() -> Dict[str, DataFrame]:
    # Consts
    file_name = ROOT_PATH / "papers" / \
                "Unambiguous_Identification_of_miRNA_Target_Site_Interactions" \ 
                "1-s2.0-S1097276514003566-mmc3.xls"

    sheets = {"celegans": "Sheet2",
              "human": "Sheet3",
              "viral": "Sheet4",
              "mouse": "Sheet5"}

    skiprows = 3
    nrows = 50
    # nrows = None

    usecols = ['miRNA ID', 'miRNA sequence', 'target sequence']
    dtype = {'miRNA ID': str,
             'miRNA sequence': str,
             'target': str}

    result: Dict[str, DataFrame] = {}
    for organism, sheet in sheets.items():
        logger.info(f"Reading file {file_name}, organism={organism} sheet={sheet}")

        inter_df = pd.read_excel(file_name, sheet_name=current_sheet, header=None)
        assert organism in inter_df.iloc[0, 0].lower(), "Read the wrong sheet. no {} in the first cell".format(organism)

        df: DataFrame = pd.read_excel(file_name, sheet_name=current_sheet, nrows=nrows, usecols=usecols,
                                      dtype=dtype, skiprows=skiprows)

        result[organism] = df

    return result


def change_columns_names(df: DataFrame) -> DataFrame:
    return df.rename(columns={"miRNA": "miRNA ID",
                              "region": "paper region"})



def add_meta_data(df: DataFrame, organism: str) -> DataFrame:
    paper_name = "Unambiguous Identification of miRNA:Target Site Interactions by Different Types of Ligation Reactions"
    df.insert(0, "paper name", paper_name)
    df.insert(0, "organism", organism)

    df["key"] = df.reset_index().index

    return df


def save(df: DataFrame, file_name: str):
    full_path = READ_PATH / file_name
    to_csv(df, full_path)


@click.command()
@click.argument('out_filename', type=str)
def run(out_filename: str):
    dfdict: Dict[str, DataFrame] = read()
    for organism in dfdict:
        dfdict[organism] = change_columns_names(dfdict[organism])
        dfdict[organism] = add_meta_data(dfdict[organism], organism)
        save(dfdict[organism], out_filename.replace(".csv", f"{organism}.csv"))


if __name__ == '__main__':
    run()
