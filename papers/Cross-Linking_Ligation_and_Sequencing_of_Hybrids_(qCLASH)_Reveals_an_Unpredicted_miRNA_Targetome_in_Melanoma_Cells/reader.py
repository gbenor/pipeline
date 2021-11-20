from typing import Tuple

import pandas as pd
from numpy import int64
from pandas import DataFrame
import click

from consts.global_consts import HUMAN_SITE_EXTENDED_LEN, ROOT_PATH
from consts.pipeline_steps import READ_PATH
from utils.logger import logger
from utils.utils import get_subsequence_by_coordinates, get_wrapper, to_csv



def read() -> DataFrame:
    # Consts
    file_name = ROOT_PATH / "papers" / \
                "Cross-Linking_Ligation_and_Sequencing_of_Hybrids_(qCLASH)_Reveals_an_Unpredicted_miRNA_Targetome_in_Melanoma_Cells" / \
                "Table_S2.xlsx"

    sheet_name="miRNA_mRNA"
    skiprows = None
    nrows = 100
    nrows = None
    validation_string = "ReadSequences"
    usecols = ['ReadSequences', 'miRNA',
               'Read_start_5', 'Read_end_5',
               'Read_start_3', 'Read_end_3',
               'Vienna', 'Seed', 'Binding_Region']

    dtype = {'ReadSequences': str,
             'miRNA': str,
             'Read_start_5': int,
             'Read_end_5': int,
             'Read_start_3': int,
             'Read_end_3': int,
             'Vienna': str,
             'Seed': str,
             'Binding_Region': str}

    # Logic
    logger.info(f"Reading file {file_name}")
    df: DataFrame = pd.read_excel(file_name, skiprows=skiprows, nrows=nrows, usecols=usecols,
                                  sheet_name=sheet_name, dtype=dtype)
    assert df.columns[0] == validation_string, f"reader validation error: {df.columns[0]}"

    return df


def chimera_split(chimera: str, mirna_start: int, mirna_end: int, target_start: int, target_end: int) -> Tuple[str, str]:
    # mirna = get_subsequence_by_coordinates(chimera, mirna_start, mirna_end-1) #(mirna-1) since later I am going to extarct the mirna from mirbase and I would like to avoid from overlapping sequence with the target
    # target = get_subsequence_by_coordinates(chimera, max(mirna_end + 1, target_start), target_end)
    mirna = get_subsequence_by_coordinates(chimera, mirna_start, mirna_end)
    target = get_subsequence_by_coordinates(chimera, target_start, target_end)

    return mirna, target


def change_columns_names(in_df: DataFrame) -> DataFrame:
    df = in_df[['miRNA', "mirna_seq_tmp", "site",
                'Vienna', 'Seed', 'Binding_Region']]
    return df.rename(columns={"miRNA": "miRNA ID",
                              "Binding_Region": "paper region"})


def add_meta_data(df: DataFrame) -> DataFrame:
    paper_name = "Mapping_the_Human_miRNA_Interactome_by_CLASH_Reveals_Frequent_Noncanonical_Binding"
    df.insert(0, "paper name", paper_name)
    df.insert(0, "organism", "human")
    df["key"] = df.reset_index().index

    return df


def save(df: DataFrame, file_name: str):
    full_path = READ_PATH / file_name
    to_csv(df, full_path)


@click.command()
@click.argument('out_filename', type=str)
def run(out_filename: str):
    df = read()
    df[["mirna_seq_tmp", "site"]] = df.apply(func=get_wrapper(
        chimera_split, 'ReadSequences', 'Read_start_5', 'Read_end_5',
        'Read_start_3', 'Read_end_3', ),
        axis=1, result_type="expand")

    df[["mirna_seq_tmp", "site"]] = df[["mirna_seq_tmp", "site"]].replace(to_replace='T', value='U', regex=True)

    df = change_columns_names(df)
    df = add_meta_data(df)
    save(df, out_filename)


if __name__ == '__main__':
    run()
