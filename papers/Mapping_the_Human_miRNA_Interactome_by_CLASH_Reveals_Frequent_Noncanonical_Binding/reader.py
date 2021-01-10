import pandas as pd
from numpy import int64
from pandas import DataFrame
import click

from consts.global_consts import HUMAN_SITE_EXTENDED_LEN, ROOT_PATH
from consts.pipeline_steps import READ_PATH
from utils.logger import logger
from utils.utils import get_wrapper, to_csv



def read() -> DataFrame:
    # Consts
    file_name = ROOT_PATH / "papers" / \
        "Mapping_the_Human_miRNA_Interactome_by_CLASH_Reveals_Frequent_Noncanonical_Binding" / \
        "1-s2.0-S009286741300439X-mmc1.txt"
    skiprows = 30
    nrows = None
    sep = "\t"
    validation_string = "microRNA_name"
    usecols = ['microRNA_name', 'miRNA_seq', 'mRNA_name', 'mRNA_seq_extended',
               '5\'UTR', 'CDS', '3\'UTR', 'mRNA_start', 'mRNA_end_extended']

    dtype = {'microRNA_name': str,
             'miRNA_seq': str,
             'mRNA_name': str,
             'mRNA_seq_extended': str,
             '5\'UTR': str,
             'CDS': str,
             '3\'UTR': str}

    # Logic
    logger.info(f"Reading file {file_name}")
    df: DataFrame = pd.read_csv(file_name, skiprows=skiprows, nrows=nrows, usecols=usecols,
                           sep=sep, dtype=dtype)
    assert df.columns[0] == validation_string, "reader validation error"

    return df


def change_columns_names(df: DataFrame) -> DataFrame:
    return df.rename(columns={"microRNA_name": "miRNA ID",
                              "miRNA_seq": "miRNA sequence",
                              "mRNA_name": "mRNA ID"})



def add_region(df: DataFrame) -> DataFrame:
    def get_region(utr5: str, cds: str, utr3: str) ->str:
        result = []
        if utr5 == "1":
            result.append("utr5")
        if cds == "1":
            result.append("cds")
        if utr3 == "1":
            result.append("utr3")
        return "+".join(result)

    df["paper region"] = df.apply(func=get_wrapper(get_region,
                                             '5\'UTR', 'CDS', '3\'UTR'),
                            axis=1)
    df.drop(['5\'UTR', 'CDS', '3\'UTR'], axis=1, inplace=True)
    return df

def remove_extended_from_target(df: DataFrame) -> DataFrame:
    df["site"] = df["mRNA_seq_extended"].apply(lambda s: s[:-1*HUMAN_SITE_EXTENDED_LEN])
    return df


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
    df = change_columns_names(df)
    df = remove_extended_from_target(df)
    df = add_region(df)
    df = add_meta_data(df)
    save(df, out_filename)


if __name__ == '__main__':
    run()
