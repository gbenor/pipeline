from pathlib import Path
from typing import List
from functools import reduce

import click
import pandas as pd
from pandas import DataFrame, Series

from consts.global_consts import SITE_EXTRA_CHARS
from consts.pipeline_steps import GAMBIAE_INFORMATION_COLUMN_NAMES, GAMBIAE_INFORMATION_DTYPE, \
    GAMBIAE_INFORMATION_FILENAME, \
    GAMBIAE_INFORMATION_USECOLS, NORMALIZATION_COLUMNS, REGION_PATH, SITE_PATH
from utils.logger import logger
from utils.utils import get_subsequence_by_coordinates, get_subsequence_by_coordinates_no_exception, \
    get_substring_index, get_wrapper, read_csv, to_csv

def extract_seed_family(m: str) -> str:
    try:
        return m[1:7]
    except TypeError:
        return ""


@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def finalize(fin: str, fout:str):
    df: DataFrame = read_csv(Path(fin))

    logger.info("extract the site")
    df["site"] = df[df["sequence"].notnull()].apply(
        func=get_wrapper(get_subsequence_by_coordinates_no_exception,
                         "sequence", "start", "end",
                         extra_chars=SITE_EXTRA_CHARS),
        axis=1)


    def eta(x):
        try:
            return int(x) - SITE_EXTRA_CHARS
        except Exception:
            print(x)
            raise Exception()

    df["start"] = df[df["start"].notnull()]["start"].apply(lambda x: int(x) - SITE_EXTRA_CHARS if int(x) > SITE_EXTRA_CHARS else 1)
    df["end"] = df[df["end"].notnull()]["end"].apply(lambda x: int(x) + SITE_EXTRA_CHARS)


    logger.info("replace T with U")
    seq_cols = ['miRNA sequence', 'site', 'sequence']
    df[seq_cols] = df[seq_cols].replace(to_replace='T', value='U', regex=True)

    logger.info("Add seed family")
    df["seed_family"] = df['miRNA sequence'].apply(extract_seed_family)

    logger.info("Add valid/invalid flag")
    invalid_conditions = [pd.isna(df["miRNA sequence"]),
                          pd.isna(df["site"]),
                          df["miRNA sequence"].str.contains('X'),
                          df["miRNA sequence"].str.contains('N'),
                          df["site"].str.contains("N"),
                          df["site"].str.contains("Error"),
                          df["sequence"].str.contains('N'),
                          df["sequence"].str.contains('X'),
                          df["sequence"].str.contains("Error"),
                          df["sequence"].str.contains("None")]
    df["valid_row"] = ~reduce((lambda x, y: x | y), invalid_conditions)

    df = df[NORMALIZATION_COLUMNS]
    to_csv(df, Path(fout))

 # df["start_end"] = df.apply(
    #     func=get_wrapper(get_substring_index,
    #                      "region_sequence", "site"),
    #     axis=1)
    # df[['start', 'end']] = pd.DataFrame([*df["start_end"]], df.index)


@click.group()
def cli():
    pass


cli.add_command(finalize)



if __name__ == '__main__':
    cli()
