from collections import namedtuple
from pathlib import Path
from typing import List

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from pandas import DataFrame, Series
import click

from consts.global_consts import SITE_EXTRA_CHARS
from consts.mirna_utils import MIRBASE_FILE
from consts.pipeline_steps import READ_PATH, MIRNA_SEQ_PATH, SITE_PATH
from utils.logger import logger
from utils.utils import get_subsequence_by_coordinates, get_wrapper, read_csv, to_csv

#
# def get_site_from_row_by_coordinates(row: Series) -> str:
#     try:
#         return get_subsequence_by_coordinates(str(row["mRNA sequence"]),
#                                               int(row["chimera_start"]),
#                                               int(row["chimera_end"]),
#                                               extra_chars=SITE_EXTRA_CHARS)
#     except Exception:
#         return "None"



@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def insert_site_by_coordinates(fin: str, fout:str):
    logger.info(f"Insert site to {fin}")
    df: DataFrame = read_csv(fin)
    df["site"] = df.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                           "mRNA sequence", "chimera_start", "chimera_end",
                                           extra_chars=SITE_EXTRA_CHARS),
                          axis=1)

    to_csv(df, fout)
    logger.info(f"finish the site sequence insertion to {fin}")


@click.group()
def cli():
    pass


cli.add_command(insert_site_by_coordinates)


if __name__ == '__main__':
    cli()

