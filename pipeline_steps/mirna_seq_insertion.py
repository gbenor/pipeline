from pathlib import Path
import pandas as pd
from pandas import DataFrame
import click

from consts.mirna_utils import MIRBASE_FILE
from consts.pipeline_steps import READ_PATH, MIRNA_SEQ_PATH
from utils.logger import logger
from utils.utils import read_csv, to_csv


@click.command()
@click.argument('fname', type=str)
def mirna_seq_insertion(fname: str):
    logger.info(f"Insert mirna sequence to {fname}")

    fin_full_path = READ_PATH / fname
    fout_full_path = MIRNA_SEQ_PATH / fname

    df: DataFrame = read_csv(fin_full_path)
    mirbase_df: DataFrame = pd.read_csv(MIRBASE_FILE, index_col=0,
                                        usecols=["miRNA ID", "miRNA sequence"])
    join_df = df.merge(mirbase_df, how="left", left_on="miRNA ID", right_on="miRNA ID")
    to_csv(join_df, fout_full_path)
    logger.info(f"Finish the mirna sequence insertion to {fname}")


if __name__ == '__main__':
    mirna_seq_insertion()