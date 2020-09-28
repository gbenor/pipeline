import pandas as pd
from Bio import SeqIO
from pandas import DataFrame

from consts.mirna_utils import GAMBIAE_FILE, MIRBASE_FILE
from utils.logger import logger


def read_aga_fasta():
    mirbase_df: DataFrame = pd.read_csv(MIRBASE_FILE, index_col=0)

    ver = "aga_paper"
    with GAMBIAE_FILE.open()  as fasta:
        logger.info(f"read fasta file {GAMBIAE_FILE}")
        cnt = 0
        for seq_record in SeqIO.parse(fasta, 'fasta'):  # (generator)
            if str(seq_record.id).startswith("aga-"):
                cnt+=1
                mirbase_df = mirbase_df.append(
                    pd.Series([ver, seq_record.id, str(seq_record.seq)], index=mirbase_df.columns), ignore_index=True)

    mirbase_df.to_csv(MIRBASE_FILE)
    logger.info("update mirbase file")
    logger.info(f"add {cnt} sequences")

if __name__ == '__main__':
    read_aga_fasta()

