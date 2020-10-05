from pathlib import Path
from typing import Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import pandas as pd
from pandas import DataFrame

from utils.logger import logger


def get_wrapper(func, *columns, **kwargs):
    def wrapper(row):
        return func(*[row[c] for c in columns], **kwargs)
    return wrapper


def get_subsequence_by_coordinates(full_sequence: str, start: int, end: int, strand="+",
                                   extra_chars: int = 0) -> str:
    assert strand=='+' or strand=='-', "strand value incorrect {}".format(strand)
    strand = 1 if strand == '+' else -1
    start = max(0, start - 1 - extra_chars) # because python is zero based
    end = min(len(full_sequence), end + extra_chars)
    if (end - start) <= 0:
        return  f"Error:no subsequence to extract. start={start}, end={end}, full_seq_len={len(full_sequence)}"
    seq = Seq(full_sequence)
    feature_loc: FeatureLocation = FeatureLocation(start, end, strand=strand)
    sub_sequence: Seq = feature_loc.extract(seq)
    return str(sub_sequence)


def to_csv(df, path: Path) -> None:
    df.reset_index(inplace=True)
    df.loc[-1] = df.dtypes
    df.index = df.index + 1
    df.sort_index(inplace=True)
    df.to_csv(path, index=False)
    logger.info(f"Saved file {path}")


def read_csv(path: Path) -> DataFrame:
    logger.info(f"read file {path}")
    # Read types first line of csv
    dtypes = pd.read_csv(path, nrows=1).iloc[0].to_dict()
    # Read the rest of the lines with the types from above
    return pd.read_csv(path, dtype=dtypes, skiprows=[1], index_col=0)

def get_substring_index(full: str, substring:str) -> Tuple[int, int]:
    start = full.find(substring)
    if start == -1:
        return -1, -1
    return start + 1, start+len(substring)


def fasta_to_dataframe(fasta_filename: Path, match: str = "") -> DataFrame:
    # df: DataFrame = DataFrame(columns=["ID", "sequence"])


    with fasta_filename.open() as fasta:
        logger.info(f"read fasta file {fasta_filename}")
        d = [{"ID": seq_record.id,
          "sequence": str(seq_record.seq)}
            for seq_record in SeqIO.parse(fasta, 'fasta') if str(seq_record.id).startswith(match)]


        #
        # cnt = 0
        # for seq_record in SeqIO.parse(fasta, 'fasta'):  # (generator)
        #     if str(seq_record.id).startswith(match):
        #         df.iloc[cnt] = [seq_record.id, str(seq_record.seq)]
        #         # df = df.append(
        #         #     pd.Series([seq_record.id, str(seq_record.seq)], index=df.columns), ignore_index=True)
        #         cnt += 1
        #         print(df)


    logger.info(f"read {len(d)} fasta sequences")
    return pd.DataFrame(d)



def filter_Sequenceunavailable_from_fasta(fasta_file: Path):
    fasta_sequences = SeqIO.parse(fasta_file.open(), 'fasta')
    valid_seq = []
    all_cnt = 0
    valid_cnt = 0
    for s in fasta_sequences:
        all_cnt += 1
        if s.seq != 'Sequenceunavailable':
            s.id=s.id[:50] # max len for blastdb id
            valid_seq.append(s)
            valid_cnt += 1

    SeqIO.write(valid_seq, fasta_file, 'fasta')
    logger.info(f"filter_Sequenceunavailable_from_fasta \n"
                f"all seq={all_cnt} \n"
                f"valid seq={valid_cnt}")