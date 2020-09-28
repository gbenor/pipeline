from pathlib import Path

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
