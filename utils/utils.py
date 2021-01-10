from multiprocessing import Pool, Process
from pathlib import Path
from typing import Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import pandas as pd
import numpy as np

from airflow import DAG
from airflow.operators.python_operator import PythonOperator
from pandas import DataFrame

from consts.biomart import BIOMART_DATA_PATH
from utils.logger import logger
import subprocess

def get_wrapper(func, *columns, **kwargs):
    def wrapper(row):
        return func(*[row[c] for c in columns], **kwargs)
    return wrapper

#
# def parallelize_apply(df, func, n_cores=4):
#     df_split = np.array_split(df, n_cores)
#     pool = Pool(n_cores)
#     df = pd.concat(pool.map(func, df_split))
#     pool.close()
#     pool.join()
#     return df


def get_subsequence_by_coordinates(full_sequence: str, start: int, end: int, strand="+",
                                   extra_chars: int = 0) -> str:
    if start == -1:
        return  f"Error:no subsequence to extract. start={start}, end={end}"

    start = int(start)
    end = int(end)

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
    logger.info(f"read file {str(path)}")
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
    with fasta_filename.open() as fasta:
        logger.info(f"read fasta file {fasta_filename}")
        d = [{"ID": seq_record.id,
              "sequence": str(seq_record.seq)}
             for seq_record in SeqIO.parse(fasta, 'fasta') if str(seq_record.id).startswith(match)]

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


def concatenate_biomart_df(organism: str):
    df_dict = {
        "utr3" : pd.read_csv(BIOMART_DATA_PATH/f"{organism}_3utr.csv"),
        "utr5": pd.read_csv(BIOMART_DATA_PATH / f"{organism}_5utr.csv"),
        "cds": pd.read_csv(BIOMART_DATA_PATH / f"{organism}_coding.csv")
    }
    for region in df_dict:
        df_dict[region]["region"] = region

    return pd.concat(df_dict.values(), ignore_index=True)

def call_wrapper(cmd: str, cwd: Path):
    # logger.info(f"running {cmd} from {cwd}")
    return subprocess.call(cmd.split(), cwd=cwd.resolve())


def DirectorySpecificBashOperator(task_id: str, cmd: str, dag: DAG, cwd: Path) -> \
        PythonOperator:
    return PythonOperator(
        task_id=task_id,
        python_callable=call_wrapper,
        op_kwargs={"cmd": cmd,
                   'cwd': cwd},
        dag=dag)
