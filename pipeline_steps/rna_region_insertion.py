from pathlib import Path
from typing import List

import click
import pandas as pd
from pandas import DataFrame, Series

from consts.pipeline_steps import GAMBIAE_INFORMATION_COLUMN_NAMES, GAMBIAE_INFORMATION_DTYPE, \
    GAMBIAE_INFORMATION_FILENAME, \
    GAMBIAE_INFORMATION_USECOLS, REGION_PATH, SITE_PATH
from utils.logger import logger
from utils.utils import get_wrapper, read_csv


def add_gambiae_region_information(fin: Path) -> DataFrame:
    logger.info(f"enter to add_gambiae_region_information")
    gambiae_region_information: DataFrame = pd.read_csv(GAMBIAE_INFORMATION_FILENAME,
                                                        names=GAMBIAE_INFORMATION_COLUMN_NAMES,
                                                        usecols=GAMBIAE_INFORMATION_USECOLS,
                                                        dtype=GAMBIAE_INFORMATION_DTYPE,
                                                        delimiter="\t")

    in_df: DataFrame = read_csv(fin)
    join_df = in_df.merge(gambiae_region_information, how="left",
                          left_on="mRNA ID",  right_on="TRANSCRIPT_ID")
    logger.info(f"Finish enter the gambiae region information")
    return join_df


def get_region_ranges(len_5utr: int, len_cds: int, len_3utr: int) ->List:
    start = 0
    end = 0
    lens = [len_5utr, len_cds, len_3utr]
    range_list = []
    for i in range(3):
        end += lens[i]
        range_list.append(range(start, end))
        start += lens[i]
    return range_list

def check_trascript_length():


# assert len(transcript) == (len_5utr + len_cds + len_3utr), \
#     f"transcript and region length error: {len(transcript)} != {(len_5utr + len_cds + len_3utr)}"
# assert chimera_start in range(0, len(transcript)), f"chimera_start is not in proper range. " \
#                                                    f"chimera_start={chimera_start}, range=0..{len(transcript)}"
# assert chimera_end in range(0, len(transcript)), f"chimera_end is not in proper range. " \
#                                                    f"chimera_end={chimera_end}, range=0..{len(transcript)}"
    return




def find_gambiae_region(len_5utr: int, len_cds: int, len_3utr: int,
                        chimera_start: int, chimera_end: int) -> str:
    chimera_start -= 1  # python is zerobase
    chimera_end -= 1    # python is zerobase

    range_list = get_region_ranges(len_5utr, len_cds, len_3utr)

    # check the belonging of the chimera to each region
    result = []
    for region in range(3):
        result.append(
            (chimera_start in range_list[region]) or
            (chimera_end in range_list[region]))

    # if chimera start is within the 5utr and chimera end is within the 3utr, than the cds have to be true as well
    result[1] = result[1] or (result[0] and result[2])
    result: Series = pd.Series(data=result, index=['utr5', 'cds', 'utr3'])
    return "+".join(result.index[result])


def get_gambiae_region_from_row(row: Series) -> Series:
    df["site"] = df.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                           "mRNA sequence", "chimera_start", "chimera_end",
                                           extra_chars=SITE_EXTRA_CHARS),
                          axis=1)

    try:
        return find_gambiae_region(transcript=str(row["seq"]),
                                   len_5utr=int(row["LEN_5UTR"]),
                                   len_cds=int(row["LEN_CDS"]),
                                   len_3utr=int(row["LEN_3UTR"]),
                                   chimera_start=int(row["chimera_start"])-1, #python is zerobase
                                   chimera_end=int(row["chimera_end"])-1)     #python is zerobase
    except Exception as e:
        return pd.Series(data=[""]*5, index=['utr5', 'cds', 'utr3', 'region', 'len_error'])


def insert_gambiae_region(df) ->DataFrame:
    logger.info(f"enter to insert_gambiae_region")
    df = pd.concat([df,
                    df.apply(func=get_wrapper(find_gambiae_region,
                                              "LEN_5UTR", "LEN_CDS", "LEN_3UTR",
                                              "chimera_start", "chimera_end"),
                             axis=1)], axis=1)
    return df


def find_gambiae_region_sequence(transcript: str, region: str, len_5utr: int, len_cds: int, len_3utr)\
        -> str:

    region_key = {'utr5': 0,
                  'cds' : 1,
                  'utr3' :2}

    range_list = get_region_ranges(len_5utr, len_cds, len_3utr)
    current_range = range_list[region_key[region]]
    return transcript[current_range.start:current_range.stop]



def get_gambiae_region_sequence_from_row(row: Series) -> str:
    try:
        return find_gambiae_region_sequence(transcript=str(row["seq"]),
                                   region=str(row["region"]),
                                   len_5utr=int(row["LEN_5UTR"]),
                                   len_cds=int(row["LEN_CDS"]),
                                   len_3utr=int(row["LEN_3UTR"]))

    except Exception as e:
        return ""


def insert_gambiae_region_sequence(df) ->DataFrame:
    logger.info(f"enter to insert_gambiae_region_sequence")
    df["region_sequence"] = df.apply(
        func=get_wrapper(find_gambiae_region_sequence,
                         "seq", "region", "LEN_5UTR", "LEN_CDS", "LEN_3UTR"),
        axis=1)
    return df

    #
    #
    #                                  get_gambiae_region_sequence_from_row, axis=1)
    # return df




@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def gambiae_run(fin: str, fout:str):
    df: DataFrame = add_gambiae_region_information(Path(fin))
    df = insert_gambiae_region(df)
    df = insert_gambiae_region_sequence(df)

    logger.info(df["region"].value_counts())
    # print(df.columns)
    df = df[['paper name', 'miRNA ID', 'mi_seq', 'Transcript ID', 'site', 'region', 'len_error', 'region_sequence']]


    #
    #
    # 'paper name', 'miRNA ID', 'mi_seq', 'chimera_start', 'chimera_end', 'Transcript ID', 'seq', 'site']]
    df.to_csv(REGION_PATH / source_file)


    #
    # insert_region_sequence()


if __name__ == '__main__':
    gambiae_run()

