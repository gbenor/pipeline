from pathlib import Path
from typing import Dict, Tuple

import click
from pandas import DataFrame, Series
import pandas as pd

from consts.global_consts import DUPLEX_DICT
from duplex.Duplex import Duplex
from duplex.ViennaDuplex import ViennaDuplex
from duplex.MirandaDuplex import MirandaDuplex
from duplex.rnaHybrid import rnaHybrid

from utils.logger import logger
from utils.utils import get_wrapper, read_csv, to_csv




def do_duplex(mirna: str, target: str, cls: Duplex) -> Series:
    # print(f"mirna: {mirna}")
    # print(target)

    if pd.isna(mirna) or pd.isna(target):
        return Series({"duplex_valid" : False,
              "mrna_bulge": "",
              "mrna_inter": "",
              "mir_inter": "",
              "mir_bulge": ""})
    dp = cls.fromChimera(mirna, target)
    return Series({"duplex_valid": dp.valid,
              "mrna_bulge": dp.mrna_bulge,
              "mrna_inter": dp.mrna_inter,
              "mir_inter": dp.mir_inter,
              "mir_bulge": dp.mir_bulge})


@click.command()
@click.argument('method', type=str)
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def duplex(method: str, fin: str, fout: str):
    duplex_cls: Duplex = DUPLEX_DICT[method]
    logger.info(f"{method} do_duplex to {fin}")
    in_df: DataFrame = read_csv(Path(fin))
    # [in_df["miRNA sequence"].notnull() & in_df.site.notnull()]
    duplex_df = in_df.query("valid_row").apply(func=get_wrapper(
        do_duplex, "miRNA sequence", "site", cls=duplex_cls),
        axis=1)


    result = pd.merge(left=in_df, right=duplex_df, left_index=True, right_index=True, how='left')

    result["duplex_method"] = method
    to_csv(result, Path(fout))


@click.group()
def cli():
    pass


cli.add_command(duplex)



if __name__ == '__main__':
    cli()
