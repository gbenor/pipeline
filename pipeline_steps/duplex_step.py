from pathlib import Path
from typing import Tuple

import click
from pandas import DataFrame
import pandas as pd

from consts.global_consts import DUPLEX_DICT
from duplex.Duplex import Duplex
from duplex.ViennaDuplex import ViennaDuplex
from duplex.MirandaDuplex import MirandaDuplex
from duplex.rnaHybrid import rnaHybrid

from utils.logger import logger
from utils.utils import get_wrapper, read_csv, to_csv




def do_duplex(mirna: str, target: str, cls: Duplex) -> Tuple[bool, str, str, str, str]:
    # print(f"mirna: {mirna}")
    # print(target)
    if pd.isna(mirna) or pd.isna(target):
        return False, "", "", "", ""
    dp = cls.fromChimera(mirna, target)
    return (True, *dp.serialize())


@click.command()
@click.argument('method', type=str)
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def duplex(method: str, fin: str, fout: str):
    duplex_cls: Duplex = DUPLEX_DICT[method]
    logger.info(f"{method} do_duplex to {fin}")
    in_df: DataFrame = read_csv(Path(fin))
    # [in_df["miRNA sequence"].notnull() & in_df.site.notnull()]
    duplex_df = in_df.apply(func=get_wrapper(
        do_duplex, "miRNA sequence", "site", cls=duplex_cls),
        axis=1)
    result = pd.concat([in_df,
                        pd.DataFrame(duplex_df.tolist(),
                                     columns=["duplex_valid", "mrna_bulge", "mrna_inter", "mir_inter", "mir_bulge"])],
                       axis=1)
    result["duplex_method"] = method
    to_csv(result, Path(fout))


@click.group()
def cli():
    pass


cli.add_command(duplex)



if __name__ == '__main__':
    cli()
