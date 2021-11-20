from pathlib import Path

import click
from pandas import DataFrame
import pandas as pd

from consts.mirna_utils import MIRBASE_FILE
from utils.utils import read_csv

SPLUNK_COL_TO_DROP = [
    "paper name",
    "organism",
    "sequence",
    "duplex_method"
]

@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def splunk(fin: str, fout: str):
    in_df: DataFrame = read_csv(Path(fin))
    in_df = in_df[in_df['valid_row']]
    in_df = in_df.astype({'duplex_valid':bool})
    in_df = in_df[in_df['duplex_valid']]
    print(in_df['duplex_valid'].unique())

    in_df.drop(columns=SPLUNK_COL_TO_DROP, inplace=True)
    Path(fout).parent.mkdir(parents=True, exist_ok=True)
    in_df.to_csv(fout)


@click.command()
@click.argument('fin', type=str)
def mirnaid_fix(fin: str):
    d = {"mouse": "mmu",
         "human": "hsa",
         "elegans": "cel",
         "cattle": "bta",
         "fly": "aga"
         }
    prefix = None
    for k,v in d.items():
        if k in fin:
            prefix = v
    if prefix is None:
        raise Exception("unrecognized mirbase prefix")

    mirbase_df: DataFrame = pd.read_csv(MIRBASE_FILE).query("prefix==@prefix")
    mirbase_df.sort_values(by="version", ascending=False, inplace=True)
    mirbase_df.drop_duplicates("miRNA sequence", keep="first", inplace=True)
    d = pd.read_csv(fin, index_col=0)
    join_df = d.merge(mirbase_df, how="left", left_on="miRNA sequence", right_on="miRNA sequence")
    d['miRNA ID'] = join_df['miRNA ID_y']
    d.to_csv(fin)



@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def splunk_lite(fin: str, fout: str):
    def con(x, ll):
        for l in ll:
            if x.find(l) > -1:
                return False
        return True

    features_to_drop = ["MRNA_Down", "Acc_", "HotPairing", "Energy_MEF", "MRNA_Up", "MRNA", "miRNAPairingCount",
                        "Seed_match"]

    d: DataFrame = pd.read_csv(fin)
    f = [x for x in d.columns if con(x, features_to_drop)]
    d = d[f]

    a = [x for x in d.columns if "miRNAMatchPosition" in x]
    d[a] = d[a].replace(["GU", "AU", "GC"], 100)
    d[a] = d[a].replace(["MM", "BB"], 0)
    d["seed"] = d["miRNA sequence"].apply(lambda x: x[:8])

    Path(fout).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(fout)





@click.group()
def cli():
    pass




cli.add_command(splunk)
cli.add_command(mirnaid_fix)
cli.add_command(splunk_lite)


if __name__ == '__main__':
    cli()

