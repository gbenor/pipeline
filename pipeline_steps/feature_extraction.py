from pathlib import Path
from features.AccessibilityFeatures import AccessibilityFeatures
from features.EnergyFeatures import EnergyFeatures
from features.MatchingFeatures import MatchingFeatures
from features.MrnaFeatures import MrnaFeatures
from features.SeedFeatures import SeedFeatures
from utils.logger import logger
from utils.utils import apply_in_chunks, get_wrapper, read_csv, to_csv
from pandas import Series, DataFrame
from duplex.Duplex import Duplex
import pandas as pd
import click

FEATURE_CLASSES = [SeedFeatures, MatchingFeatures, MrnaFeatures, EnergyFeatures, AccessibilityFeatures]


def row_feature_extractor(miRNA: str, site: str, start: int, end: int, sequence: str,
                          mrna_bulge: str, mrna_inter: str, mir_inter: str, mir_bulge: str) -> Series:

    dp = Duplex.fromStrings(mrna_bulge, mrna_inter, mir_inter, mir_bulge)
    f = [feature_cls(dp, miRNA, site, start, end, sequence).get_features() for feature_cls in FEATURE_CLASSES]
    return pd.concat(f)


def df_feature_extractor(valid_df: DataFrame) -> DataFrame:
    return apply_in_chunks(df=valid_df,
                           func=get_wrapper(row_feature_extractor,
                                           "miRNA sequence", "site", "start", "end", "sequence",
                                           'mrna_bulge', 'mrna_inter', 'mir_inter', 'mir_bulge')
                           )

     # return valid_df.apply(func=get_wrapper(row_feature_extractor,
     #                                           "miRNA sequence", "site", "start", "end", "sequence",
     #                                           'mrna_bulge', 'mrna_inter', 'mir_inter', 'mir_bulge'),
     #                          axis=1)




@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def feature_extraction(fin: str, fout: str):
    in_df: DataFrame = read_csv(Path(fin))
    valid_df = in_df.query("valid_row & duplex_valid=='True'")
    feature_df = df_feature_extractor(valid_df)
    result = pd.merge(left=in_df, right=feature_df, left_index=True, right_index=True, how='left')
    to_csv(result, Path(fout))

@click.group()
def cli():
    pass


cli.add_command(feature_extraction)


if __name__ == '__main__':
    cli()

