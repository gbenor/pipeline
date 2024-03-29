from consts.global_consts import ROOT_PATH
from numpy import int64

READ_PATH = ROOT_PATH / "data/pipeline_steps/read"
MIRNA_SEQ_PATH = ROOT_PATH / "data/pipeline_steps/mirna_sequence"
SITE_PATH = ROOT_PATH / "data/pipeline_steps/site"
REGION_PATH = ROOT_PATH / "data/pipeline_steps/region"
CONCAT_BLAST = ROOT_PATH / "data/pipeline_steps/concat_blast"
NORMALIZATION_PATH = ROOT_PATH / "data/pipeline_steps/normalization_final"

NORMALIZATION_COLUMNS = ['key', 'paper name', 'organism', 'miRNA ID', 'miRNA sequence', 'seed_family',
                         'site',
                         'region', 'region count', 'paper region',
                         'start', 'end', 'sequence', 'Gene_ID',
                         'identity', 'coverage',
                         'valid_row']

GAMBIAE_INFORMATION_FILENAME = ROOT_PATH / "data/gambiae/AgamTransP44updateGF.txt"
GAMBIAE_INFORMATION_COLUMN_NAMES = ["TRANSCRIPT_FAMILY", "TRANSCRIPT_ID", "SEQ_TYPE",
                                    "LEN_5UTR", "LEN_CDS", "LEN_3UTR", "COMMENTS"]

GAMBIAE_INFORMATION_USECOLS = ["TRANSCRIPT_ID", "LEN_5UTR", "LEN_CDS", "LEN_3UTR"]
GAMBIAE_INFORMATION_DTYPE = {'TRANSCRIPT_ID': str,
                             'LEN_5UTR': int64,
                             'LEN_CDS': int64,
                             'LEN_3UTR': int64}


