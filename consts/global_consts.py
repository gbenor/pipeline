from pathlib import Path
from duplex.ViennaDuplex import ViennaDuplex
from duplex.MirandaDuplex import MirandaDuplex
from duplex.rnaHybrid import rnaHybrid

ROOT_PATH = Path ("/home/local/BGU-USERS/benorgi/pipeline/")
log_file = ROOT_PATH / "pipeline_log.log"

SITE_EXTRA_CHARS: int = 3
HUMAN_SITE_EXTENDED_LEN = 25

MINIMAL_BLAST_COVERAGE = 95
MINIMAL_BLAST_IDENTITY = 95

MINIMAL_LENGTH_TO_BLAST = 10

DUPLEX_DICT = {"ViennaDuplex": ViennaDuplex,
               "MirandaDuplex" : MirandaDuplex,
               "rnaHybrid": rnaHybrid}

