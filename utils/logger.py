import sys
from logbook import Logger, RotatingFileHandler, StreamHandler

from consts.global_consts import log_file


class PipelineLogger(Logger):
    def __init__(self, log_file):
        super().__init__()
        self.handlers.append(RotatingFileHandler(log_file, bubble=True))
        self.handlers.append(StreamHandler(sys.stdout))

logger = PipelineLogger(log_file)
