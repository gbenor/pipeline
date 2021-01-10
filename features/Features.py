from abc import ABC, abstractmethod
from functools import reduce
from typing import Dict

from pandas import Series

from duplex.Duplex import Duplex
import pandas as pd


class Features(ABC):
    def __init__(self, duplex: Duplex, miRNA_sequence: str, site: str,
                 start: int, end: int, region_sequence: str):
        self._duplex = duplex
        self._miRNA_sequence = miRNA_sequence
        self._site = site
        self._start = int(start)
        self._end = int(end)
        self._region_sequence = region_sequence

        self._features_dict: Dict = {}

    def get_features(self) -> Series:
        def u(x, y):
            x.update(y)
            return x

        self.extract_features()
        united_dicts = reduce(u, self._features_dict.values())
        return pd.Series(united_dicts)

    @abstractmethod
    def extract_features(self):
        pass

