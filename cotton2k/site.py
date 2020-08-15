from dataclasses import dataclass
from pathlib import Path
from cotton2k.io import parse_parameter
from typing import List


@dataclass
class Site:
    _site_par: List[float]
    wind_start_hours_after_sunraise: float

    @classmethod
    def from_dat(cls, site_dat_path: Path):
        site_par = parse_parameter(site_dat_path.read_text(), 1)
        return cls(_site_par=site_par,wind_start_hours_after_sunraise=site_par[0])

    def __getitem__(self, item):
        return self._site_par[item-1]
