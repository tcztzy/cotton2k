from dataclasses import dataclass
from pathlib import Path
from cotton2k.io import parse_parameter
from typing import List


@dataclass
class Site:
    _site_par: List[float]
    wind_start_hours_after_sunraise: float
    wind_max_hours_after_noon: float
    wind_stop_hours_after_sunset: float

    @classmethod
    def from_dat(cls, site_dat_path: Path):
        site_par = parse_parameter(site_dat_path.read_text(), 16)
        return cls(
            _site_par=site_par,
            wind_start_hours_after_sunraise=site_par[0],
            wind_max_hours_after_noon=site_par[1],
            wind_stop_hours_after_sunset=site_par[2],
        )

    def __getitem__(self, item):
        return self._site_par[item - 1]
