"""Site-specific data processing"""
from __future__ import annotations

from functools import wraps
from operator import itemgetter
from pathlib import Path
from typing import List, Optional, Tuple

from cotton2k.io import parse_parameter


def parameter(item):
    return property(itemgetter(item))


PARAMETER_RANGES = (
    (0, 2),
    (2, 4),
    (0, 2.5),
    (0.0025, 0.0080),
    (6.264, 20),
    (14.554, 24.302),
    (40, 120),
    (1, 3.5),
    (20, 24),
    (3, 6.5),
    (180, 180),
)


class Site:
    _site_par: List[float]
    wind_start_hours_after_sunrise: float = parameter(1)
    wind_max_hours_after_noon: float = parameter(2)
    wind_stop_hours_after_sunset: float = parameter(3)
    night_wind_factor: float = parameter(4)
    dew_temperature_when_max_temperture_over_40: float = parameter(5)
    dew_temperature_when_max_temperture_below_20: float = parameter(6)
    cloud_correction: float = parameter(7)
    hours_after_noon_reatch_max_temperature: float = parameter(8)
    deep_soil_temperature_offset: float = parameter(9)
    deep_soil_temperature_gain: float = parameter(10)
    deep_soil_temperature_shift: float = parameter(11)

    @property
    def parameters(self):
        class _:
            def __getitem__(myself, item):
                if item <= 0:
                    raise KeyError
                return self._site_par[item - 1]

            def __setitem__(myself, key, value):
                i = key - 1
                lower, upper = PARAMETER_RANGES[i]
                if value > upper or value < lower:
                    raise ValueError
                self._site_par[i] = value

        return _()

    @parameters.setter
    def parameters(self, value):
        for i, v in enumerate(value):
            lower, upper = PARAMETER_RANGES[i]
            if v > upper or v < lower:
                raise ValueError(f"Site[{i+1}] is {v}, not in range ({lower}, {upper})")
        self._site_par = value

    @classmethod
    def from_dat(cls, site_dat_path: Path) -> Site:
        site_par = parse_parameter(site_dat_path.read_text(), 16)
        return cls(site_par)

    def __init__(self, _site_par):
        self.parameters = _site_par

    def __getitem__(self, item):
        return self.parameters[item]

    def __setitem__(self, key, value):
        self.parameters[key] = value
