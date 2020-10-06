"""Site-specific data processing"""
from __future__ import annotations

from math import pi, sin
from operator import itemgetter
from pathlib import Path
from typing import List

from cotton2k.io import parse_parameter


def parameter(item):
    """Property for get data via human readable form"""
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
    """Site specific parameters"""

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
        """Property paramters for literally easy-understanding"""
        return self

    def deep_soil_temperature(self, daynum):
        """
        The temperature of the last soil layer (lower boundary) is computed as
        a sinusoidal function of the date (Day of year), with site-specific
        parameters.
        """
        return self[9] + self[10] * sin(2 * pi * (daynum - self[11]) / 365)

    @classmethod
    def from_dat(cls, site_dat_path: Path) -> Site:
        """Read site parameters from data file"""
        site_par = parse_parameter(site_dat_path.read_text(), 16)
        return cls(site_par)

    def __init__(self, _site_par):
        for i, value in enumerate(_site_par):
            lower, upper = PARAMETER_RANGES[i]
            if value > upper or value < lower:
                raise ValueError(
                    f"Site[{i+1}] is {value}, not in range ({lower}, {upper})"
                )
        self._site_par = _site_par

    def __getitem__(self, item):
        if item <= 0:
            raise KeyError
        return self._site_par[item - 1]

    def __setitem__(self, key, value):
        i = key - 1
        lower, upper = PARAMETER_RANGES[i]
        if value > upper or value < lower:
            raise ValueError
        self._site_par[i] = value
