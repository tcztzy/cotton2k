"""Climate module"""
from __future__ import annotations

from dataclasses import dataclass
from locale import atof, atoi
from pathlib import Path
from re import findall
from typing import Any, Union


def read_climate_data(climate_file: Path) -> Climate:
    """Read climate data

    :param climate_file: the climate file path
    :type climate_file: Path
    """
    return Climate.from_file(climate_file)


def parse_weather(content: str) -> dict:
    """Parse weather file"""
    lines = content.splitlines()
    head_line, *daily_climate_lines = lines
    result: dict[str, Any] = dict()  # type: ignore
    n_length = len(head_line)
    if n_length >= 31:
        result["isw_rad"] = bool(atoi(head_line[31:34]))
    if n_length >= 34:
        result["isw_tmp"] = bool(atoi(head_line[34:37]))
    if n_length >= 37:
        result["isw_rain"] = bool(atoi(head_line[37:40]))
    if n_length >= 40:
        result["isw_wind"] = bool(atoi(head_line[40:43]))
    if n_length >= 43:
        result["isw_dewt"] = bool(atoi(head_line[43:46]))
    if n_length >= 61:
        result["average_wind"] = atof(head_line[61:71])
    clim = list()
    for line in daily_climate_lines:
        kwargs: dict[str, Union[bool, float]] = {  # type: ignore
            "_" + k: v for k, v in result.items() if k.startswith("isw")
        }
        (
            kwargs["_rad"],
            kwargs["_tmax"],
            kwargs["_tmin"],
            kwargs["_rain"],
            kwargs["_wind"],
            kwargs["_dewt"],
        ) = map(atof, findall(".{7}", line[21:]))
        c = DailyClimate(**kwargs)  # type: ignore
        clim.append(c)
    result["_climate"] = clim
    return result


@dataclass
class DailyClimate:  # pylint: disable=too-many-instance-attributes
    """Class represent daily climate"""

    _rad: float
    _tmax: float
    _tmin: float
    _rain: float
    _wind: float
    _dewt: float

    _isw_rad: bool
    _isw_tmp: bool
    _isw_rain: bool
    _isw_wind: bool
    _isw_dewt: bool

    @property
    def radiation(self):
        """Radiation"""
        return self._rad if not self._isw_rad else self._rad * 23.884

    @property
    def max_temperature(self):
        """Max temperature"""
        return self._tmax if self._isw_tmp else (self._tmax - 32) / 1.8

    @property
    def min_temperature(self):
        """Min temperature"""
        return self._tmin if self._isw_tmp else (self._tmin - 32) / 1.8

    @property
    def rain(self):
        """Rainfall"""
        return self._rain if self._isw_rain else self._rain * 25.4

    @property
    def wind(self):
        """Wind speed"""
        return self._wind if self._isw_wind else self._wind * 1.609

    @property
    def dew_temperature(self):
        """Dew point temperature"""
        return self._dewt if self._isw_dewt else (self._dewt - 32) / 1.8


@dataclass
class Climate:
    """Climate class"""

    _climate: list[DailyClimate]  # type: ignore
    isw_rad: bool
    isw_tmp: bool
    isw_rain: bool
    isw_wind: bool
    isw_dewt: bool
    average_wind: float

    @classmethod
    def from_file(cls, path: Path) -> Climate:
        """Create new Climate from file"""
        return cls(**parse_weather(path.read_text()))

    def __getitem__(
        self, item: Union[int, slice]
    ) -> Union[DailyClimate, list[DailyClimate]]:  # type: ignore
        return self._climate[item]
