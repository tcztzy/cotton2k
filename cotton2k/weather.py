from __future__ import annotations

from dataclasses import dataclass
from locale import atof, atoi
from pathlib import Path
from re import findall
from typing import Any, Dict, List


def read_climate_data(weather_file):
    return Weather.from_file(weather_file)


def parse_weather(content: str):
    lines = content.splitlines()
    head_line, *daily_climate_lines = lines
    result: Dict[str, Any] = dict()
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
        radiation, *rest = map(atof, findall(".{7}", line[21:]))
        clim.append(Climate(radiation, *rest))
    result["clim"] = clim
    return result


@dataclass
class Climate:
    radiation: float
    max_temperature: float
    min_temperature: float
    rain: float
    wind: float
    dew_temperature: float


@dataclass
class Weather:
    clim: List
    isw_rad: bool
    isw_tmp: bool
    isw_rain: bool
    isw_wind: bool
    isw_dewt: bool
    average_wind: float

    @classmethod
    def from_file(cls, path: Path):
        return cls(**parse_weather(path.read_text()))

    def __getitem__(self, item):
        return self.clim[item]
