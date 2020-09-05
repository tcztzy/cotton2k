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
    if len(head_line) >= 31:
        result["isw_rad"] = bool(atoi(head_line[31:34]))
    if len(head_line) >= 34:
        result["isw_tmp"] = bool(atoi(head_line[34:37]))
    if len(head_line) >= 37:
        result["isw_rain"] = bool(atoi(head_line[37:40]))
    if len(head_line) >= 40:
        result["isw_wind"] = bool(atoi(head_line[40:43]))
    if len(head_line) >= 43:
        result["isw_dewt"] = bool(atoi(head_line[43:46]))
    if len(head_line) >= 61:
        result["average_wind"] = atof(head_line[61:71])
    result["clim"] = list(
        map(
            lambda line: Climate(*map(atof, findall(".{7}", line[21:]))),
            daily_climate_lines,
        )
    )
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
