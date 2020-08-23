from dataclasses import dataclass
from typing import List
from pathlib import Path
from re import findall
from locale import atof


def read_climate_data():
    pass


def parse_weather(content: str):
    lines = content.splitlines()
    clim = map(lambda line: Climate(*map(atof, findall(".{7}", line[21:]))), lines[1:])
    return dict(clim=clim)


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

    @classmethod
    def from_file(cls, path: Path):
        return cls(**parse_weather(path.read_text()))
