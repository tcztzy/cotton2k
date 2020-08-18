from dataclasses import dataclass
from typing import List
from pathlib import Path
from re import findall


def parse_weather(content: str):
    lines = content.splitlines()
    clim = map(lambda line: findall('.{7}', line[21:]), lines[1:])
    return dict(
        clim=clim
    )

@dataclass
class Weather:
    clim: List

    @classmethod
    def from_file(cls, path: Path):
        return cls(**parse_weather(path.read_text()))
