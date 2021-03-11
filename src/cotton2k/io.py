"""Input/Output"""
import json
from pathlib import Path
from typing import Optional

from _cotton2k import _Simulation  # pylint: disable=import-error# noqa: F401


def read_input(path: Path) -> tuple[_Simulation, dict]:
    sim = _Simulation()
    kwargs = json.loads(path.read_text())
    for attr in [
        "start_date",
        "stop_date",
        "emerge_date",
        "plant_date",
        "latitude",
        "longitude",
        "elevation",
        "site_parameters",
        "cultivar_parameters",
        "output_flags",
        "row_space",
        "skip_row_width",
        "plants_per_meter",
    ]:
        if attr in kwargs:
            setattr(sim, attr, kwargs.get(attr))
    sim.year = int(kwargs["start_date"][:4])
    sim.read_input(**kwargs)
    return sim, kwargs


def write_output(sim: _Simulation, path: Optional[Path] = None) -> None:
    states = list(sim.states)
    with open(path or "UNKNOWN.cotton2k-output.json", "w") as f:
        json.dump(states, f)
