"""Input/Output"""
import json
from pathlib import Path

from _cotton2k import _Simulation  # pylint: disable=import-error# noqa: F401


def read_input(path: Path) -> _Simulation:
    sim = _Simulation()
    kwargs = json.loads(path.read_text())
    for attr in [
        "start_date",
        "stop_date",
        "emerge_date",
        "plant_date",
    ]:
        if attr in kwargs:
            setattr(sim, attr, kwargs.get(attr))
    sim.year = int(kwargs["start_date"][:4])
    sim.read_input(**kwargs)
    return sim
