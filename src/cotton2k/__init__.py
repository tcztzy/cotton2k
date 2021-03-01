"""Cotton2k model."""
import json
from importlib.metadata import metadata, version
from pathlib import Path

from _cotton2k import _Simulation  # pylint: disable=import-error# noqa: F401

__all__ = ("run",)

__version__ = version(__name__)
meta = metadata(__name__)
__author__: str = meta["Author"]
__license__: str = meta["License"]


def run(profile_path: Path):
    sim = _Simulation()
    sim.read_input(**json.loads(profile_path.read_text()))
    sim.run()
    return sim
