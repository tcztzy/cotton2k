"""Cotton2k model."""
from importlib.metadata import metadata, version
from pathlib import Path

from .io import read_input, write_output

__all__ = ("run",)

__version__ = version(__name__)
meta = metadata(__name__)
__author__: str = meta["Author"]
__license__: str = meta["License"]


def run(profile_path: Path):
    sim = read_input(profile_path)
    sim.run()
    write_output(sim, profile_path.parent / (profile_path.stem + "-output.json"))
    return sim
