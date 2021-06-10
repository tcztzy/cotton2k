"""Cotton2k model."""
from importlib.metadata import metadata, version
from pathlib import Path
from typing import Union

from .io import read_input, write_output

__all__ = ("run",)

__version__ = version(__name__)
meta = metadata(__name__)
__author__: str = meta["Author"]
__license__: str = meta["License"]


def run(
    profile_path: Union[Path, str, dict], *, output: bool = False, name: str = "default"
):
    sim = read_input(profile_path)
    sim.run()
    if output:
        if isinstance(profile_path, dict):
            output_dir = Path(".")
        else:
            output_dir = Path(profile_path).parent
            name = Path(profile_path).stem
        write_output(sim, output_dir, name)
    return sim
