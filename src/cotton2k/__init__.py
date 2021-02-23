"""Cotton2k model."""
from importlib.metadata import metadata, version

from _cotton2k import run  # pylint: disable=import-error# noqa: F401

__all__ = ("run",)

__version__ = version(__name__)
meta = metadata(__name__)
__author__: str = meta["Author"]
__license__: str = meta["License"]
