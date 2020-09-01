"""Cotton2k model."""
from importlib.metadata import metadata, version

__version__ = version(__name__)
meta = metadata(__name__)
__author__: str = meta["Author"]
__license__: str = meta["License"]
