from importlib.metadata import PackageNotFoundError, metadata, version
from pathlib import Path

from appdirs import user_data_dir

__version__ = version(__name__)
meta = metadata(__name__)
__author__ = meta["Author"]
__license__ = meta["License"]
ROOT_DIR = Path(user_data_dir(__name__, __author__, __version__, True))
