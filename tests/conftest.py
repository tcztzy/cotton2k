import json
from os import linesep
from pathlib import Path

from pytest import fixture


@fixture
def empty_json(tmp_path: Path) -> Path:
    empty = tmp_path / "empty.cotton2k.json"
    empty.write_text(json.dumps({}))
    return empty
