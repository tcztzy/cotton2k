import json
from os import linesep
from pathlib import Path

from pytest import fixture


@fixture
def weather_file(tmp_path: Path) -> Path:
    file_path = tmp_path / "test.act"
    lines = [
        f"{'Whatever':<30}{1:>3}{1:>3}{1:>3}{1:>3}{1:>3}{' '*15}{0:>10.2f}",
        f"{92:>4}{'01-APR-2020':^17}{20:>7.2f}{13.4:>7.2f}{4:>7.2f}{0:>7.2f}{110:>7.2f}{10:>7.2f}",
    ]
    file_path.write_bytes(linesep.join(lines).encode())
    return file_path


@fixture
def empty_json(tmp_path: Path) -> Path:
    empty = tmp_path / "empty.cotton2k.json"
    empty.write_text(json.dumps({}))
    return empty
