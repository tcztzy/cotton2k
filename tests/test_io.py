import json

import pytest

from _cotton2k import Simulation
from cotton2k.io import read_input, write_output


def test_read_input(empty_json, test_json):
    with pytest.raises((KeyError, TypeError)):
        read_input(empty_json)
    with pytest.raises((KeyError, TypeError)):
        read_input({})
    with pytest.raises((KeyError, TypeError)):
        read_input(str(empty_json))
    read_input(test_json)


def test_write_output(sim: Simulation, test_json):
    write_output(sim, test_json.parent, "alaer")
    alaer = test_json.parent / "alaer"
    assert alaer.exists()
    files = list(alaer.glob("**/1984-09-28/index.json"))
    assert len(files) == 1
    with open(files[0]) as fp:
        result = json.load(fp)
        assert result["lint_yield"] == 2205.223254512975
