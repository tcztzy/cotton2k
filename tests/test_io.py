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


def test_write_output(sim: Simulation):
    session, simulation = write_output(sim)
    assert len(simulation.states) == 181
    assert simulation.states[-1].lint_yield == 2205.223254512975
    assert simulation.states[-1].number_of_fruiting_branches == 20
