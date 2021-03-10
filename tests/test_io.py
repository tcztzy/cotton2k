import pytest

from cotton2k.io import read_input


def test_read_input(empty_json):
    with pytest.raises(KeyError, match="start_date"):
        read_input(empty_json)
