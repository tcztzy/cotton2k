from datetime import date

import pytest

from cotton2k.utils import strptime


def test_strptime():
    assert strptime("01-OCT-2020") == date(2020, 10, 1)
    with pytest.raises(ValueError):
        strptime("Whatever")
