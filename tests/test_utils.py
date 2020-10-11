from datetime import date

import pytest

from cotton2k.utils import date_to_day_of_year, strptime


def test_strptime():
    assert strptime("01-OCT-2020") == date(2020, 10, 1)
    with pytest.raises(ValueError):
        strptime("Whatever")


def test_date_to_day_of_year():
    assert date_to_day_of_year(date(2020, 10, 1)) == 275
    assert date_to_day_of_year(date(2020, 1, 1), 2019) == 366
