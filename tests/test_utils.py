from datetime import date

from _cotton2k.utils import date2doy, doy2date


def test_date2doy():
    assert date2doy(date(2020, 10, 1)) == 275
    assert date2doy("2020-10-1") == 275
    assert date2doy(275) == 275
    assert date2doy(0) == 0


def test_doy2date():
    assert doy2date(2020, 275) == date(2020, 10, 1)
    assert doy2date(2020, 0) is None
