from datetime import date

from _cotton2k.utils import date2doy


def test_date_to_day_of_year():
    assert date2doy(date(2020, 10, 1)) == 275
    assert date2doy("2020-10-1") == 275
    assert date2doy(275) == 275
    assert date2doy(0) == 0
