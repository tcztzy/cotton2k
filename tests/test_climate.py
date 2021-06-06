import datetime

from _cotton2k.climate import compute_day_length


def test_compute_day_length():
    lat = 40.54778
    lon = 81.29
    tz = datetime.timezone(datetime.timedelta(hours=lon // 15))
    date = datetime.date(2020, 10, 20)
    zero_o_clock = datetime.datetime.combine(date, datetime.time(0, 0, 0), tzinfo=tz)
    result = compute_day_length((lat, lon), date)
    sunrise = zero_o_clock + datetime.timedelta(hours=result["sunr"])
    solar_noon = zero_o_clock + datetime.timedelta(hours=result["solar_noon"])
    assert (
        abs(
            solar_noon - datetime.datetime.fromisoformat("2020-10-20T06:19:30+00:00")
        ).seconds
        <= 120
    )
    assert (
        abs(
            sunrise - datetime.datetime.fromisoformat("2020-10-20T00:52:07+00:00")
        ).seconds
        <= 240
    )
