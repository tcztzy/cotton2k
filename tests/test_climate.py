from datetime import date, datetime, timedelta

from cotton2k.climate import (
    compute_day_length,
    parse_weather,
    read_climate_data,
    tdewest,
    vapor_pressure,
)


def test_read_climate_data(weather_file):
    result = read_climate_data(weather_file)
    assert result.isw_rad
    assert result.isw_tmp
    assert result.isw_rain
    assert result.isw_wind
    assert result.isw_dewt
    assert result.average_wind == 0
    climate = result[0]
    assert climate.radiation == 20 * 23.884
    assert climate.max_temperature == 13.4
    assert climate.min_temperature == 4
    assert climate.rain == 0
    assert climate.wind == 110
    assert climate.dew_temperature == 10


def test_tdewest():
    assert tdewest(16, 10, 20) == 10
    assert tdewest(40, 10, 20) == 20
    assert tdewest(30, 10, 20) == 15


def test_parse_weather():
    parse_weather("Whatever")  # for partial coverage


def test_vapor_pressure():
    def good_enough(a, b):
        return abs(1 - b / a) < 0.0009

    for t, p in ((0, 0.6113), (20, 2.3388), (35, 5.6267), (50, 12.344)):
        assert good_enough(p, vapor_pressure(t))


def test_compute_day_length():
    day_length, sunrise, solar_noon, sunset, *rest = compute_day_length(
        date(2020, 10, 20), 40.54778, 81.29
    )
    results = {
        "sunrise": "2020-10-20T00:52:07+00:00",
        "sunset": "2020-10-20T11:46:53+00:00",
        "solar_noon": "2020-10-20T06:19:30+00:00",
        "day_length": 39286,
        "civil_twilight_begin": "2020-10-20T00:24:28+00:00",
        "civil_twilight_end": "2020-10-20T12:14:31+00:00",
        "nautical_twilight_begin": "2020-10-19T23:52:41+00:00",
        "nautical_twilight_end": "2020-10-20T12:46:19+00:00",
        "astronomical_twilight_begin": "2020-10-19T23:21:05+00:00",
        "astronomical_twilight_end": "2020-10-20T13:17:55+00:00",
    }
    threshold = timedelta(minutes=3)
    assert abs(sunrise - datetime.fromisoformat(results["sunrise"])) <= threshold
    assert abs(sunset - datetime.fromisoformat(results["sunset"])) <= threshold
    assert abs(solar_noon - datetime.fromisoformat(results["solar_noon"])) <= threshold
    assert abs(day_length - timedelta(seconds=39286)) <= 2 * threshold
