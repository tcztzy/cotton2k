from datetime import date, datetime, timedelta

from cotton2k.climate import (
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
