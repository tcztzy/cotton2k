from cotton2k.climate import (
    parse_weather,
    read_climate_data,
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


def test_parse_weather():
    parse_weather("Whatever")  # for partial coverage
