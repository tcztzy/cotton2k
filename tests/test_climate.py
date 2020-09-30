from os import linesep
from pathlib import Path

import pytest

from cotton2k.climate import read_climate_data, tdewest


@pytest.fixture
def weather_file(tmp_path: Path) -> Path:
    file_path = tmp_path / "test.act"
    lines = [
        f"{'Whatever':<30}{1:>3}{1:>3}{1:>3}{1:>3}{1:>3}{' '*15}{0:>10.2f}",
        f"{92:>4}{'01-APR-2020':^17}{20:>7.2f}{13.4:>7.2f}{4:>7.2f}{0:>7.2f}{110:>7.2f}{10:>7.2f}",
    ]
    file_path.write_bytes(linesep.join(lines).encode())
    return file_path


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
