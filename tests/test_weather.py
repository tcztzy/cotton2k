from os import linesep
from pathlib import Path

import pytest

from cotton2k.weather import read_climate_data


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
    assert result[0].radiation == 20 * 23.884
    assert result[0].max_temperature == 13.4
    assert result[0].min_temperature == 4
    assert result[0].rain == 0
    assert result[0].wind == 110
    assert result[0].dew_temperature == 10
