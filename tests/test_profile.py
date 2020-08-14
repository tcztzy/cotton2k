import datetime
from os import unlink
from pathlib import Path
from tempfile import mkstemp

import pytest

from cotton2k.io import ROOT_DIR, read_profile_file
from cotton2k.profile import (
    Profile,
    parse_profile,
    parse_profile_carbon_dioxide,
    parse_profile_description,
    parse_profile_field,
    parse_profile_geometry,
    parse_profile_output_flags,
    parse_profile_output_options,
    parse_profile_parameter_files,
    parse_profile_simulation_dates,
    parse_profile_soil_mulch,
    parse_profile_weather,
)

CONTENT = """test.pro            Test profile
01-MAY-2020    20-APR-2020    15-OCT-2020                        1.000  100  101
test.act                                         1     0.000     0.000  122    0
test.hyd            test.int            test.agi
    40.548    81.296  1013.000         0
    75.000     0.000    40.000         0
        10    10-APR-2020    20-OCT-2020        10    01-JUN-2020    20-OCT-2020
  0  0  1  1  0  1  0  1  1  1  1  1  0  0  0  0  0  1  0  0  0  0  0"""


def test_parse_profile():
    result = parse_profile(CONTENT)
    assert result["dayEndMulch"] == 289


def test_description():
    result = parse_profile_description("test.pro            Test profile")
    assert result["description"] == "Test profile"
    assert result["profile_file_name"] == "test.pro"
    with pytest.warns(Warning, match="file extension is not 'pro'"):
        parse_profile_description("test.json           Test profile")


def test_simulation_dates():
    line = "08-APR-1984    01-APR-1984    28-SEP-1984"
    result = parse_profile_simulation_dates(line)
    assert result.get("dateEmerge") == datetime.date(1984, 4, 8)
    assert result.get("dateSimStart") == datetime.date(1984, 4, 1)
    assert result.get("dateSimEnd") == datetime.date(1984, 9, 28)
    assert result.get("datePlant") is None
    line = "               01-APR-1984    28-SEP-1984"
    with pytest.raises(TypeError):
        parse_profile_simulation_dates(line)
    with pytest.raises(ValueError):
        parse_profile_simulation_dates('01-JAN-2020    01-FEB-2020    01-JAN-2020')


def test_carbon_dioxide():
    result = parse_profile_carbon_dioxide("")
    assert result.get("CO2EnrichmentFactor") == 0
    result = parse_profile_carbon_dioxide("     1.000  100  101")
    assert result.get("CO2EnrichmentFactor") == 1.000
    assert result.get("DayStartCO2") == 100
    assert result.get("DayEndCO2") == 101
    with pytest.raises(ValueError):
        parse_profile_carbon_dioxide("     1.000  101  100")


def test_weather():
    result = parse_profile_weather("REHA84.ACT")
    assert result.get("actualWeatherFileName") == "REHA84.ACT"


def test_soil_mulch():
    result = parse_profile_soil_mulch("")
    assert result["mulchIndicator"] == 0


def test_parameter_files():
    result = parse_profile_parameter_files(
        "HAZOR.HYD           HAZOR.INT           HAKB1.AGI"
    )
    assert result.get("soilHydraulicFileName") == "HAZOR.HYD"
    assert result.get("soilInitFileName") == "HAZOR.INT"
    assert result.get("agriculturalInputFileName") == "HAKB1.AGI"
    assert result.get("plantmapFileName") == ""


def test_geometry():
    result = parse_profile_geometry("    32.000    35.000    50.000         1")
    assert result.get("latitude") == 32.000
    assert result.get("longitude") == 35.000
    assert result.get("elevation") == 50.000
    assert result.get("siteNumber") == 1


def test_field():
    result = parse_profile_field("    96.520     0.000    10.000         4")
    assert result["rowSpace"] == 96.520
    assert result["skipRowWidth"] == 0.000
    assert result["plantsPerMeter"] == 10.000
    assert result["varNumber"] == 4


def test_output_options():
    result = parse_profile_output_options(
        "        10    20-APR-1984    20-SEP-1984        10    01-JUN-1984    20-SEP-1984"
    )
    assert result["soilMapFrequency"] == 10
    assert result["soilMapStartDate"], datetime.date(1984, 4 == 20)
    assert result["soilMapEndDate"], datetime.date(1984, 9 == 20)
    assert result["plantMapFrequency"] == 10
    assert result["plantMapStartDate"], datetime.date(1984, 6 == 1)
    assert result["plantMapEndDate"], datetime.date(1984, 9 == 20)


def test_output_flags():
    flags = "  0  0  1  1  0  1  1  1  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0"
    result = parse_profile_output_flags(flags)
    assert not result["UnitedStatesCustomarySystemOfUnitsOrInternationalSystemOfUnits"]
    assert not result["perSquareMeterOrPerPlant"]
    assert result["outputDryWeight"]


@pytest.fixture
def pro_file():
    fd, path = mkstemp(suffix=".pro", dir=ROOT_DIR / "profiles", text=True)
    with open(path, "w") as fp:
        fp.write(CONTENT)
    return Path(path)


@pytest.fixture
def tmp_file():
    fd, path = mkstemp()
    return Path(path)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_from_pro(pro_file, tmp_file):
    profile = Profile.from_pro(pro_file)
    assert profile.description == "Test profile"
    with pytest.raises(TypeError):
        Profile.from_pro(tmp_file)
    unlink(pro_file)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_read_profile_file(pro_file):
    result = read_profile_file(pro_file.name)
    assert result.description == "Test profile"
    unlink(pro_file)
    pro_file_name = 'does not exist.pro'
    with pytest.raises(FileNotFoundError):
        read_profile_file(pro_file_name)
