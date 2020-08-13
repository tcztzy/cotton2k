from dataclasses import dataclass
from json import JSONDecoder
from locale import atof, atoi
from os.path import splitext
from pathlib import Path
from warnings import warn

from cotton2k.utils import date_to_day_of_year, strptime


def parse_profile(content):
    lines = content.splitlines()
    result = {}
    result.update(parse_profile_description(lines[0]))
    line2 = lines[1].rstrip()
    result.update(parse_profile_simulation_dates(line2[:60]))
    if len(line2) > 76:
        result.update(parse_profile_carbon_dioxide(line2[60:]))
    line3 = lines[2].rstrip()
    result.update(parse_profile_weather(line3[:40]))
    if len(line3) > 41:
        result.update(parse_profile_soil_mulch(line3[40:]))
        if result["dayEndMulch"] <= 0:
            result["dayEndMulch"] = date_to_day_of_year(result["dateSimEnd"])
    result.update(parse_profile_parameter_files(lines[3]))
    result.update(parse_profile_geometry(lines[4]))
    result.update(parse_profile_field(lines[5]))
    result.update(parse_profile_output_options(lines[6]))
    result.update(parse_profile_output_flags(lines[7]))
    return result


def parse_profile_description(line):
    """Read file description"""
    profile_file_name = line[:20].strip()
    _, ext = splitext(profile_file_name)
    if ext.lower() != ".pro":
        warn("file extension is not 'pro'")
    return dict(profile_file_name=profile_file_name, description=line[20:].strip())


def parse_profile_simulation_dates(line):
    """Read dates of emergence, start and end of simulation, and planting date."""
    line = line.rstrip()
    start = strptime(line[15:26])
    end = strptime(line[30:41])
    if start >= end:
        raise ValueError("Start day should be greater than end day")
    dateEmerge = strptime(line[:11])
    datePlant = strptime(line[45:56])
    if dateEmerge is None and datePlant is None:
        raise TypeError("Planting date or emergence date must be given in the profile!")
    return dict(
        dateEmerge=dateEmerge, dateSimStart=start, dateSimEnd=end, datePlant=datePlant
    )


def parse_profile_carbon_dioxide(line):
    """For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates for start and stop of enrichment (there are left blank if there is no CO2 enrichment)."""
    if not line:
        return dict(CO2EnrichmentFactor=0)
    start = atoi(line[10:15])
    end = atoi(line[15:])
    if start >= end:
        raise ValueError("Start day should be greater than end day")
    return dict(CO2EnrichmentFactor=atof(line[:10]), DayStartCO2=start, DayEndCO2=end,)


def parse_profile_weather(line):
    result = dict(
        actualWeatherFileName=line[:20].strip(),
        predictedWeatherFileName=line[20:40].strip(),
    )
    return result


def parse_profile_soil_mulch(line):
    MulchIndicator = atoi(line[:10]) if line else 0
    result = {"mulchIndicator": MulchIndicator}
    if MulchIndicator > 0:
        result["mulchTranSW"] = atof(line[10:20])
        result["mulchTranLW"] = atof(line[20:30])
        result["dayStartMulch"] = atoi(line[30:35])
        result["dayEndMulch"] = atoi(line[35:40])
    return result


def parse_profile_parameter_files(line):
    """Names of files of soil hydraulic data, soil initial conditions, agricultural input, and plant map adjustment."""
    return dict(
        soilHydraulicFileName=line[:20].strip(),
        soilInitFileName=line[20:40].strip(),
        agriculturalInputFileName=line[40:60].strip(),
        plantmapFileName=line[60:].strip(),
    )


def parse_profile_geometry(line):
    """Latitude and longitude of this site, elevation (in m above sea level), and the number for this geographic site."""
    return dict(
        latitude=atof(line[:10]),
        longitude=atof(line[10:20]),
        elevation=atof(line[20:30]),
        siteNumber=atoi(line[30:]),
    )


def parse_profile_field(line):
    """Row spacing in cm, skip-row spacing in cm (blank or 0 for no skip rows), number of plants per meter of row, and index number of the cultivar."""
    return dict(
        rowSpace=atof(line[:10]),
        skipRowWidth=atof(line[10:20]),
        plantsPerMeter=atof(line[20:30]),
        varNumber=atoi(line[30:]),
    )


def parse_profile_output_options(line):
    """Frequency in days for output of soil maps, and dates for start and stop of this output (blank or 0 if no such output is required. Same is repeated for output of plant maps."""
    return dict(
        soilMapFrequency=atoi(line[:10]),
        soilMapStartDate=strptime(line[14:25]),
        soilMapEndDate=strptime(line[29:40]),
        plantMapFrequency=atoi(line[40:50]),
        plantMapStartDate=strptime(line[54:65]),
        plantMapEndDate=strptime(line[69:80]),
    )


def parse_profile_output_flags(line):
    result = {}
    (
        result["UnitedStatesCustomarySystemOfUnitsOrInternationalSystemOfUnits"],
        result["perSquareMeterOrPerPlant"],
        result["outputDryWeight"],
        *rest,
    ) = map(lambda x: bool(int(x)), line.split())
    return result


class ProfileJSONDecoder(JSONDecoder):
    pass


@dataclass
class Profile:
    profile_file_path: Path
    description: str

    @classmethod
    def from_pro(cls, profile_file_path: Path) -> None:
        warn(
            "Old-style profile file format is deprecated, "
            "and will be remove in stable release.",
            DeprecationWarning,
        )
        _, ext = splitext(profile_file_path.name)
        if ext.lower() != ".pro":
            raise TypeError(f"Expect 'pro' file, got {profile_file_path.name}")
        return cls(**parse_profile(profile_file_path.read_text()))
