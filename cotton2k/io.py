from locale import atof, atoi
from pathlib import Path
from typing import Sequence
from os.path import splitext
import warnings

from appdirs import user_data_dir

from cotton2k.utils import date_to_day_of_year, strptime

ROOT_DIR = Path(user_data_dir("cotton2k", "Tang Ziya"))


def read_profile_file(profile_file_name):
    """
    TODO: Maybe profile should be a object instead of dict
    TODO: Profile is self defined file format and not readable for human,
    maybe JSON and TOML are better alternatives"""
    path = ROOT_DIR / "profiles" / profile_file_name
    if not path.exists():
        raise FileNotFoundError(f"{path} not found!")
    return parse_profile(path.read_text())


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
        warnings.warn("file extension is not 'pro'")
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


def read_calibration_data(var_number: int, site_number: int):
    """
    This function reads the values of the calibration parameters
    from input files. It is called from ReadInput(). It calls GetLineData().

    The following global or file-scope variables are set here:
    SiteName, SitePar, VarName, VarPar
    
    TODO: Maybe JSON or CSV is more suitable file format than self defined DAT files.
    """
    data_dir = ROOT_DIR / "data"
    vars_dir = data_dir / "vars"
    varlist = vars_dir / "varlist.dat"
    var_name, var_file = parse_list_dat(varlist.read_text())[var_number]
    if not (var_file_path := vars_dir / var_file).exist():
        raise FileNotFoundError(f"{var_file_path} not found!")
    var_par = parse_parameter(var_file_path.read_text(), 60)
    site_dir = data_dir / "site"
    sitelist = site_dir / "sitelist.dat"
    site_name, site_file = parse_list_dat(sitelist.read_text())[site_number]
    if not (site_file_path := site_dir / site_file).exist():
        raise FileNotFoundError(f"{site_file_path} not found!")
    site_par = parse_parameter(site_file_path.read_text(), 20)
    return dict(siteName=site_name, sitePar=site_par, varName=var_name, varPar=var_par,)


def parse_list_dat(content: str) -> dict:
    result = {}
    for var in content.splitlines():
        k = int(var[:4].strip())
        v = (var[5:25].strip(), var[40:60].strip())
        result[k] = v
    return result


def parse_parameter(content: str, number: int) -> Sequence:
    lines = content.splitlines()
    return list(map(lambda line: atof(line[:15]), lines[1 : number + 1]))
