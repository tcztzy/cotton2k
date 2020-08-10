from locale import atof, atoi
from datetime import datetime
from typing import Optional


def strptime(t) -> Optional[datetime.date]:
    t = t.strip()
    if len(t) == 11:
        return datetime.strptime(t, "%d-%b-%Y").date()


def parse_profile(content):
    lines = content.splitlines()
    result = {}
    result.update(parse_profile_description(lines[0]))
    line2 = lines[1].rstrip()
    result.update(parse_profile_simulation_dates(line2[:60]))
    if len(line2) > 76:
        result.update(parse_profile_carbon_dioxide(line2[60:]))
    else:
        result['CO2EnrichmentFactor'] = 0
    line3 = lines[2].rstrip()
    result.update(parse_profile_weather(line3[:40]))
    if len(line3) > 41:
        result.update(parse_profile_soil_mulch(line3[40:]))
    line4 = lines[3]
    result['soilHydraulicFileName'] = line4[:20].strip()
    result['soilInitFileName'] = line4[20:40].strip()
    result['agriculturalInputFileName'] = line4[40:60].strip()
    result['plantmapFileName'] = line4[60:].strip()
    result.update(parse_profile_geometry(lines[4]))
    result.update(parse_profile_field(lines[5]))
    line7 = lines[6]
    result['soilMapFrequency'] = atoi(line7[:10])
    result['soilMapStartDate'] = strptime(line7[14:25])
    result['soilMapEndDate'] = strptime(line7[29:40])
    result['plantMapFrequency'] = atoi(line7[40:50])
    result['plantMapStartDate'] = strptime(line7[54:65])
    result['plantMapEndDate'] = strptime(line7[69:80])
    result.update(parse_profile_output_flags(lines[7]))
    return result


def parse_profile_description(line):
    """Read file description"""
    return dict(
        profile_file_name=line[:20].strip(),
        description=line[20:].strip()
    )


def parse_profile_simulation_dates(line):
    """Read dates of emergence, start and end of simulation, and planting date."""
    line = line.rstrip()
    start = strptime(line[15:26])
    end = strptime(line[30:41])
    if start >= end:
        raise ValueError('Start day should be greater than end day')
    dateEmerge = strptime(line[:11])
    datePlant = strptime(line[45:56])
    if dateEmerge is None and datePlant is None:
        raise TypeError('Planting date or emergence date must be given in the profile!')
    return dict(
        dateEmerge=dateEmerge,
        dateSimStart=start,
        dateSimEnd=end,
        datePlant=datePlant
    )


def parse_profile_carbon_dioxide(line):
    """For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates for start and stop of enrichment (there are left blank if there is no CO2 enrichment)."""
    if not line:
        return dict(CO2EnrichmentFactor=0)
    start = atoi(line[10:15])
    end = atoi(line[15:])
    if start >= end:
        raise ValueError('Start day should be greater than end day')
    return dict(
        CO2EnrichmentFactor=atof(line[:10]),
        DayStartCO2=start,
        DayEndCO2=end,
    )


def parse_profile_weather(line):
    result = dict(
        actualWeatherFileName=line[:20].strip(),
        predictedWeatherFileName=line[20:40].strip(),
    )
    return result


def parse_profile_soil_mulch(line):
    pass


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
        varNumber=atoi(line[30:])
    )


def parse_profile_output_flags(line):
    result = {}
    result['UnitedStatesCustomarySystemOfUnitsOrInternationalSystemOfUnits'], result['perSquareMeterOrPerPlant'], result['outputDryWeight'], *rest = map(lambda x: bool(int(x)), line.split())
    return result


def parse_varlist(content):
    result = {}
    for var in content.splitlines():
        k = int(var[:4].strip())
        v = (var[5:25].strip(), var[40:20].strip())
        result[k] = v
