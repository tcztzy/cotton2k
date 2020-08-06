import os
from locale import atof, atoi
from pathlib import Path
from datetime import datetime

from .utils import get_line_data


def parse_profile(content):
    lines = content.splitlines()
    line = lines[1].strip()
    # For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates 
    # for start and stop of enrichment (these are left blank if there is no CO2 enrichment).
    if len(line) > 76:
        ttt = line[60:70].strip()
        CO2EnrichmentFactor = atof(ttt)
        ttt = line[70:75].strip()
        DayStartCO2 = atoi(ttt)
        ttt = line[75:].strip()
        DayEndCO2 = atoi(ttt)
    else:
        CO2EnrichmentFactor = 0
    return {
        'description': lines[0][20:].strip(),
        'dateEmerge': datetime.strptime(lines[1][:11], "%d-%b-%Y").date(),
        'dateSimStart': datetime.strptime(lines[1][15:26], "%d-%b-%Y").date(),
        'dateSimEnd': datetime.strptime(lines[1][30:41], "%d-%b-%Y").date(),
        'CO2EnrichmentFactor': CO2EnrichmentFactor,
        'actualWeatherFileName': lines[2][:20].strip(),
    }
