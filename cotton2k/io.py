from locale import atof, atoi
from datetime import datetime

def strptime(time):
    return datetime.strptime(time, "%d-%b-%Y").date()


def parse_profile(content):
    lines = content.splitlines()
    result = {}
    result['description'] = lines[0][20:].strip()
    line = lines[1].strip()
    result['dateEmerge'] = strptime(line[:11])
    result['dateSimStart'] = strptime(line[15:26])
    result['dateSimEnd'] = strptime(line[30:41])
    # For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates 
    # for start and stop of enrichment (these are left blank if there is no CO2 enrichment).
    if len(line) > 76:
        ttt = line[60:70].strip()
        result['CO2EnrichmentFactor'] = atof(ttt)
        ttt = line[70:75].strip()
        result['DayStartCO2'] = atoi(ttt)
        ttt = line[75:].strip()
        result['DayEndCO2'] = atoi(ttt)
    else:
        result['CO2EnrichmentFactor'] = 0
    result['actualWeatherFileName'] = lines[2][:20].strip()
    line4 = lines[3]
    result['soilHydraulicFileName'] = line4[:20].strip()
    result['soilInitFileName'] = line4[20:40].strip()
    result['agriculturalInputFileName'] = line4[40:60].strip()
    result['plantmapFileName'] = line4[60:].strip()
    line5 = lines[4]
    result['latitude'] = atof(line5[:10])
    result['longitude'] = atof(line5[10:20])
    result['elevation'] = atof(line5[20:30])
    result['siteNumber'] = atoi(line5[30:])
    line6 = lines[5]
    result['rowSpace'] = atof(line6[:10])
    result['skipRowWidth'] = atof(line6[10:20])
    result['plantsPerMeter'] = atof(line6[20:30])
    result['varNumber'] = atoi(line6[30:])
    line7 = lines[6]
    result['soilMapFrequency'] = atoi(line7[:10])
    result['soilMapStartDate'] = strptime(line7[14:25])
    result['soilMapEndDate'] = strptime(line7[29:40])
    result['plantMapFrequency'] = atoi(line7[40:50])
    result['plantMapStartDate'] = strptime(line7[54:65])
    result['plantMapEndDate'] = strptime(line7[69:80])
    return result
