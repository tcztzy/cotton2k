import datetime
from calendar import isleap
from math import acos, cos, floor, pi, sin, tan

from _cotton2k.utils import date2doy


def compute_day_length(coordinate: tuple[float, float], date: datetime.date):
    lat, lon = coordinate
    xday: float = 2 * pi * date2doy(date) / (365 + int(isleap(date.year)))
    declination = (
        0.006918
        - 0.399912 * cos(xday)
        + 0.070257 * sin(xday)
        - 0.006758 * cos(2 * xday)
        + 0.000907 * sin(2 * xday)
        - 0.002697 * cos(3 * xday)
        + 0.001480 * sin(3 * xday)
    )
    exday = (
        (
            0.000075
            + 0.001868 * cos(xday)
            - 0.032077 * sin(xday)
            - 0.014615 * cos(2 * xday)
            - 0.04089 * sin(2 * xday)
        )
        * 12.0
        / pi
    )
    st = 15 * floor(lon / 15)
    f = (lon - st) / 15
    if lon > 0:
        if f > 0.5:
            f -= 1
    elif f < -0.5:
        f += 1

    solar_noon = 12 - f - exday
    xlat = lat * pi / 180
    ht = -tan(xlat) * tan(declination)
    if ht > 1:
        ht = 1
    elif ht < -1:
        ht = -1

    day_length = 2 * acos(ht) * 12 / pi
    sunr = solar_noon - day_length / 2
    return {
        "day_length": day_length,
        "solar_noon": solar_noon,
        "sunr": sunr,
        "suns": sunr + day_length,
        "declination": declination,
        "tmpisr": 1367
        * (
            1.00011
            + 0.034221 * cos(xday)
            + 0.00128 * sin(xday)
            + 0.000719 * cos(2 * xday)
            + 0.000077 * sin(2 * xday)
        ),
    }
