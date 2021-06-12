import datetime
from calendar import isleap
from math import acos, cos, degrees, radians, sin, tan

from _cotton2k.utils import date2doy

TZ_WIDTH = 15  # timezone width in degree


def compute_day_length(coordinate: tuple[float, float], date: datetime.date) -> dict:
    lat, lon = coordinate
    xday: float = radians(360 * date2doy(date) / (365 + int(isleap(date.year))))
    declination: float = (
        0.006918
        - 0.399912 * cos(xday)
        + 0.070257 * sin(xday)
        - 0.006758 * cos(2 * xday)
        + 0.000907 * sin(2 * xday)
        - 0.002697 * cos(3 * xday)
        + 0.001480 * sin(3 * xday)
    )
    exday: datetime.timedelta = datetime.timedelta(
        hours=degrees(
            0.000075
            + 0.001868 * cos(xday)
            - 0.032077 * sin(xday)
            - 0.014615 * cos(2 * xday)
            - 0.04089 * sin(2 * xday)
        )
        / TZ_WIDTH
    )

    solar_noon: datetime.datetime = (
        datetime.datetime.combine(
            date,
            datetime.time(12),
            datetime.timezone(datetime.timedelta(hours=lon / TZ_WIDTH)),
        ).astimezone(
            datetime.timezone(
                datetime.timedelta(hours=(lon + TZ_WIDTH / 2) // TZ_WIDTH)
            )
        )
        - exday
    )
    ht: float = -tan(radians(lat)) * tan(declination)
    if ht > 1:
        ht = 1
    elif ht < -1:
        ht = -1

    day_length = datetime.timedelta(hours=degrees(2 * acos(ht)) / TZ_WIDTH)
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


def radiation(radsum: float, sinb: float, c11: float) -> float:
    """
    Function radiation() computes the hourly values of global radiation, in W m-2,
    using the measured daily total global radiation.

    The algorithm follows the paper of Spitters et al. (1986). It assumes
    that atmospheric transmission of radiation is lower near the margins of
    the daylight period, because of an increase in the path length through
    the atmosphere at lower solar heights. Radiation is therefore assumed to be
    proportional to sinb * (1 + c11 * sinb), where the value of c11 is set as 0.4 .

    Input arguments:
    radsum - daily radiation integral.
    sinb - sine of the solar elevation.
    c11 - constant parameter (0.4).

    References:

    Spitters, C.J.T., Toussaint, H.A.J.M. and Goudriaan, J. 1986.
    Separating the diffuse and direct component of global radiation and
    its implications for modeling canopy photosynthesis. Part I.
    Components of incoming radiation. Agric. For. Meteorol. 38:217-229.

    Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    diurnal patterns of air temperature, radiation, wind speed and
    relative humidity by equations from daily characteristics.
    Agricultural Systems 51:377-393.
    """
    return 0 if sinb <= 0 else radsum * sinb * (1 + c11 * sinb)
