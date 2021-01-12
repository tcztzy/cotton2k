"""Climate module"""
from __future__ import annotations

import datetime
from calendar import isleap
from dataclasses import dataclass
from locale import atof, atoi
from math import acos, cos, exp, pi, radians, sin, tan
from pathlib import Path
from re import findall
from typing import Any, Union

from cotton2k.utils import date_to_day_of_year


def read_climate_data(climate_file: Path) -> Climate:
    """Read climate data

    :param climate_file: the climate file path
    :type climate_file: Path
    """
    return Climate.from_file(climate_file)


def compute_day_length(date: datetime.date, latitude: float, longtitude: float, /):
    """Can get improved by importing astral package"""
    # pylint: disable=C0401
    # Convert day of year to corresponding angle in radians (xday).
    xday = 2 * pi * (date_to_day_of_year(date) - 1) / (365 + int(isleap(date.year)))
    # Compute declination angle for this day. The equation used here for
    # computing it is taken from the CIMIS algorithm.
    declination = (
        0.006918
        - 0.399912 * cos(xday)
        + 0.070257 * sin(xday)
        - 0.006758 * cos(2 * xday)
        + 0.000907 * sin(2 * xday)
        - 0.002697 * cos(3 * xday)
        + 0.001480 * sin(3 * xday)
    )
    # Compute extraterrestrial radiation in W  m-2. The 'solar constant'
    # (average value = 1367 W m-2) is corrected for this day's distance
    # between earth and the sun. The equation used here is from the CIMIS
    # algorithm, which is based on the work of Iqbal (1983).
    tmpisr = 1367 * (
        1.00011
        + 0.034221 * cos(xday)
        + 0.00128 * sin(xday)
        + 0.000719 * cos(2 * xday)
        + 0.000077 * sin(2 * xday)
    )

    # Time of solar noon (SolarNoon) is computed by the CIMIS algorithm,
    # using a correction for longitude (f), and the date correction (exday).
    # It is assumed that the time zone is "geographically correct". For
    # example,longitude between 22.5 and 37.5 East is in time zone GMT+2, and
    # longitude between 112.5 and 127.5 West is in time zone GMT-8. All daily
    # times in the model are computed by this method.
    local_timezone = datetime.timezone(datetime.timedelta(hours=longtitude / 15))
    exday = datetime.timedelta(
        hours=(
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
    )

    solar_noon = (
        datetime.datetime.combine(
            date,
            datetime.time(12, 0, 0),
            local_timezone,
        )
        - exday
    ).astimezone(datetime.timezone.utc)
    xlat = radians(latitude)
    ht = -tan(xlat) * tan(declination)
    ht = 1 if ht > 1 else -1 if ht < -1 else ht
    day_length = datetime.timedelta(hours=2 * acos(ht) * 12 / pi)
    sunr = solar_noon - day_length / 2
    suns = sunr + day_length
    return day_length, sunr, solar_noon, suns, declination, tmpisr


def parse_weather(content: str) -> dict:
    """Parse weather file"""
    lines = content.splitlines()
    head_line, *daily_climate_lines = lines
    result: dict[str, Any] = dict()  # type: ignore
    n_length = len(head_line)
    if n_length >= 31:
        result["isw_rad"] = bool(atoi(head_line[31:34]))
    if n_length >= 34:
        result["isw_tmp"] = bool(atoi(head_line[34:37]))
    if n_length >= 37:
        result["isw_rain"] = bool(atoi(head_line[37:40]))
    if n_length >= 40:
        result["isw_wind"] = bool(atoi(head_line[40:43]))
    if n_length >= 43:
        result["isw_dewt"] = bool(atoi(head_line[43:46]))
    if n_length >= 61:
        result["average_wind"] = atof(head_line[61:71])
    clim = list()
    for line in daily_climate_lines:
        kwargs: dict[str, Union[bool, float]] = {  # type: ignore
            "_" + k: v for k, v in result.items() if k.startswith("isw")
        }
        (
            kwargs["_rad"],
            kwargs["_tmax"],
            kwargs["_tmin"],
            kwargs["_rain"],
            kwargs["_wind"],
            kwargs["_dewt"],
        ) = map(atof, findall(".{7}", line[21:]))
        c = DailyClimate(**kwargs)  # type: ignore
        clim.append(c)
    result["_climate"] = clim
    return result


def tdewest(t: float, site_parameter5: float, site_parameter6: float) -> float:
    """
    This function estimates the approximate daily average dewpoint temperature
    when it is not available.

    It is called by ReadClimateData().

    Global variables referenced: SitePar[5] and SitePar[6]

    Argument used: t = maximum temperature of this day.
    """
    if t <= 20:
        return site_parameter5
    if t >= 40:
        return site_parameter6
    return ((40 - t) * site_parameter5 + (t - 20) * site_parameter6) / 20


def vapor_pressure(temperature: float) -> float:
    """
    Compute vapor pressure in the air (in KPa units) function of the air at
    temperature (C).

    Tetens, O. 1930. Uber einige meteorologische Begriffe. Z. Geophys.. 6. 297–309.

    Buck, A. L., 1981: New Equations for Computing Vapor Pressure and
    Enhancement Factor. J. Appl. Meteor., 20, 1527–1532,
    https://doi.org/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2.
    """
    return 0.61078 * exp(17.269 * temperature / (temperature + 237.3))


@dataclass
class DailyClimate:  # pylint: disable=too-many-instance-attributes
    """Class represent daily climate"""

    _rad: float
    _tmax: float
    _tmin: float
    _rain: float
    _wind: float
    _dewt: float

    _isw_rad: bool
    _isw_tmp: bool
    _isw_rain: bool
    _isw_wind: bool
    _isw_dewt: bool

    @property
    def radiation(self):
        """Radiation"""
        return self._rad if not self._isw_rad else self._rad * 23.884

    @property
    def max_temperature(self):
        """Max temperature"""
        return self._tmax if self._isw_tmp else (self._tmax - 32) / 1.8

    @property
    def min_temperature(self):
        """Min temperature"""
        return self._tmin if self._isw_tmp else (self._tmin - 32) / 1.8

    @property
    def rain(self):
        """Rainfall"""
        return self._rain if self._isw_rain else self._rain * 25.4

    @property
    def wind(self):
        """Wind speed"""
        return self._wind if self._isw_wind else self._wind * 1.609

    @property
    def dew_temperature(self):
        """Dew point temperature"""
        return self._dewt if self._isw_dewt else (self._dewt - 32) / 1.8


@dataclass
class Climate:
    """Climate class"""

    _climate: list[DailyClimate]  # type: ignore
    isw_rad: bool
    isw_tmp: bool
    isw_rain: bool
    isw_wind: bool
    isw_dewt: bool
    average_wind: float

    @classmethod
    def from_file(cls, path: Path) -> Climate:
        """Create new Climate from file"""
        return cls(**parse_weather(path.read_text()))

    def __getitem__(
        self, item: Union[int, slice]
    ) -> Union[DailyClimate, list[DailyClimate]]:  # type: ignore
        return self._climate[item]