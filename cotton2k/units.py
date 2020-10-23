"""Unit Conversion"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

Number = Union[float, int]  # pylint: disable=unsubscriptable-object


@dataclass
class Quantity:
    """Quantity class"""

    value: Number
    unit: TemperatureUnit


@dataclass
class TemperatureUnit:
    """TemperatureUnit class"""

    gain: Number
    offset: Number

    def __rmul__(self, other):
        return Quantity(other, self)


degree_Celsius = TemperatureUnit(1, 273.15)
