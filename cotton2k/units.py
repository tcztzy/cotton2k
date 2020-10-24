"""Unit Conversion"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from typing import Any, Union

Number = Union[int, Decimal, Fraction]  # pylint: disable=unsubscriptable-object
AnyNumber = Union[float, Number]  # pylint: disable=unsubscriptable-object


@dataclass
class Quantity:
    """Quantity class"""

    value: Number
    unit: Unit

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Quantity):
            raise TypeError
        if self.unit == o.unit:
            return self.value == o.value
        self_value = (self.value - self.unit.offset) / self.unit.gain
        o_value = (o.value - o.unit.offset) / o.unit.gain
        return self_value == o_value


@dataclass
class Unit:
    """Unit class"""

    gain: Number = 1
    offset: Number = 0

    def __setattr__(self, name: str, value: Any) -> None:
        if name not in ["gain", "offset"]:
            ...
        elif not isinstance(value, (int, float, Decimal)):
            raise TypeError
        elif isinstance(value, float):
            value = Decimal(str(value))
        super().__setattr__(name, value)

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Unit):
            raise TypeError
        if self.__class__ != o.__class__:
            raise TypeError(f"{type(self)} can't compare with {type(o)}")
        return self.gain == o.gain and self.offset == o.offset

    def __rmul__(self, other) -> Quantity:
        if isinstance(other, float):
            value = Decimal(str(other))
        else:
            value = other
        return Quantity(value, self)


class TemperatureUnit(Unit):  # pylint: disable=too-few-public-methods
    """TemperatureUnit class"""


degree_Kelvin = K = TemperatureUnit()
degree_Celsius = C = TemperatureUnit(1, Decimal("-273.15"))
degree_Fahrenheit = F = TemperatureUnit(Decimal("1.8"), Decimal("-459.67"))


class LengthUnit(Unit):  # pylint: disable=too-few-public-methods
    """LengthUnit class"""


meter = m = LengthUnit()
