"""Unit Conversion"""

from __future__ import annotations

from decimal import Decimal
from dataclasses import dataclass
from typing import Union

Number = Union[float, int, Decimal]  # pylint: disable=unsubscriptable-object


@dataclass
class Quantity:
    """Quantity class"""

    value: Number
    unit: TemperatureUnit

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Quantity):
            raise TypeError
        if self.unit == o.unit:
            return self.value == o.value
        else:
            self_value = (self.value - self.unit.offset) / self.unit.gain
            o_value = (o.value - o.unit.offset) / o.unit.gain
            return self_value == o_value


@dataclass
class TemperatureUnit:
    """TemperatureUnit class"""

    gain: Number
    offset: Number

    def __rmul__(self, other):
        return Quantity(other, self)


    def __eq__(self, o: object) -> bool:
        if not isinstance(o, TemperatureUnit):
            raise TypeError
        return self.gain == o.gain and self.offset == o.offset

degree_Celsius = TemperatureUnit(1, Decimal("-273.15"))
degree_Fahrenheit = TemperatureUnit(Decimal("1.8"), Decimal("-459.67"))
