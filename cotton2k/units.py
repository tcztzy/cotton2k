"""Unit Conversion"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from decimal import Decimal
from typing import Any, Type, Union

Number = Union[int, Decimal]  # pylint: disable=unsubscriptable-object
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
    # https://github.com/python/mypy/pull/9564
    base_units: Counter[Type[Unit]] = field(init=False)  # type: ignore

    def __init_subclass__(cls) -> None:
        cls.base_units = Counter({cls: 1})

    def __setattr__(self, name: str, value: Any) -> None:
        if name not in ["gain", "offset"]:
            ...
        elif isinstance(value, float):
            value = Decimal(str(value))
        elif not isinstance(value, (int, float, Decimal)):
            raise TypeError
        super().__setattr__(name, value)

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Unit):
            raise TypeError
        if self.__class__ != o.__class__:
            raise TypeError(f"{type(self)} can't compare with {type(o)}")
        return (
            self.gain == o.gain
            and self.offset == o.offset
            and self.base_units == o.base_units
        )

    @staticmethod
    def __simplify__(base_units):
        values = list(base_units.values())
        if values.count(1) == 1 and values.count(0) == len(base_units) - 1:
            return base_units.most_common(1)[0][0]
        return None

    def __mul__(self, other):
        if not isinstance(other, (Unit, int, float, Decimal)):
            raise TypeError
        if isinstance(other, (int, float, Decimal)):
            if isinstance(other, float):
                other = Decimal(str(other))
            unit = Unit(self.gain * other)
            unit.base_units = self.base_units
            return unit
        base_units = self.base_units.copy()
        base_units.update(other.base_units)
        if self.offset != 0 or other.offset != 0:
            raise ValueError
        gain = self.gain * other.gain
        if cls := self.__simplify__(base_units):
            return cls(gain)
        unit = Unit(gain)
        unit.base_units = base_units
        return unit

    def __truediv__(self, other):
        if not isinstance(other, (Unit, int, float, Decimal)):
            raise TypeError
        if isinstance(other, (int, float, Decimal)):
            if isinstance(other, float):
                other = Decimal(str(other))
            unit = Unit(self.gain / other)
            unit.base_units = self.base_units
            return unit
        base_units = self.base_units.copy()
        base_units.subtract(other.base_units)
        gain = self.gain * other.gain
        if cls := self.__simplify__(base_units):
            return cls(gain)
        unit = Unit(gain)
        unit.base_units = base_units
        return unit

    def __rmul__(self, other) -> Quantity:
        if isinstance(other, float):
            value = Decimal(str(other))
        else:
            value = other
        return Quantity(value, self)

    def __pow__(self, other):
        base_units = Counter({u: v.__mul__(other) for u, v in self.base_units.items()})
        if self.offset != 0 or not isinstance(other, int):
            raise ValueError
        gain = self.gain.__pow__(other)
        if cls := self.__simplify__(base_units):
            return cls(gain)
        unit = Unit(gain)
        unit.base_units = base_units
        return unit


class LengthUnit(Unit):  # pylint: disable=too-few-public-methods
    """LengthUnit class"""


class MassUnit(Unit):  # pylint: disable=too-few-public-methods
    """WeightUnit class"""


class TimeUnit(Unit):  # pylint: disable=too-few-public-methods
    """TimeUnit class"""


class TemperatureUnit(Unit):  # pylint: disable=too-few-public-methods
    """TemperatureUnit class"""


meter = m = LengthUnit()
kilo_gram = kg = MassUnit()
second = s = TimeUnit()

degree_Kelvin = K = TemperatureUnit()
degree_Celsius = C = TemperatureUnit(1, Decimal("-273.15"))
degree_Fahrenheit = F = TemperatureUnit(Decimal("1.8"), Decimal("-459.67"))

# Force
N = kg * m / s ** 2

# Energy
J = Joule = N * m
Cal = Calorie = J * 0.238846
