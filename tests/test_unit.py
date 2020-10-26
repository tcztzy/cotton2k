from decimal import Decimal

import pytest

from cotton2k.units import (
    MJ,
    C,
    J,
    Langley,
    LengthUnit,
    MassUnit,
    N,
    TimeUnit,
    cm,
    degree_Celsius,
    degree_Fahrenheit,
    degree_Kelvin,
    kg,
    km,
    m,
    mile,
    s,
)


def test_quantity():
    assert 1 * m * 1000 == 1 * km
    assert 1 * (m * 1000) == 1 * km
    assert 1 * km * 1.609344 == 1 * mile
    assert 1 * m * m == 1 * m ** 2
    assert 1 * mile / Decimal("1.609344") == 1 * km
    assert 1 * km / 1000 == 1 * m
    with pytest.raises(TypeError):
        1 * m * object()
    with pytest.raises(TypeError):
        1 * m / object()


def test_unit():
    with pytest.raises(TypeError, match="can't compare"):
        assert m == degree_Celsius
    唐 = LengthUnit(1.78, 0)  # author's height
    唐.name = "唐梓涯的身高"
    with pytest.raises(TypeError):
        唐.gain = None
    assert (
        N.base_units[TimeUnit] == -2
        and N.base_units[MassUnit] == 1
        and N.base_units[LengthUnit] == 1
    )
    assert (N * 2).gain == 0.5
    assert (N / 2).gain == 2
    assert (N / 0.5).gain == 0.5
    with pytest.raises(TypeError):
        N * object()
    with pytest.raises(TypeError):
        N / object()
    with pytest.raises(ValueError):
        C * C
    with pytest.raises(ValueError):
        C ** 2
    assert isinstance((m ** -1) ** -1, LengthUnit)
    assert isinstance(s ** 2 * s ** -1, TimeUnit)
    assert isinstance(N * s ** 2 / kg, LengthUnit)


def test_temperature_unit():
    zero_degree_Celsius = 0 * degree_Celsius
    assert zero_degree_Celsius.value == 0
    assert zero_degree_Celsius.unit == degree_Celsius
    with pytest.raises(TypeError):
        assert 0 == zero_degree_Celsius
    with pytest.raises(TypeError):
        assert zero_degree_Celsius == 0
    with pytest.raises(TypeError):
        assert degree_Celsius == object()
    assert zero_degree_Celsius == 0 * degree_Celsius
    assert zero_degree_Celsius == 32 * degree_Fahrenheit
    assert 100 * degree_Celsius == 212 * degree_Fahrenheit
    assert 0 * degree_Kelvin == -273.15 * degree_Celsius


def test_solar_radiation_unit():
    assert 1 * MJ / m ** 2 == 1000000 * J / m ** 2
    assert 1 * J / cm ** 2 == 10000 * J / m ** 2
    assert 4.184 * (MJ / m ** 2) == 100 * Langley
