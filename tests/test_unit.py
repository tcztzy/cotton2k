import pytest

from cotton2k.units import (
    C,
    LengthUnit,
    MassUnit,
    N,
    TimeUnit,
    degree_Celsius,
    degree_Fahrenheit,
    degree_Kelvin,
    kg,
    meter,
    s,
)


def test_unit():
    with pytest.raises(TypeError, match="can't compare"):
        assert meter == degree_Celsius
    唐 = LengthUnit(1.78, 0)  # author's height
    唐.name = "唐梓涯的身高"
    with pytest.raises(TypeError):
        唐.gain = None
    assert (
        N.base_units[TimeUnit] == -2
        and N.base_units[MassUnit] == 1
        and N.base_units[LengthUnit] == 1
    )
    with pytest.raises(TypeError):
        N * 2
    with pytest.raises(TypeError):
        N / 2
    with pytest.raises(ValueError):
        C * C
    with pytest.raises(ValueError):
        C ** 2
    assert isinstance((meter ** -1) ** -1, LengthUnit)
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
