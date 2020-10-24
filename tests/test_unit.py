import pytest

from cotton2k.units import (
    LengthUnit,
    degree_Celsius,
    degree_Fahrenheit,
    degree_Kelvin,
    meter,
)


def test_unit():
    with pytest.raises(TypeError, match="can't compare"):
        assert meter == degree_Celsius
    唐 = LengthUnit(1.78, 0)  # author's height
    唐.name = "唐梓涯的身高"
    with pytest.raises(TypeError):
        唐.gain = None


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
