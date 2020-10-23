import pytest

from cotton2k.units import degree_Celsius, degree_Fahrenheit


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
