from cotton2k.units import degree_Celsius


def test_temperature_unit():
    zero_degree_Celsius = 0 * degree_Celsius
    assert zero_degree_Celsius.value == 0
    assert zero_degree_Celsius.unit == degree_Celsius
