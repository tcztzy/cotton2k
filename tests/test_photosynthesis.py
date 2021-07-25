from _cotton2k.photosynthesis import ambient_co2_factor


def test_ambient_co2_factor():
    assert ambient_co2_factor(1959) == 1
    assert ambient_co2_factor(2019) == 1.3052599999999999
    assert ambient_co2_factor(2020) == 1.310124
    assert ambient_co2_factor(2021) == 1.314988
