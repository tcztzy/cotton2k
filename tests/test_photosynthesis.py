from _cotton2k.photosynthesis import ambient_co2_factor
from cotton2k.simulation import compute_light_interception


def test_ambient_co2_factor():
    assert ambient_co2_factor(1959) == 1
    assert ambient_co2_factor(2019) == 1.3052599999999999
    assert ambient_co2_factor(2020) == 1.310124
    assert ambient_co2_factor(2021) == 1.314988


def test_compute_light_interception():
    nan = float("nan")
    assert (
        compute_light_interception(1, nan, nan, nan, version=0x0500)
        == 0.6865138191173947
    )
