from pathlib import Path

import pytest

from cotton2k.site import Site

from .fixtures import data_dir, invalid_site_dat, site_dat, site_dir


def test_site(site_dat):
    site = Site.from_dat(site_dat)
    assert site[1] == site.wind_start_hours_after_sunrise
    assert site.parameters[1] == site[1]
    assert site[1] == 1.00
    assert site[2] == 2.00
    assert site[3] == 2.5
    assert site[4] == 0.0060
    with pytest.raises(KeyError):
        site[0]
    site[1] = 2.00
    assert site._site_par[0] == 2.00 and site[1] == 2.00
    with pytest.raises(ValueError):
        site[1] = 3.00
    site.deep_soil_temperature(10)


def test_invalid_site(invalid_site_dat):
    with pytest.raises(ValueError):
        site = Site.from_dat(invalid_site_dat)
