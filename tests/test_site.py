from pathlib import Path

import pytest

from cotton2k.site import Site


@pytest.fixture
def site_dat(tmp_path: Path):
    dat = tmp_path / "test_site.dat"
    lines = [
        "Test site",
        "{:>8.1f}".format(1),
        "{:>8.1f}".format(2),
        "{:>8.1f}".format(2.5),
        "{:>8.4f}".format(0.0060),
    ]
    dat.write_text("\n".join(lines))
    return dat


def test_site(site_dat):
    site = Site.from_dat(site_dat)
    assert site[1] == site.wind_start_hours_after_sunrise
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


@pytest.fixture
def invalid_site_dat(tmp_path: Path):
    dat = tmp_path / "test_site.dat"
    lines = [
        "Invalid site",
        "{:>8.1f}".format(3),
    ]
    dat.write_text("\n".join(lines))
    return dat


def test_invalid_site(invalid_site_dat):
    with pytest.raises(ValueError):
        site = Site.from_dat(invalid_site_dat)
