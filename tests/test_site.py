from pathlib import Path

import pytest

from cotton2k.site import Site


@pytest.fixture
def site_dat(tmp_path: Path):
    dat = tmp_path / "test_site.dat"
    lines = [
        "Test site",
        "{:>8.1f}".format(1),
        ""
    ]
    dat.write_text("\n".join(lines))
    return dat


def test_site(site_dat):
    site = Site.from_dat(site_dat)
    assert site[1] == site.wind_start_hours_after_sunraise
    assert site[1] == 1.00
