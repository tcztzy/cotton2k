from os import linesep
from pathlib import Path

from pytest import fixture

CONTENT = """test.pro            Test profile
01-MAY-2020    20-APR-2020    15-OCT-2020                        1.000  100  101
test.act                                         1     0.000     0.000  122    0
test.hyd            test.int            test.agi
    40.548    81.296  1013.000         0
    75.000     0.000    40.000         0
        10    10-APR-2020    20-OCT-2020        10    01-JUN-2020    20-OCT-2020
  0  0  1  1  0  1  0  1  1  1  1  1  0  0  0  0  0  1  0  0  0  0  0"""


@fixture
def pro_file(tmp_path: Path):
    pro_dir = tmp_path / "profiles"
    pro_dir.mkdir()
    path = pro_dir / "test.pro"
    path.write_text(CONTENT)
    return path


@fixture
def tmp_file(tmp_path):
    return tmp_path / "tmp"


@fixture
def varlist(tmp_path, test_var):
    content = "1".rjust(4) + " Test var".ljust(36) + test_var.name + "\n"
    content += "2".rjust(4) + " Test var".ljust(36) + "not exist.dat"
    path = tmp_path / "varlist.dat"
    path.write_text(content)
    return path


@fixture
def test_var(tmp_path: Path):
    path = tmp_path / "test_var.dat"
    headline = "Test var"
    lines = (f"{i}".rjust(8) + " " * 7 + f"VARPAR({i+1})" for i in range(60))
    path.write_text(headline + "\n" + "\n".join(lines))
    return path


@fixture
def sitelist(tmp_path: Path, test_site: Path):
    content = "1".rjust(4) + " Test site".ljust(36) + test_site.name + "\n"
    content += "2".rjust(4) + " Not exist file".ljust(36) + "not exist.dat"
    path = tmp_path / "sitelist.dat"
    path.write_text(content)
    return path


@fixture
def test_site(tmp_path: Path):
    path = tmp_path / "test_site.dat"
    headline = "Test site"
    lines = (f"{i}".rjust(8) + " " * 7 + f"SITEPAR({i+1})" for i in range(20))
    path.write_text(headline + "\n" + "\n".join(lines))
    return path


@fixture
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


@fixture
def invalid_site_dat(tmp_path: Path):
    dat = tmp_path / "test_site.dat"
    lines = [
        "Invalid site",
        "{:>8.1f}".format(3),
    ]
    dat.write_text("\n".join(lines))
    return dat


@fixture
def weather_file(tmp_path: Path) -> Path:
    file_path = tmp_path / "test.act"
    lines = [
        f"{'Whatever':<30}{1:>3}{1:>3}{1:>3}{1:>3}{1:>3}{' '*15}{0:>10.2f}",
        f"{92:>4}{'01-APR-2020':^17}{20:>7.2f}{13.4:>7.2f}{4:>7.2f}{0:>7.2f}{110:>7.2f}{10:>7.2f}",
    ]
    file_path.write_bytes(linesep.join(lines).encode())
    return file_path
