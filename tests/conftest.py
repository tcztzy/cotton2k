from os import linesep
from pathlib import Path

from pytest import fixture


@fixture
def CONTENT() -> str:
    return """test.pro            Test profile
01-MAY-2020    20-APR-2020    15-OCT-2020                        1.000  100  101
test.act                                         1     0.000     0.000  122    0
test.hyd            test.int            test.agi
    40.548    81.296  1013.000         1
    75.000     0.000    40.000         0
        10    10-APR-2020    20-OCT-2020        10    01-JUN-2020    20-OCT-2020
  0  0  1  1  0  1  0  1  1  1  1  1  0  0  0  0  0  1  0  0  0  0  0"""


@fixture
def pro_file(tmp_path: Path, CONTENT: str, sitelist: Path) -> Path:
    pro_dir = tmp_path / "profiles"
    pro_dir.mkdir()
    path = pro_dir / "test.pro"
    path.write_text(CONTENT)
    return path


@fixture
def tmp_file(tmp_path: Path) -> Path:
    return tmp_path / "tmp"


@fixture
def data_dir(tmp_path: Path) -> Path:
    path = tmp_path / "data"
    path.mkdir()
    return path


@fixture
def vars_dir(data_dir: Path) -> Path:
    path = data_dir / "vars"
    path.mkdir()
    return path


@fixture
def site_dir(data_dir: Path) -> Path:
    path = data_dir / "site"
    path.mkdir()
    return path


@fixture
def varlist(test_var: Path) -> Path:
    content = "1".rjust(4) + " Test var".ljust(36) + test_var.name + "\n"
    content += "2".rjust(4) + " Test var".ljust(36) + "not exist.dat"
    path = test_var.parent / "varlist.dat"
    path.write_text(content)
    return path


@fixture
def test_var(vars_dir: Path) -> Path:
    path = vars_dir / "test_var.dat"
    headline = "Test var"
    lines = (f"{i}".rjust(8) + " " * 7 + f"VARPAR({i+1})" for i in range(60))
    path.write_text(headline + "\n" + "\n".join(lines))
    return path


@fixture
def sitelist(site_dat: Path) -> Path:
    content = "1".rjust(4) + " Test site".ljust(36) + site_dat.name + "\n"
    content += "2".rjust(4) + " Not exist file".ljust(36) + "not exist.dat"
    path = site_dat.parent / "sitelist.dat"
    path.write_text(content)
    return path


@fixture
def site_dat(site_dir: Path) -> Path:
    dat = site_dir / "test_site.dat"
    lines = [
        "Test site",
        "{:>8.1f}".format(1),
        "{:>8.1f}".format(2),
        "{:>8.1f}".format(2.5),
        "{:>8.4f}".format(0.0060),
        "{:>8.3f}".format(6.264),
        "{:>8.3f}".format(14.554),
        "{:>8.0f}".format(60),
        "{:>8.1f}".format(3.5),
        "{:>8f}".format(20),
        "{:>8f}".format(4),
        "{:>8.0f}".format(180),
    ]
    dat.write_text("\n".join(lines))
    return dat


@fixture
def invalid_site_dat(tmp_path: Path) -> Path:
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