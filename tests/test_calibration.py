from os import unlink
from pathlib import Path
from tempfile import mkstemp

import pytest

from cotton2k.io import ROOT_DIR, parse_list_dat, parse_parameter, read_calibration_data


def test_list_dat():
    result = parse_list_dat(
        """   2 GC510                              VPARGC510.DAT
   3 MAXXA                              VPARMAXXA.DAT
   4 ACALA SJ2                          VPARSJ2.DAT
   5 SIVON                              VPARSIVON.DAT
  11 DELTAPINE 61                       VPARDP61.DAT
  12 DELTAPINE 77                       VPARDP77.DAT
  17 新陆早8号                              VPARxlz8.DAT
  18 China Test1                        VPARTest1_xlz8.DAT
  19 China Test2                        VPARTest2_xlz8.DAT
  20 China Test3                        VPARTest3_xlz8.DAT"""
    )
    assert result[2] == ("GC510", "VPARGC510.DAT")
    assert result[3] == ("MAXXA", "VPARMAXXA.DAT")
    assert result[4] == ("ACALA SJ2", "VPARSJ2.DAT")
    assert result[5] == ("SIVON", "VPARSIVON.DAT")
    assert result[11] == ("DELTAPINE 61", "VPARDP61.DAT")
    assert result[12] == ("DELTAPINE 77", "VPARDP77.DAT")
    assert result[17] == ("新陆早8号", "VPARxlz8.DAT")
    assert result[18] == ("China Test1", "VPARTest1_xlz8.DAT")
    assert result[19] == ("China Test2", "VPARTest2_xlz8.DAT")
    assert result[20] == ("China Test3", "VPARTest3_xlz8.DAT")
    result = parse_list_dat(
        """   0 CALIFORNIA, West SJ                USACALSJWEST.DAT   
   1 ISRAEL, Coastal                    ILCOASTAL.DAT
   2 ARIZONA, Phoenix                   USAAZCENTRAL.DAT
   3 ISRAEL, Galil                      ILGALIL.DAT 
   4 ISRAEL, Avdat                      ILAVDAT.DAT 
"""
    )
    assert result[0] == ("CALIFORNIA, West SJ", "USACALSJWEST.DAT")
    assert result[1] == ("ISRAEL, Coastal", "ILCOASTAL.DAT")
    assert result[2] == ("ARIZONA, Phoenix", "USAAZCENTRAL.DAT")
    assert result[3] == ("ISRAEL, Galil", "ILGALIL.DAT")
    assert result[4] == ("ISRAEL, Avdat", "ILAVDAT.DAT")


def test_parameter():
    result = parse_parameter(
        """Test
   0.050       Whatever
   0.049       Yet another""",
        2,
    )
    assert result[0] == 0.050
    assert result[1] == 0.049


@pytest.fixture
def varlist(test_var):
    content = "1".rjust(4) + " Test var".ljust(36) + test_var.name + "\n"
    content +=  "2".rjust(4) + " Test var".ljust(36) + "not exist.dat"
    _, path = mkstemp(text=True)
    path = Path(path)
    path.write_text(content)
    return path


@pytest.fixture
def test_var():
    _, path = mkstemp(dir=ROOT_DIR / "data" / "vars", text=True)
    path = Path(path)
    headline = "Test var"
    lines = (f"{i}".rjust(8) + " " * 7 + f"VARPAR({i+1})" for i in range(60))
    path.write_text(headline + "\n" + "\n".join(lines))
    return path


@pytest.fixture
def sitelist(test_site):
    content = "1".rjust(4) + " Test site".ljust(36) + test_site.name + "\n"
    content += "2".rjust(4) + " Not exist file".ljust(36) + "not exist.dat"
    _, path = mkstemp(text=True)
    path = Path(path)
    path.write_text(content)
    return path


@pytest.fixture
def test_site():
    _, path = mkstemp(dir=ROOT_DIR / "data" / "site", text=True)
    path = Path(path)
    headline = "Test site"
    lines = (f"{i}".rjust(8) + " " * 7 + f"SITEPAR({i+1})" for i in range(20))
    path.write_text(headline + "\n" + "\n".join(lines))
    return path


def test_read_calibration_data(
    varlist: Path, test_var: Path, sitelist: Path, test_site: Path
):
    read_calibration_data(1, 1, varlist, sitelist)
    with pytest.raises(FileNotFoundError):
        read_calibration_data(2, 1, varlist, sitelist)
    with pytest.raises(FileNotFoundError):
        read_calibration_data(1, 2, varlist, sitelist)
    unlink(varlist)
    unlink(test_var)
    unlink(sitelist)
    unlink(test_site)
