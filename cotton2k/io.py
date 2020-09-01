from locale import atof
from pathlib import Path
from typing import List, Optional

from cotton2k.profile import Profile


def read_profile_file(path: Path) -> Profile:
    """Read the profile file.

    TODO: Maybe profile should be a object instead of dict
    TODO: Profile is self defined file format and not readable for human,
    maybe JSON and TOML are better alternatives
    """
    if not path.exists():
        raise FileNotFoundError(f"{path} not found!")
    return Profile.from_pro(path)


def read_calibration_data(
    var_number: int,
    site_number: int,
    varlist: Optional[Path] = None,  # pylint: disable=E1136
    sitelist: Optional[Path] = None,  # pylint: disable=E1136
):
    """Read the values of the calibration parameters from input files.

    It is called from ReadInput(). It calls GetLineData().

    The following global or file-scope variables are set here:
    SiteName, SitePar, VarName, VarPar

    TODO: Maybe JSON or CSV is more suitable file format than self defined DAT files.
    """
    vars_dir = Path(".")
    site_dir = Path(".")
    varlist = varlist or (vars_dir / "varlist.dat")
    sitelist = sitelist or (site_dir / "sitelist.dat")
    vars_dir = varlist.parent
    site_dir = sitelist.parent
    var_name, var_file = parse_list_dat(varlist.read_text())[var_number]
    if not (var_file_path := vars_dir / var_file).exists():
        raise FileNotFoundError(f"{var_file_path} not found!")
    var_par = parse_parameter(var_file_path.read_text(), 60)
    site_name, site_file = parse_list_dat(sitelist.read_text())[site_number]
    if not (site_file_path := site_dir / site_file).exists():
        raise FileNotFoundError(f"{site_file_path} not found!")
    site_par = parse_parameter(site_file_path.read_text(), 20)
    return dict(siteName=site_name, sitePar=site_par, varName=var_name, varPar=var_par,)


def parse_list_dat(content: str) -> dict:
    result = {}
    for var in content.splitlines():
        k = int(var[:4].strip())
        v = (var[5:25].strip(), var[40:60].strip())
        result[k] = v
    return result


def parse_parameter(content: str, number: int) -> List[float]:
    lines = content.splitlines()
    return list(map(lambda line: atof(line[:15]), lines[1 : number + 1]))  # noqa: E203
