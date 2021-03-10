"""Input/Output"""
import json
from locale import atof
from pathlib import Path

from _cotton2k import _Simulation  # pylint: disable=import-error# noqa: F401


def read_calibration_data(var_number: int, varlist: Path, type_: str = "var"):
    """Read the values of the calibration parameters from input files.

    It is called from ReadInput(). It calls GetLineData().

    The following global or file-scope variables are set here:
    SiteName, SitePar, VarName, VarPar

    TODO: Maybe JSON or CSV is more suitable file format than self defined DAT files.
    """
    vars_dir = varlist.parent
    var_name, var_file = parse_list_dat(varlist.read_text())[var_number]
    if not (var_file_path := vars_dir / var_file).exists():
        raise FileNotFoundError(f"{var_file_path} not found!")
    var_par = parse_parameter(var_file_path.read_text(), 60 if type_ == "var" else 16)
    return {type_ + "Name": var_name, type_ + "Par": var_par}


def parse_list_dat(content: str) -> dict:
    """Parse the culvar and site list .dat file"""
    result = {}
    for var in content.splitlines():
        k = int(var[:4].strip())
        values = (var[5:25].strip(), var[40:60].strip())
        result[k] = values
    return result


def parse_parameter(content: str, number: int) -> list[float]:  # type: ignore
    """Parse the parameter from .dat file"""
    lines = content.splitlines()
    return list(map(lambda line: atof(line[:15]), lines[1 : number + 1]))  # noqa: E203


def read_input(path: Path) -> _Simulation:
    sim = _Simulation()
    kwargs = json.loads(path.read_text())
    for attr in [
        "start_date",
        "stop_date",
        "emerge_date",
        "plant_date",
    ]:
        if attr in kwargs:
            setattr(sim, attr, kwargs.get(attr))
    sim.year = int(kwargs["start_date"][:4])
    sim.read_input(**kwargs)
    return sim
