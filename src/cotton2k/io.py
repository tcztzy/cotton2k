"""Input/Output"""
import csv
import json
from pathlib import Path
from typing import Optional, Union

from _cotton2k import (  # pylint: disable=import-error, no-name-in-module
    Climate,
    FruitingBranch,
    Simulation,
    SoilImpedance,
    SoilInit,
    State,
    VegetativeBranch,
)

SOIL_IMPEDANCE = SoilImpedance()
with open(Path(__file__).parent / "soil_imp.csv") as csvfile:
    reader = csv.DictReader(csvfile)
    SOIL_IMPEDANCE.curves = list(
        map(
            lambda row: {
                (k if k == "water" else float(k)): float(v) for k, v in row.items()
            },
            reader,
        )
    )


def read_input(path: Union[Path, str, dict]) -> Simulation:
    sim = Simulation()
    if isinstance(path, dict):
        kwargs = path
    else:
        kwargs = json.loads(Path(path).read_text())
    for attr in [
        "start_date",
        "stop_date",
        "emerge_date",
        "plant_date",
        "topping_date",
        "latitude",
        "longitude",
        "elevation",
        "site_parameters",
        "cultivar_parameters",
        "row_space",
        "skip_row_width",
        "plants_per_meter",
    ]:
        if attr in kwargs:
            setattr(sim, attr, kwargs.get(attr))
    soil = SoilInit(**kwargs.get("soil", {}))
    sim.year = int(kwargs["start_date"][:4])
    sim.read_input(lyrsol=soil.lyrsol, **kwargs)
    sim.climate = Climate(kwargs.get("climate_start_date", 0), kwargs.get("climate"))[
        sim.start_date :
    ]
    return sim


class Cotton2KJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (State, VegetativeBranch, FruitingBranch)):
            return dict(o)
        return super.default(o)


def write_output(sim: Simulation, path: Optional[Path] = None) -> str:
    content = json.dumps(sim.states, cls=Cotton2KJSONEncoder)
    if path is not None:
        path.write_text(content)
    return content
