"""Input/Output"""
import csv
import datetime
import json
from pathlib import Path
from typing import Union

from _cotton2k import (  # pylint: disable=import-error, no-name-in-module
    Climate,
    Simulation,
    SoilImpedance,
    SoilInit,
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
    if isinstance(path, dict):
        kwargs = path
    else:
        kwargs = json.loads(Path(path).read_text())
    sim = Simulation(kwargs.get("version", 0x0400))
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


def write_output(  # pylint: disable=too-many-locals
    sim: Simulation, output_dir: Path = Path("."), name: str = "default"
) -> None:
    path = (
        output_dir
        / name
        / (
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            + f" v{sim.version // 0x100:x}.{sim.version % 0x100:x}"
        )
    )
    if not path.exists():
        path.mkdir(parents=True)
    for state in sim.states:
        doy = state["daynum"]
        state_dir = path / datetime.datetime.strptime(
            f"{sim.year} {doy}", "%Y %j"
        ).strftime("%Y-%m-%d")
        if not state_dir.exists():
            state_dir.mkdir(parents=True)
        vegetative_branches = state["vegetative_branches"]
        with open(state_dir / "index.json", "w") as f:
            json.dump(
                {
                    key: state[key]
                    for key in state.keys()
                    if "vegetative_branches" not in key
                },
                f,
            )
        for i in range(state["number_of_vegetative_branches"]):
            vb_dir = state_dir / f"vegetative_branch{i}"
            fruiting_branches = vegetative_branches[i]["fruiting_branches"]
            if not vb_dir.exists():
                vb_dir.mkdir(parents=True)
            with open(vb_dir / "index.json", "w") as f:
                json.dump(
                    {
                        key: vegetative_branches[i][key]
                        for key in vegetative_branches[i]
                        if "fruiting_branches" not in key
                    },
                    f,
                )
            for j in range(vegetative_branches[i]["number_of_fruiting_branches"]):
                fb_dir = vb_dir / f"fruiting_branch{j}"
                if not fb_dir.exists():
                    fb_dir.mkdir(parents=True)
                with open(fb_dir / "index.json", "w") as f:
                    json.dump(
                        {
                            key: fruiting_branches[j][key]
                            for key in fruiting_branches[j]
                            if "nodes" not in key
                        },
                        f,
                    )
                for k in range(fruiting_branches[j]["number_of_fruiting_nodes"]):
                    fruiting_nodes = fruiting_branches[j]["nodes"]
                    with open(fb_dir / f"node{k}.json", "w") as f:
                        json.dump(fruiting_nodes[k], f)
