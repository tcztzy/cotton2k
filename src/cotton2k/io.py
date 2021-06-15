"""Input/Output"""
import csv
import datetime
import json
from operator import attrgetter
from pathlib import Path
from typing import Union

from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import Session

from _cotton2k import (  # pylint: disable=import-error, no-name-in-module
    Climate,
    Simulation,
    SoilImpedance,
    SoilInit,
)
from cotton2k.models import Simulation as SimulationModel
from cotton2k.models import State as StateModel

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
    sim = Simulation(kwargs.get("name", "default"), kwargs.get("version", 0x0400))
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


def prepare_database(engine_url):
    engine = create_engine(engine_url)
    for model in (SimulationModel, StateModel):
        if not inspect(engine).has_table(model.__tablename__):
            model.__table__.create(bind=engine, checkfirst=True)
    session = Session(bind=engine)
    return session


def write_output(
    simulation: Simulation,
    engine_url: str = "sqlite+pysqlite:///:memory:",
):
    session = prepare_database(engine_url)
    simulation_model = SimulationModel(
        version=simulation.version,
        execute_time=datetime.datetime.now(),
        name=simulation.name,
    )
    session.add(simulation_model)
    session.commit()
    for state in simulation.states:
        state_model = StateModel(
            simulation_id=simulation_model.id,
            date=datetime.datetime.strptime(
                f"{simulation.year} {state.daynum}", "%Y %j"
            ).date(),
            **{
                key: state[key]
                for key in (
                    "kday",
                    "lint_yield",
                    "ginning_percent",
                    "plant_height",
                    "plant_weight",
                    "stem_weight",
                    "square_weight",
                    "green_bolls_weight",
                    "open_bolls_weight",
                    "leaf_weight",
                    "leaf_area_index",
                    "light_interception",
                    "number_of_squares",
                    "number_of_green_bolls",
                    "number_of_open_bolls",
                )
            },
            number_of_fruiting_branches=sum(
                map(
                    attrgetter("number_of_fruiting_branches"), state.vegetative_branches
                )
            ),
        )
        session.add(state_model)
        session.commit()
    return session, simulation_model
