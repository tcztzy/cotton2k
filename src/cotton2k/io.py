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
    SoilImpedance,
    SoilInit,
)
from cotton2k.models import AgronomyOperation, AgronomyOperationType
from cotton2k.models import Climate as ClimateModel
from cotton2k.models import Cultivar, Profile
from cotton2k.models import Simulation as SimulationModel
from cotton2k.models import Site, Soil, SoilHydrology
from cotton2k.models import State as StateModel
from cotton2k.models import association_table
from cotton2k.simulation import Simulation

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


def prepare_database(engine_url):
    engine = create_engine(engine_url)
    for model in (
        AgronomyOperation,
        ClimateModel,
        Cultivar,
        Profile,
        Soil,
        SoilHydrology,
        SimulationModel,
        Site,
        StateModel,
    ):
        if not inspect(engine).has_table(model.__tablename__):
            model.__table__.create(bind=engine, checkfirst=True)
    if not inspect(engine).has_table("profile_agronomy_operation"):
        association_table.create(bind=engine, checkfirst=True)
    return engine


def read_input(path: Union[Path, str, dict, Profile], session=None) -> Simulation:
    # pylint: disable=line-too-long
    if isinstance(path, Profile):
        profile = path
        topping = (
            session.query(AgronomyOperation)
            .filter(AgronomyOperation.profiles.any(Profile.id == profile.id))  # type: ignore
            .filter(
                AgronomyOperation.profiles.any(  # type: ignore
                    AgronomyOperation.type == AgronomyOperationType.topping
                )
            )
            .one_or_none()
        )
        plant = (
            session.query(AgronomyOperation)
            .filter(AgronomyOperation.profiles.any(Profile.id == profile.id))  # type: ignore
            .filter(
                AgronomyOperation.profiles.any(  # type: ignore
                    AgronomyOperation.type == AgronomyOperationType.planting
                )
            )
            .one()
        )
        kwargs = {
            "id": profile.id,
            "start_date": profile.start_date,
            "stop_date": profile.stop_date,
            "emerge_date": profile.emerge_date,
            "plant_date": plant.date,
            "latitude": profile.site.latitude,
            "longitude": profile.site.longitude,
            "elevation": profile.site.elevation,
            "site_parameters": [
                attrgetter(f"param{i + 1:0>2}")(profile.site) for i in range(16)
            ],
            "cultivar_parameters": [
                attrgetter(f"param{i + 1:0>2}")(profile.cultivar) for i in range(50)
            ],
            "row_space": plant.planting_row_width,
            "skip_row_width": plant.planting_skip_row_width,
            "plants_per_meter": plant.planting_plants_per_meter,
            "soil": {
                "initial": [
                    {
                        "ammonium_nitrogen": layer.ammonium_nitrogen,
                        "nitrate_nitrogen": layer.nitrate_nitrogen,
                        "organic_matter": layer.organic_matter,
                        "water": layer.water,
                    }
                    for layer in session.query(Soil)
                    .where(Site.id == profile.site.id)
                    .order_by(Soil.depth)
                ],
                "hydrology": {
                    "ratio_implicit": profile.site.soil_implicit_solution_ratio,
                    "max_conductivity": profile.site.max_conductivity,
                    "field_capacity_water_potential": profile.site.field_capacity_water_potential,
                    "immediate_drainage_water_potential": profile.site.immediate_drainage_water_potential,
                    "layers": [
                        {
                            "depth": layer.depth,
                            "air_dry": layer.air_dry_water_content,
                            "theta": layer.saturated_water_content,
                            "alpha": layer.alpha,
                            "beta": layer.beta,
                            "saturated_hydraulic_conductivity": layer.saturated_hydraulic_conductivity,
                            "field_capacity_hydraulic_conductivity": layer.field_capacity_hydraulic_conductivity,
                            "bulk_density": layer.bulk_density,
                            "clay": layer.clay,
                            "sand": layer.sand,
                        }
                        for layer in session.query(SoilHydrology)
                        .where(Site.id == profile.site.id)
                        .order_by(SoilHydrology.depth)
                    ],
                },
            },
            "climate_start_date": profile.start_date,
            "climate": [
                {
                    "min": c.min_temperature,
                    "max": c.max_temperature,
                    "rain": c.precipitation,
                    "wind": c.wind_speed,
                    "radiation": c.radiation,
                    "dewpoint": c.dewpoint,
                }
                for c in profile.site.climate
            ],
        }
        if topping is not None:
            kwargs["topping_date"] = topping.date
    elif isinstance(path, dict):
        kwargs = path
    else:
        kwargs = json.loads(Path(path).read_text())
    sim = Simulation(kwargs.get("id", 0), kwargs.get("version", 0x0400), **kwargs)  # type: ignore[arg-type]
    soil = SoilInit(**kwargs.get("soil", {}))  # type: ignore[arg-type]
    start_date = kwargs["start_date"]
    if not isinstance(start_date, (datetime.date, str)):
        raise ValueError
    sim.year = (
        start_date.year
        if isinstance(start_date, datetime.date)
        else int(start_date[:4])
    )
    sim.read_input(lyrsol=soil.lyrsol, **kwargs)
    climate_start_date = kwargs.get("climate_start_date", 0)
    sim.climate = Climate(climate_start_date, kwargs.get("climate"))[sim.start_date :]  # type: ignore[misc]
    return sim


def write_output(
    simulation: Simulation,
    engine_url: str = "sqlite+pysqlite:///:memory:",
    session=None,
    start_at=None,
    execute_time=None,
):
    engine = prepare_database(engine_url)
    if session is None:
        session = Session(engine)
    simulation_model = SimulationModel(
        version=simulation.version,
        start_at=start_at,
        execute_time=execute_time,
        profile_id=simulation.profile_id or None,
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
