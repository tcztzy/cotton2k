# pylint: disable=too-few-public-methods
import enum

from sqlalchemy import (
    Column,
    Date,
    DateTime,
    Enum,
    Float,
    ForeignKey,
    Integer,
    String,
    Table,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Model = declarative_base(name="Model")


class Climate(Model):
    __tablename__ = "climate"
    id = Column(Integer, primary_key=True)
    site_id = Column(Integer, ForeignKey("site.id"))
    date = Column(Date)
    radiation = Column("ALLSKY_SFC_SW_DWN", Float)
    temperature = Column("T2M", Float)
    max_temperature = Column("T2M_MAX", Float)
    min_temperature = Column("T2M_MIN", Float)
    precipitation = Column("PRECTOT", Float)
    wind_speed = Column("WS2M", Float)
    dewpoint = Column("T2MDEW", Float)
    humidity = Column("RH2M", Float)


class Cultivar(Model):
    __tablename__ = "cultivar"
    id = Column(Integer, primary_key=True)
    name = Column(String)
    param01 = Column(Float)
    param02 = Column(Float)
    param03 = Column(Float)
    param04 = Column(Float)
    param05 = Column(Float)
    param06 = Column(Float)
    param07 = Column(Float)
    param08 = Column(Float)
    param09 = Column(Float)
    param10 = Column(Float)
    param11 = Column(Float)
    param12 = Column(Float)
    param13 = Column(Float)
    param14 = Column(Float)
    param15 = Column(Float)
    param16 = Column(Float)
    param17 = Column(Float)
    param18 = Column(Float)
    param19 = Column(Float)
    param20 = Column(Float)
    param21 = Column(Float)
    param22 = Column(Float)
    param23 = Column(Float)
    param24 = Column(Float)
    param25 = Column(Float)
    param26 = Column(Float)
    param27 = Column(Float)
    param28 = Column(Float)
    param29 = Column(Float)
    param30 = Column(Float)
    param31 = Column(Float)
    param32 = Column(Float)
    param33 = Column(Float)
    param34 = Column(Float)
    param35 = Column(Float)
    param36 = Column(Float)
    param37 = Column(Float)
    param38 = Column(Float)
    param39 = Column(Float)
    param40 = Column(Float)
    param41 = Column(Float)
    param42 = Column(Float)
    param43 = Column(Float)
    param44 = Column(Float)
    param45 = Column(Float)
    param46 = Column(Float)
    param47 = Column(Float)
    param48 = Column(Float)
    param49 = Column(Float)
    param50 = Column(Float)


class Site(Model):
    __tablename__ = "site"
    id = Column(Integer, primary_key=True)
    name = Column(String)
    latitude = Column(Float)
    longitude = Column(Float)
    elevation = Column(Float)
    climate: list["Climate"] = relationship("Climate", backref="site")
    soil_implicit_solution_ratio = Column(Float)
    max_conductivity = Column(Float)
    field_capacity_water_potential = Column(Float)
    immediate_drainage_water_potential = Column(Float)
    param01 = Column(Float)
    param02 = Column(Float)
    param03 = Column(Float)
    param04 = Column(Float)
    param05 = Column(Float)
    param06 = Column(Float)
    param07 = Column(Float)
    param08 = Column(Float)
    param09 = Column(Float)
    param10 = Column(Float)
    param11 = Column(Float)
    param12 = Column(Float)
    param13 = Column(Float)
    param14 = Column(Float)
    param15 = Column(Float)
    param16 = Column(Float)


class AgronomyOperationType(enum.Enum):
    planting = 0
    irrigation = 1
    fertilization = 2
    topping = 3
    defoliation = 4


class IrrigationMethod(enum.Enum):
    sprinkler = 0
    furrow = 1
    drip = 2


class FertilizationMethod(enum.Enum):
    broadcast = 0
    sidedress = 1
    foliar = 2
    drip = 3


association_table = Table(
    "profile_agronomy_operation",
    Model.metadata,
    Column("profile_id", Integer, ForeignKey("profile.id")),
    Column("agronomy_operation_id", Integer, ForeignKey("agronomy_operation.id")),
)


class AgronomyOperation(Model):
    __tablename__ = "agronomy_operation"
    id = Column(Integer, primary_key=True)
    date = Column(Date)
    type = Column(Enum(AgronomyOperationType))
    planting_row_width = Column(Float)
    planting_skip_row_width = Column(Float)
    planting_plants_per_meter = Column(Float)
    irrigation_amount = Column(Float)
    irrigation_method = Column(Enum(IrrigationMethod))
    drip_depth = Column(Float)
    drip_horizontal_place = Column(Float)
    fertilization_ammonium = Column(Float)
    fertilization_nitrate = Column(Float)
    fertilization_urea = Column(Float)
    fertilization_method = Column(Enum(FertilizationMethod))
    defoliation_threshold = Column(Float)


class Profile(Model):
    __tablename__ = "profile"
    id = Column(Integer, primary_key=True)
    name = Column(String)
    start_date = Column(Date)
    stop_date = Column(Date)
    emerge_date = Column(Date)
    site_id = Column(Integer, ForeignKey("site.id"))
    site: "Site" = relationship("Site", backref="profiles")
    cultivar_id = Column(Integer, ForeignKey("cultivar.id"))
    cultivar: "Cultivar" = relationship("Cultivar", backref="profiles")
    agronomy_operations: list["AgronomyOperation"] = relationship(
        "AgronomyOperation", secondary="profile_agronomy_operation", backref="profiles"
    )
    simulations: list["Simulation"] = relationship("Simulation", backref="profile")


class Soil(Model):
    __tablename__ = "soil"
    id = Column(Integer, primary_key=True)
    site_id = Column(Integer, ForeignKey("site.id"))
    site: "Site" = relationship("Site", backref="soil")
    depth = Column(Float)
    ammonium_nitrogen = Column(Float)
    nitrate_nitrogen = Column(Float)
    organic_matter = Column(Float)
    water = Column(Float)


class SoilHydrology(Model):
    __tablename__ = "soil_hydrology"
    id = Column(Integer, primary_key=True)
    site_id = Column(Integer, ForeignKey("site.id"))
    site: "Site" = relationship("Site", backref="soil_hydrology")
    depth = Column(Float)
    air_dry_water_content = Column(Float)
    saturated_water_content = Column(Float)
    alpha = Column(Float)
    beta = Column(Float)
    saturated_hydraulic_conductivity = Column(Float)
    field_capacity_hydraulic_conductivity = Column(Float)
    bulk_density = Column(Float)
    clay = Column(Float)
    sand = Column(Float)


class Simulation(Model):
    __tablename__ = "simulation"
    id = Column(Integer, primary_key=True)
    profile_id = Column(Integer, ForeignKey("profile.id"))
    version = Column(Integer, default=0x0400)
    start_at = Column(DateTime)
    execute_time = Column(Float)
    states: list["State"] = relationship("State", backref="simulation")


class State(Model):
    __tablename__ = "state"
    id = Column(Integer, primary_key=True)
    simulation_id = Column(Integer, ForeignKey("simulation.id"))
    date = Column(Date)
    kday = Column(Integer)
    lint_yield = Column(Float)
    ginning_percent = Column(Float)
    plant_height = Column(Float)
    plant_weight = Column(Float)
    stem_weight = Column(Float)
    square_weight = Column(Float)
    green_bolls_weight = Column(Float)
    open_bolls_weight = Column(Float)
    leaf_weight = Column(Float)
    leaf_area_index = Column(Float)
    light_interception = Column(Float)
    number_of_squares = Column(Float)
    number_of_green_bolls = Column(Float)
    number_of_open_bolls = Column(Float)
    number_of_fruiting_branches = Column(Integer)
