from sqlalchemy import Column, Date, DateTime, Float, ForeignKey, Integer
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relationship

Model = declarative_base(name="Model")


class Simulation(Model):  # pylint: disable=too-few-public-methods
    __tablename__ = "simulations"
    id = Column(Integer, primary_key=True)
    version = Column(Integer, default=0x0400)
    execute_time = Column(DateTime)
    states: list["State"] = relationship("State", backref=backref("simulations"))


class State(Model):  # pylint: disable=too-few-public-methods
    __tablename__ = "states"
    id = Column(Integer, primary_key=True)
    simulation_id = Column(Integer, ForeignKey("simulations.id"))
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
