# distutils: language=c++
# cython: language_level=3
from datetime import datetime, date, timedelta
from math import sin, cos, acos, sqrt, pi

from libc.math cimport exp
from libc.stdlib cimport malloc
from libcpp cimport bool as bool_t
from libcpp.string cimport string

import numpy as np

from _cotton2k.climate import compute_day_length, radiation
from _cotton2k.leaf import temperature_on_leaf_growth_rate
from _cotton2k.phenology import days_to_first_square
from _cotton2k.photosynthesis import ambient_co2_factor, compute_light_interception
from _cotton2k.soil import compute_soil_surface_albedo, compute_incoming_short_wave_radiation, root_psi
from _cotton2k.utils import date2doy, doy2date
from .climate cimport ClimateStruct
from .cxx cimport (
    cSimulation,
    SandVolumeFraction,
    ClayVolumeFraction,
    DayTimeTemp,
    FoliageTemp,
    NightTimeTemp,
    LeafWaterPotential,
    StemWeight,
    AverageLwpMin,
    AverageLwp,
    LwpX,
    LwpMin,
    LwpMax,
    LwpMinX,
    NetPhotosynthesis,
    CumNetPhotosynth,
    ActualStemGrowth,
    PotGroStem,
    PotGroAllRoots,
    PotGroAllSquares,
    PotGroAllBolls,
    PotGroAllBurrs,
    PotGroAllLeaves,
    PotGroAllPetioles,
    PotGroLeafAreaPreFru,
    PotGroLeafWeightPreFru,
    PotGroPetioleWeightPreFru,
    SoilTempDailyAvrg,
    PoreSpace,
    SoilPsi,
    RootImpede,
    ReserveC,
    CultivationDepth,
    CultivationDate,
    LateralRootFlag,
    LastTaprootLayer,
    DepthLastRootLayer,
    TapRootLength,
    CarbonAllocatedForRootGrowth,
    RootNitrogen,
    RootWeightLoss,
)
from .irrigation cimport Irrigation
from .rs cimport SlabLoc, tdewest, dl, wk, daywnd, PotentialStemGrowth, AddPlantHeight, TemperatureOnFruitGrowthRate, VaporPressure, clearskyemiss, SoilNitrateOnRootGrowth, SoilAirOnRootGrowth, SoilMechanicResistance, SoilTemOnRootGrowth
from .state cimport cState, cVegetativeBranch, cFruitingBranch, cMainStemLeaf
from .fruiting_site cimport FruitingSite, Leaf


class SimulationEnd(RuntimeError):
    pass


cdef CopyState(cSimulation & sim, uint32_t i):
    cdef cState state = sim.states[i]
    state.daynum = sim.day_start + i + 1
    sim.states[i + 1] = state

cdef void initialize_switch(cSimulation & sim):
    global isw
    cdef cState state0 = sim.states[0]
    # If the date of emergence has not been given, emergence will be simulated
    # by the model. In this case, isw = 0, and a check is performed to make
    # sure that the date of planting has been given.
    if sim.day_emerge <= 0:
        if sim.day_plant <= 0:
            raise Exception(
                " planting date or emergence date must be given in the profile file !!")
        isw = 0
    # If the date of emergence has been given in the input: isw = 1 if
    # simulation starts before emergence, or isw = 2 if simulation starts at emergence.
    elif sim.day_emerge > sim.day_start:
        isw = 1
    else:
        isw = 2
        sim.states[0].kday = 1

cdef class SoilInit:
    cdef unsigned int number_of_layers
    def __init__(self, initial, hydrology, layer_depth=None):
        if layer_depth is not None:
            self.layer_depth = layer_depth
        self.initial = initial
        self.hydrology = hydrology
        self.number_of_layers = len(hydrology["layers"])

    @property
    def lyrsol(self):
        return self.number_of_layers

    @property
    def layer_depth(self):
        return LayerDepth

    @layer_depth.setter
    def layer_depth(self, value):
        global LayerDepth
        LayerDepth = value

    @property
    def initial(self):
        return [
            {
                "ammonium_nitrogen": rnnh4[i],
                "nitrate_nitrogen": rnno3[i],
                "organic_matter": oma[i],
                "water": h2oint[i]
            } for i in range(14)
        ]

    @initial.setter
    def initial(self, init_soil):
        for i, layer in enumerate(init_soil):
            rnnh4[i] = layer["ammonium_nitrogen"]
            rnno3[i] = layer["nitrate_nitrogen"]
            oma[i] = layer["organic_matter"]
            h2oint[i] = layer["water"]

    @property
    def hydrology(self):
        return {
            "ratio_implicit": RatioImplicit,
            "max_conductivity": conmax,
            "field_capacity_water_potential": psisfc,
            "immediate_drainage_water_potential": psidra,
            "layers": [
                {
                    "depth": ldepth[i],
                    "air_dry": airdr[i],
                    "theta": thetas[i],
                    "alpha": alpha[i],
                    "beta": vanGenuchtenBeta,
                    "saturated_hydraulic_conductivity": SaturatedHydCond[i],
                    "field_capacity_hydraulic_conductivity": condfc[i],
                    "bulk_density": BulkDensity[i],
                    "clay": pclay[i],
                    "sand": psand[i],
                }
                for i in range(self.number_of_layers)
            ]
        }

    @hydrology.setter
    def hydrology(self, soil_hydrology):
        global RatioImplicit, conmax, psisfc, psidra
        RatioImplicit = soil_hydrology["ratio_implicit"]
        conmax = soil_hydrology["max_conductivity"]
        psisfc = soil_hydrology["field_capacity_water_potential"]
        psidra = soil_hydrology["immediate_drainage_water_potential"]
        for i, layer in enumerate(soil_hydrology["layers"]):
            ldepth[i] = layer["depth"]
            airdr[i] = layer["air_dry"]
            thetas[i] = layer["theta"]
            alpha[i] = layer["alpha"]
            vanGenuchtenBeta[i] = layer["beta"]
            SaturatedHydCond[i] = layer["saturated_hydraulic_conductivity"]
            condfc[i] = layer["field_capacity_hydraulic_conductivity"]
            BulkDensity[i] = layer["bulk_density"]
            pclay[i] = layer["clay"]
            psand[i] = layer["sand"]

cdef class Climate:
    cdef ClimateStruct *climate
    cdef unsigned int start_day
    cdef unsigned int days
    cdef unsigned int current

    def __init__(self, start_date, climate):
        global LastDayWeatherData
        self.start_day = date2doy(start_date)
        self.current = self.start_day
        self.days = len(climate)
        self.climate = <ClimateStruct *> malloc(sizeof(ClimateStruct) * len(climate))
        LastDayWeatherData = len(climate) + self.start_day - 1
        for i, daily_climate in enumerate(climate):
            self.climate[i].Rad = daily_climate["radiation"]
            self.climate[i].Tmax = daily_climate["max"]
            self.climate[i].Tmin = daily_climate["min"]
            self.climate[i].Wind = daily_climate["wind"]
            self.climate[i].Rain = daily_climate["rain"]
            self.climate[i].Tdew = daily_climate.get("dewpoint",
                                                     tdewest(daily_climate["max"],
                                                             SitePar[5], SitePar[6]))

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start or 0
            stop = key.stop or self.days
            step = key.step or 1
            if isinstance(start, date):
                start = date2doy(start) - self.start_day
            if isinstance(stop, date):
                stop = date2doy(stop) - self.start_day
            return [{
                "radiation": self.climate[i].Rad,
                "max": self.climate[i].Tmax,
                "min": self.climate[i].Tmin,
                "wind": self.climate[i].Wind,
                "rain": self.climate[i].Rain,
                "dewpoint": self.climate[i].Tdew,
            } for i in range(start, stop, step)]
        else:
            if not isinstance(key, int):
                key = date2doy(key) - self.start_day
            climate = self.climate[key]
            return {
                "radiation": climate["Rad"],
                "max": climate["Tmax"],
                "min": climate["Tmin"],
                "wind": climate["Wind"],
                "rain": climate["Rain"],
                "dewpoint": climate["Tdew"],
            }

    def __iter__(self):
        return self

    def __next__(self):
        if self.current < self.start_day + self.days:
            self.current += 1
            return self[self.current - 1]
        else:
            raise StopIteration

cdef read_agricultural_input(cSimulation & sim, inputs):
    global NumNitApps, NumIrrigations
    NumNitApps = 0
    idef = 0
    cdef Irrigation irrigation
    cdef NitrogenFertilizer nf
    for i in inputs:
        if i["type"] == "irrigation":
            irrigation.day = date2doy(i["date"])  # day of year of this irrigation
            irrigation.amount = i["amount"]  # net amount of water applied, mm
            irrigation.method = i.get("method", 0)  # method of irrigation: 1=  2=drip
            isdhrz = i.get("drip_horizontal_place", 0)  # horizontal placement cm
            isddph = i.get("drip_depth", 0)  # vertical placement cm
            # If this is a drip irrigation, convert distances to soil
            # layer and column numbers by calling SlabLoc.
            if irrigation.method == 2:
                irrigation.LocationColumnDrip = SlabLoc(isdhrz, sim.row_space)
                irrigation.LocationLayerDrip = SlabLoc(isddph, 0)
            sim.irrigation[NumIrrigations] = irrigation
            NumIrrigations += 1
        elif i["type"] == "fertilization":
            nf.day = date2doy(i["date"])
            nf.amtamm = i.get("ammonium", 0)
            nf.amtnit = i.get("nitrate", 0)
            nf.amtura = i.get("urea", 0)
            nf.mthfrt = i.get("method", 0)
            isdhrz = i.get("drip_horizontal_place",
                           0)  # horizontal placement of DRIP, cm from left edge of soil slab.
            isddph = i.get("drip_depth",
                           0)  # vertical placement of DRIP, cm from soil surface.
            if nf.mthfrt == 1 or nf.mthfrt == 3:
                nf.ksdr = SlabLoc(isdhrz, sim.row_space)
                nf.lsdr = SlabLoc(isddph, 0)
            else:
                nf.ksdr = 0
                nf.lsdr = 0
            NFertilizer[NumNitApps] = nf
            NumNitApps += 1
        elif i["type"] == "defoliation prediction":
            DefoliationDate[idef] = date2doy(i["date"])
            DefoliantAppRate[idef] = -99.9
            if idef == 0:
                sim.day_defoliate = DefoliationDate[0]
            DefoliationMethod[idef] = i.get("method", 0)
            idef += 1


cdef class SoilCell:
    cdef cSoilCell *_
    cdef unsigned int l
    cdef unsigned int k

    @staticmethod
    cdef SoilCell from_ptr(cSoilCell *_ptr, unsigned int l, unsigned int k):
        cdef SoilCell cell = SoilCell.__new__(SoilCell)
        cell._ = _ptr
        cell.l = l
        cell.k = k
        return cell


cdef class NodeLeaf:
    cdef Leaf *_

    @property
    def area(self):
        return self._[0].area

    @area.setter
    def area(self, value):
        self._[0].area = value

    @property
    def weight(self):
        return self._[0].weight

    @weight.setter
    def weight(self, value):
        self._[0].weight = value

    @staticmethod
    cdef NodeLeaf from_ptr(Leaf *_ptr):
        cdef NodeLeaf leaf = NodeLeaf.__new__(NodeLeaf)
        leaf._ = _ptr
        return leaf


cdef class FruitingNode:
    cdef FruitingSite *_

    @property
    def average_temperature(self):
        return self._[0].average_temperature

    @average_temperature.setter
    def average_temperature(self, value):
        self._[0].average_temperature = value

    @property
    def age(self):
        return self._[0].age

    @property
    def fraction(self):
        return self._[0].fraction

    @fraction.setter
    def fraction(self, value):
        self._[0].fraction = value

    @property
    def leaf(self):
        return NodeLeaf.from_ptr(&self._.leaf)

    @property
    def stage(self):
        return self._[0].stage

    @stage.setter
    def stage(self, value):
        self._[0].stage = value

    @staticmethod
    cdef FruitingNode from_ptr(FruitingSite *_ptr):
        cdef FruitingNode node = FruitingNode.__new__(FruitingNode)
        node._ = _ptr
        return node


cdef class MainStemLeaf:
    cdef cMainStemLeaf *_

    @property
    def area(self):
        return self._.leaf_area

    @area.setter
    def area(self, value):
        self._.leaf_area = value

    @property
    def weight(self):
        return self._.leaf_weight

    @weight.setter
    def weight(self, value):
        self._.leaf_weight = value

    @staticmethod
    cdef MainStemLeaf from_ptr(cMainStemLeaf *_ptr):
        """Factory function to create WrapperClass objects from
        given my_c_struct pointer.

        Setting ``owner`` flag to ``True`` causes
        the extension type to ``free`` the structure pointed to by ``_ptr``
        when the wrapper object is deallocated."""
        # Call to __new__ bypasses __init__ constructor
        cdef MainStemLeaf main_stem_leaf = MainStemLeaf.__new__(MainStemLeaf)
        main_stem_leaf._ = _ptr
        return main_stem_leaf


cdef class FruitingBranch:
    cdef cFruitingBranch *_

    @property
    def number_of_fruiting_nodes(self):
        return self._[0].number_of_fruiting_nodes

    @number_of_fruiting_nodes.setter
    def number_of_fruiting_nodes(self, value):
        self._[0].number_of_fruiting_nodes = value

    @property
    def delay_for_new_node(self):
        return self._[0].delay_for_new_node

    @delay_for_new_node.setter
    def delay_for_new_node(self, value):
        self._[0].delay_for_new_node = value

    @property
    def main_stem_leaf(self):
        return MainStemLeaf.from_ptr(&self._[0].main_stem_leaf)

    @property
    def nodes(self):
        return [FruitingNode.from_ptr(&self._.nodes[i]) for i in
                range(self._.number_of_fruiting_nodes)]

    @staticmethod
    cdef FruitingBranch from_ptr(cFruitingBranch *_ptr):
        """Factory function to create WrapperClass objects from
        given my_c_struct pointer.

        Setting ``owner`` flag to ``True`` causes
        the extension type to ``free`` the structure pointed to by ``_ptr``
        when the wrapper object is deallocated."""
        # Call to __new__ bypasses __init__ constructor
        cdef FruitingBranch fruiting_branch = FruitingBranch.__new__(FruitingBranch)
        fruiting_branch._ = _ptr
        return fruiting_branch


cdef class VegetativeBranch:
    cdef cVegetativeBranch *_

    @property
    def number_of_fruiting_branches(self):
        return self._[0].number_of_fruiting_branches

    @number_of_fruiting_branches.setter
    def number_of_fruiting_branches(self, value):
        self._[0].number_of_fruiting_branches = value

    @property
    def delay_for_new_fruiting_branch(self):
        return self._[0].delay_for_new_fruiting_branch

    @delay_for_new_fruiting_branch.setter
    def delay_for_new_fruiting_branch(self, value):
        self._[0].delay_for_new_fruiting_branch = value

    @property
    def fruiting_branches(self):
        return [FruitingBranch.from_ptr(&self._[0].fruiting_branches[i]) for i in
                range(self._[0].number_of_fruiting_branches)]

    @staticmethod
    cdef VegetativeBranch from_ptr(cVegetativeBranch *_ptr):
        """Factory function to create WrapperClass objects from
        given my_c_struct pointer.

        Setting ``owner`` flag to ``True`` causes
        the extension type to ``free`` the structure pointed to by ``_ptr``
        when the wrapper object is deallocated."""
        # Call to __new__ bypasses __init__ constructor
        cdef VegetativeBranch vegetative_branch = VegetativeBranch.__new__(VegetativeBranch)
        vegetative_branch._ = _ptr
        return vegetative_branch


cdef class Hour:
    cdef cHour *_

    @staticmethod
    cdef Hour from_ptr(cHour *_ptr):
        """Factory function to create WrapperClass objects from
        given my_c_struct pointer.

        Setting ``owner`` flag to ``True`` causes
        the extension type to ``free`` the structure pointed to by ``_ptr``
        when the wrapper object is deallocated."""
        # Call to __new__ bypasses __init__ constructor
        cdef Hour hour = Hour.__new__(Hour)
        hour._ = _ptr
        return hour

cdef double[3] cgind = [1, 1, 0.10]  # the index for the capability of growth of class I roots (0 to 1).


cdef class Soil:
    cdef cSoil *_

    @property
    def number_of_layers_with_root(self):
        return self._[0].number_of_layers_with_root

    @number_of_layers_with_root.setter
    def number_of_layers_with_root(self, value):
        self._[0].number_of_layers_with_root = value

    @staticmethod
    cdef Soil from_ptr(cSoil *_ptr):
        cdef Soil soil = Soil.__new__(Soil)
        soil._ = _ptr
        return soil


cdef class State:
    cdef cState *_
    cdef public unsigned int year
    cdef public unsigned int version
    cdef double pavail  # residual available carbon for root growth from previous day.

    def keys(self):
        return dict(self._[0]).keys()

    @property
    def daynum(self):
        return self._[0].daynum

    @daynum.setter
    def daynum(self, value):
        self._[0].daynum = value

    @property
    def date(self):
        return doy2date(self.year, self._[0].daynum)

    @date.setter
    def date(self, value):
        if not isinstance(value, date):
            raise TypeError
        self._[0].daynum = date2doy(value)

    @property
    def day_inc(self):
        return self._[0].day_inc

    @day_inc.setter
    def day_inc(self, value):
        self._[0].day_inc = value

    @property
    def kday(self):
        return self._[0].kday

    @kday.setter
    def kday(self, value):
        self._[0].kday = value

    @property
    def average_temperature(self):
        return self._[0].average_temperature

    @average_temperature.setter
    def average_temperature(self, value):
        self._[0].average_temperature = value

    @property
    def lint_yield(self):
        return self._[0].lint_yield

    @lint_yield.setter
    def lint_yield(self, value):
        self._[0].lint_yield = value

    @property
    def plant_height(self):
        return self._[0].plant_height

    @plant_height.setter
    def plant_height(self, value):
        self._[0].plant_height = value

    @property
    def plant_weight(self):
        return self._[0].plant_weight

    @plant_weight.setter
    def plant_weight(self, value):
        self._[0].plant_weight = value

    @property
    def stem_weight(self):
        return self._[0].stem_weight

    @stem_weight.setter
    def stem_weight(self, value):
        self._[0].stem_weight = value

    @property
    def petiole_weight(self):
        return self._[0].petiole_weight

    @petiole_weight.setter
    def petiole_weight(self, value):
        self._[0].petiole_weight = value

    @property
    def square_weight(self):
        return self._[0].square_weight

    @square_weight.setter
    def square_weight(self, value):
        self._[0].square_weight = value

    @property
    def green_bolls_weight(self):
        return self._[0].green_bolls_weight

    @green_bolls_weight.setter
    def green_bolls_weight(self, value):
        self._[0].green_bolls_weight = value

    @property
    def green_bolls_burr_weight(self):
        return self._[0].green_bolls_burr_weight

    @green_bolls_burr_weight.setter
    def green_bolls_burr_weight(self, value):
        self._[0].green_bolls_burr_weight = value

    @property
    def open_bolls_weight(self):
        return self._[0].open_bolls_weight

    @open_bolls_weight.setter
    def open_bolls_weight(self, value):
        self._[0].open_bolls_weight = value

    @property
    def open_bolls_burr_weight(self):
        return self._[0].open_bolls_burr_weight

    @open_bolls_burr_weight.setter
    def open_bolls_burr_weight(self, value):
        self._[0].open_bolls_burr_weight = value

    @property
    def bloom_weight_loss(self):
        return self._[0].bloom_weight_loss

    @bloom_weight_loss.setter
    def bloom_weight_loss(self, value):
        self._[0].bloom_weight_loss = value

    @property
    def abscised_fruit_sites(self):
        return self._[0].abscised_fruit_sites

    @abscised_fruit_sites.setter
    def abscised_fruit_sites(self, value):
        self._[0].abscised_fruit_sites = value

    @property
    def abscised_leaf_weight(self):
        return self._[0].abscised_leaf_weight

    @abscised_leaf_weight.setter
    def abscised_leaf_weight(self, value):
        self._[0].abscised_leaf_weight = value

    @property
    def cumulative_nitrogen_loss(self):
        return self._[0].cumulative_nitrogen_loss

    @cumulative_nitrogen_loss.setter
    def cumulative_nitrogen_loss(self, value):
        self._[0].cumulative_nitrogen_loss = value

    @property
    def cumulative_transpiration(self):
        return self._[0].cumulative_transpiration

    @cumulative_transpiration.setter
    def cumulative_transpiration(self, value):
        self._[0].cumulative_transpiration = value

    @property
    def cumulative_evaporation(self):
        return self._[0].cumulative_evaporation

    @cumulative_evaporation.setter
    def cumulative_evaporation(self, value):
        self._[0].cumulative_evaporation = value

    @property
    def applied_water(self):
        return self._[0].applied_water

    @applied_water.setter
    def applied_water(self, value):
        self._[0].applied_water = value

    @property
    def water_stress(self):
        return self._[0].water_stress

    @water_stress.setter
    def water_stress(self, value):
        self._[0].water_stress = value

    @property
    def water_stress_stem(self):
        return self._[0].water_stress_stem

    @water_stress_stem.setter
    def water_stress_stem(self, value):
        self._[0].water_stress_stem = value

    @property
    def carbon_stress(self):
        return self._[0].carbon_stress

    @carbon_stress.setter
    def carbon_stress(self, value):
        self._[0].carbon_stress = value

    @property
    def extra_carbon(self):
        return self._[0].extra_carbon

    @extra_carbon.setter
    def extra_carbon(self, value):
        self._[0].extra_carbon = value

    @property
    def light_interception(self):
        return self._[0].light_interception

    @light_interception.setter
    def light_interception(self, value):
        self._[0].light_interception = value

    @property
    def leaf_area_index(self):
        return self._[0].leaf_area_index

    @leaf_area_index.setter
    def leaf_area_index(self, value):
        self._[0].leaf_area_index = value

    @property
    def leaf_area(self):
        return self._[0].leaf_area

    @leaf_area.setter
    def leaf_area(self, value):
        self._[0].leaf_area = value

    @property
    def leaf_weight(self):
        return self._[0].leaf_weight

    @leaf_weight.setter
    def leaf_weight(self, value):
        self._[0].leaf_weight = value

    @property
    def leaf_weight_area_ratio(self):
        return self._[0].leaf_weight_area_ratio

    @leaf_weight_area_ratio.setter
    def leaf_weight_area_ratio(self, value):
        self._[0].leaf_weight_area_ratio = value

    @property
    def leaf_nitrogen(self):
        return self._[0].leaf_nitrogen

    @leaf_nitrogen.setter
    def leaf_nitrogen(self, value):
        self._[0].leaf_nitrogen = value

    @property
    def leaf_nitrogen_concentration(self):
        return self._[0].leaf_nitrogen_concentration

    @leaf_nitrogen_concentration.setter
    def leaf_nitrogen_concentration(self, value):
        self._[0].leaf_nitrogen_concentration = value

    @property
    def number_of_vegetative_branches(self):
        return self._[0].number_of_vegetative_branches

    @number_of_vegetative_branches.setter
    def number_of_vegetative_branches(self, value):
        self._[0].number_of_vegetative_branches = value

    @property
    def number_of_squares(self):
        return self._[0].number_of_squares

    @number_of_squares.setter
    def number_of_squares(self, value):
        self._[0].number_of_squares = value

    @property
    def number_of_green_bolls(self):
        return self._[0].number_of_green_bolls

    @number_of_green_bolls.setter
    def number_of_green_bolls(self, value):
        self._[0].number_of_green_bolls = value

    @property
    def number_of_open_bolls(self):
        return self._[0].number_of_open_bolls

    @number_of_open_bolls.setter
    def number_of_open_bolls(self, value):
        self._[0].number_of_open_bolls = value

    @property
    def nitrogen_stress(self):
        return self._[0].nitrogen_stress

    @nitrogen_stress.setter
    def nitrogen_stress(self, value):
        self._[0].nitrogen_stress = value

    @property
    def nitrogen_stress_vegetative(self):
        return self._[0].nitrogen_stress_vegetative

    @nitrogen_stress_vegetative.setter
    def nitrogen_stress_vegetative(self, value):
        self._[0].nitrogen_stress_vegetative = value

    @property
    def nitrogen_stress_fruiting(self):
        return self._[0].nitrogen_stress_fruiting

    @nitrogen_stress_fruiting.setter
    def nitrogen_stress_fruiting(self, value):
        self._[0].nitrogen_stress_fruiting = value

    @property
    def nitrogen_stress_root(self):
        return self._[0].nitrogen_stress_root

    @nitrogen_stress_root.setter
    def nitrogen_stress_root(self, value):
        self._[0].nitrogen_stress_root = value

    @property
    def total_required_nitrogen(self):
        return self._[0].total_required_nitrogen

    @total_required_nitrogen.setter
    def total_required_nitrogen(self, value):
        self._[0].total_required_nitrogen = value

    @property
    def petiole_nitrogen_concentration(self):
        return self._[0].petiole_nitrogen_concentration

    @petiole_nitrogen_concentration.setter
    def petiole_nitrogen_concentration(self, value):
        self._[0].petiole_nitrogen_concentration = value

    @property
    def seed_nitrogen(self):
        return self._[0].seed_nitrogen

    @seed_nitrogen.setter
    def seed_nitrogen(self, value):
        self._[0].seed_nitrogen = value

    @property
    def seed_nitrogen_concentration(self):
        return self._[0].seed_nitrogen_concentration

    @seed_nitrogen_concentration.setter
    def seed_nitrogen_concentration(self, value):
        self._[0].seed_nitrogen_concentration = value

    @property
    def burr_nitrogen(self):
        return self._[0].burr_nitrogen

    @burr_nitrogen.setter
    def burr_nitrogen(self, value):
        self._[0].burr_nitrogen = value

    @property
    def root_nitrogen_concentration(self):
        return self._[0].root_nitrogen_concentration

    @root_nitrogen_concentration.setter
    def root_nitrogen_concentration(self, value):
        self._[0].root_nitrogen_concentration = value

    @property
    def square_nitrogen(self):
        return self._[0].square_nitrogen

    @square_nitrogen.setter
    def square_nitrogen(self, value):
        self._[0].square_nitrogen = value

    @property
    def square_nitrogen_concentration(self):
        return self._[0].square_nitrogen_concentration

    @square_nitrogen_concentration.setter
    def square_nitrogen_concentration(self, value):
        self._[0].square_nitrogen_concentration = value

    @property
    def stem_nitrogen(self):
        return self._[0].stem_nitrogen

    @stem_nitrogen.setter
    def stem_nitrogen(self, value):
        self._[0].stem_nitrogen = value

    @property
    def fruit_growth_ratio(self):
        return self._[0].fruit_growth_ratio

    @fruit_growth_ratio.setter
    def fruit_growth_ratio(self, value):
        self._[0].fruit_growth_ratio = value

    @property
    def ginning_percent(self):
        return self._[0].ginning_percent

    @ginning_percent.setter
    def ginning_percent(self, value):
        self._[0].ginning_percent = value

    @property
    def number_of_pre_fruiting_nodes(self):
        return self._[0].number_of_pre_fruiting_nodes

    @number_of_pre_fruiting_nodes.setter
    def number_of_pre_fruiting_nodes(self, value):
        self._[0].number_of_pre_fruiting_nodes = value

    @property
    def age_of_pre_fruiting_nodes(self):
        return self._[0].age_of_pre_fruiting_nodes

    @property
    def leaf_area_pre_fruiting(self):
        return self._[0].leaf_area_pre_fruiting

    @property
    def leaf_weight_pre_fruiting(self):
        return self._[0].leaf_weight_pre_fruiting

    @property
    def soil(self):
        return Soil.from_ptr(&self._[0].soil)

    def __getitem__(self, item):
        return getattr(self, item)

    @staticmethod
    cdef State from_ptr(cState *_ptr, unsigned int year, unsigned int version):
        """Factory function to create WrapperClass objects from
        given my_c_struct pointer.

        Setting ``owner`` flag to ``True`` causes
        the extension type to ``free`` the structure pointed to by ``_ptr``
        when the wrapper object is deallocated."""
        # Call to __new__ bypasses __init__ constructor
        cdef State state = State.__new__(State)
        state._ = _ptr
        state.year = year
        state.version = version
        return state

    @property
    def vegetative_branches(self):
        return [VegetativeBranch.from_ptr(&self._[0].vegetative_branches[k]) for k in range(self.number_of_vegetative_branches)]

    @property
    def hours(self):
        return (Hour.from_ptr(&self._[0].hours[i]) for i in range(24))

    def create_first_square(self, stemNRatio, first_square_leaf_area):
        """
        This function initiates the first square. It is called from function CottonPhenology().
        """
        # FruitFraction and FruitingCode are assigned 1 for the first fruiting site.
        self.vegetative_branches[0].number_of_fruiting_branches = 1
        self.vegetative_branches[0].fruiting_branches[0].number_of_fruiting_nodes = 1
        first_node = self.vegetative_branches[0].fruiting_branches[0].nodes[0]
        first_node.stage = Stage.Square
        first_node.fraction = 1
        # Initialize a new leaf at this position. define its initial weight and area. VarPar[34] is the initial area of a new leaf. The mass and nitrogen of the new leaf are substacted from the stem.
        first_node.leaf.area = first_square_leaf_area
        first_node.leaf.weight = first_square_leaf_area * self.leaf_weight_area_ratio
        self.stem_weight -= first_node.leaf.weight
        self.leaf_weight += first_node.leaf.weight
        self.leaf_nitrogen += first_node.leaf.weight * stemNRatio
        self.stem_nitrogen -= first_node.leaf.weight * stemNRatio
        first_node.average_temperature = self.average_temperature
        # Define the initial values of NumFruitBranches, NumNodes, state.fruit_growth_ratio, and AvrgNodeTemper.
        self.fruit_growth_ratio = 1
        # It is assumed that the cotyledons are dropped at time of first square. compute changes in AbscisedLeafWeight, state.leaf_weight, state.leaf_nitrogen and CumPlantNLoss caused by the abscission of the cotyledons.
        cdef double cotylwt = 0.20  # cotylwt is the leaf weight of the cotyledons.
        self.abscised_leaf_weight += cotylwt
        self.leaf_weight -= cotylwt
        self.cumulative_nitrogen_loss += cotylwt * self.leaf_nitrogen / self.leaf_weight
        self.leaf_nitrogen -= cotylwt * self.leaf_nitrogen / self.leaf_weight

    def add_fruiting_branch(self, k, density_factor, delayVegByCStress, delayFrtByCStress, stemNRatio, PhenDelayByNStress, time_to_new_fruiting_branch, new_node_initial_leaf_area, topping_date=None):
        """
        This function decides if a new fruiting branch is to be added to a vegetative branch, and forms it. It is called from function CottonPhenology().
        """
        if topping_date is not None and self.date >= topping_date:
            return
        # The following constant parameters are used:
        cdef double[8] vfrtbr = [0.8, 0.95, 33.0, 4.461, -0.1912, 0.00265, 1.8, -1.32]
        # Compute the cumulative delay for the appearance of the next caused by carbohydrate, nitrogen, and water stresses.
        vegetative_branch = self.vegetative_branches[k]
        vegetative_branch.delay_for_new_fruiting_branch += delayVegByCStress + vfrtbr[0] * PhenDelayByNStress
        vegetative_branch.delay_for_new_fruiting_branch += vfrtbr[1] * (1 - self.water_stress)
        # Define nbrch and compute TimeToNextFruBranch, the time (in physiological days) needed for the formation of each successive fruiting branch, as a function of the average temperature. This function is derived from data of K. R. Reddy, CSRU, adjusted for age expressed in physiological days.
        # It is different for the main stem (k = 0) than for the other vegetative branches. TimeToNextFruNode is modified for plant density. Add DelayNewFruBranch to TimeToNextFruNode.
        last_fruiting_branch = vegetative_branch.fruiting_branches[-1]
        cdef double tav = last_fruiting_branch.nodes[0].average_temperature  # modified average daily temperature.
        if tav > vfrtbr[2]:
            tav = vfrtbr[2]
        # TimeToNextFruBranch is the time, in physiological days, for the next fruiting branch to be formed.
        cdef double TimeToNextFruBranch = time_to_new_fruiting_branch + tav * (vfrtbr[3] + tav * (vfrtbr[4] + tav * vfrtbr[5]))
        if k > 0:
            TimeToNextFruBranch = TimeToNextFruBranch * vfrtbr[6]
        TimeToNextFruBranch = TimeToNextFruBranch * (1 + vfrtbr[7] * (1 - density_factor)) + vegetative_branch.delay_for_new_fruiting_branch
        # Check if the the age of the last fruiting branch exceeds TimeToNextFruBranch. If so, form the new fruiting branch:
        if last_fruiting_branch.nodes[0].age < TimeToNextFruBranch:
            return
        # Increment NumFruitBranches, define newbr, and assign 1 to NumNodes, FruitFraction and FruitingCode.
        vegetative_branch.number_of_fruiting_branches += 1
        if vegetative_branch.number_of_fruiting_branches > 30:
            vegetative_branch.number_of_fruiting_branches = 30
            return
        cdef int newbr  # the index number of the new fruiting branch on this vegetative branch, after a new branch has been added.
        newbr = vegetative_branch.number_of_fruiting_branches - 1
        new_branch = vegetative_branch.fruiting_branches[-1]
        new_branch.number_of_fruiting_nodes = 1
        new_node = new_branch.nodes[0]
        new_node.fraction = 1
        new_node.stage = Stage.Square
        # Initiate new leaves at the first node of the new fruiting branch, and at the corresponding main stem node. The mass and nitrogen in the new leaves is substacted from the stem.
        new_node.leaf.area = new_node_initial_leaf_area
        new_node.leaf.weight = new_node_initial_leaf_area * self.leaf_weight_area_ratio
        main_stem_leaf = new_branch.main_stem_leaf
        main_stem_leaf.area = new_node_initial_leaf_area
        main_stem_leaf.weight = main_stem_leaf.area * self.leaf_weight_area_ratio
        self.stem_weight -= main_stem_leaf.weight + new_node.leaf.weight
        self.leaf_weight += main_stem_leaf.weight + new_node.leaf.weight
        # addlfn is the nitrogen added to new leaves from stem.
        cdef double addlfn = (main_stem_leaf.weight + new_node.leaf.weight) * stemNRatio
        self.leaf_nitrogen += addlfn
        self.stem_nitrogen -= addlfn
        # Begin computing AvrgNodeTemper of the new node and assign zero to DelayNewFruBranch.
        new_node.average_temperature = self.average_temperature
        vegetative_branch.delay_for_new_fruiting_branch = 0

    def potential_root_growth(self, NumRootAgeGroups, NumLayersWithRoots, per_plant_area):
        """
        This function calculates the potential root growth rate.
        The return value is the sum of potential root growth rates for the whole slab (sumpdr).
        It is called from PlantGrowth().
        It calls: RootImpedance(), SoilNitrateOnRootGrowth(), SoilAirOnRootGrowth(), SoilMechanicResistance(), SoilTemOnRootGrowth() and root_psi().
        """
        # The following constant parameter is used:
        cdef double rgfac = 0.36 if self.version < 0x500 else 0.2  # potential relative growth rate of the roots (g/g/day).
        # Initialize to zero the PotGroRoots array.
        for l in range(NumLayersWithRoots):
            for k in range(nk):
                self._[0].soil.cells[l][k].root.potential_growth = 0
        RootImpedance(self._[0].soil.cells)
        cdef double sumpdr = 0  # sum of potential root growth rate for the whole slab
        for l in range(NumLayersWithRoots):
            for k in range(nk):
                # Check if this soil cell contains roots (if RootAge is greater than 0), and execute the following if this is true.
                # In each soil cell with roots, the root weight capable of growth rtwtcg is computed as the sum of RootWeight[l][k][i] * cgind[i] for all root classes.
                if self._[0].soil.cells[l][k].root.age > 0:
                    rtwtcg = 0  # root weight capable of growth in a soil soil cell.
                    for i in range(NumRootAgeGroups):
                        rtwtcg += self._[0].soil.cells[l][k].root.weight[i] * cgind[i]
                    # Compute the temperature factor for root growth by calling function SoilTemOnRootGrowth() for this layer.
                    stday = SoilTempDailyAvrg[l][k] - 273.161  # soil temperature, C, this day's average for this cell.
                    temprg = SoilTemOnRootGrowth(stday)  # effect of soil temperature on root growth.
                    # Compute soil mechanical resistance for each soil cell by calling SoilMechanicResistance{}, the effect of soil aeration on root growth by calling SoilAirOnRootGrowth(), and the effect of soil nitrate on root growth by calling SoilNitrateOnRootGrowth().
                    rtpct # effect of soil mechanical resistance on root growth (returned from SoilMechanicResistance).

                    lp1 = l if l == nl - 1 else l + 1  # layer below l.

                    # columns to the left and to the right of k.
                    kp1 = min(k + 1, nk)
                    km1 = max(k - 1, 0)

                    rtimpd0 = RootImpede[l][k]
                    rtimpdkm1 = RootImpede[l][km1]
                    rtimpdkp1 = RootImpede[l][kp1]
                    rtimpdlp1 = RootImpede[lp1][k]
                    rtimpdmin = min(rtimpd0, rtimpdkm1, rtimpdkp1, rtimpdlp1)  # minimum value of rtimpd
                    rtpct = SoilMechanicResistance(rtimpdmin)
                    # effect of oxygen deficiency on root growth (returned from SoilAirOnRootGrowth).
                    rtrdo = SoilAirOnRootGrowth(SoilPsi[l][k], PoreSpace[l], self._[0].soil.cells[l][k].water_content)
                    # effect of nitrate deficiency on root growth (returned from SoilNitrateOnRootGrowth).
                    rtrdn = SoilNitrateOnRootGrowth(self._[0].soil.cells[l][k].nitrate_nitrogen_content)
                    # The root growth resistance factor RootGroFactor(l,k), which can take a value between 0 and 1, is computed as the minimum of these resistance factors. It is further modified by multiplying it by the soil moisture function root_psi().
                    # Potential root growth PotGroRoots(l,k) in each cell is computed as a product of rtwtcg, rgfac, the temperature function temprg, and RootGroFactor(l,k). It is also multiplied by per_plant_area / 19.6, for the effect of plant population density on root growth: it is made comparable to a population of 5 plants per m in 38" rows.
                    # The sum of the potential growth for the whole slab is computed as sumpdr.
                    self._[0].soil.cells[l][k].root.growth_factor = root_psi(SoilPsi[l][k]) * min(rtrdo, rtpct, rtrdn)
                    self._[0].soil.cells[l][k].root.potential_growth = rtwtcg * rgfac * temprg * self._[0].soil.cells[l][k].root.growth_factor * per_plant_area / 19.6
                    sumpdr += self._[0].soil.cells[l][k].root.potential_growth
        return sumpdr

    def compute_actual_root_growth(self, double sumpdr, double row_space, double per_plant_area, int NumRootAgeGroups, unsigned int day_emerge, unsigned int plant_row_column):
        # The following constant parameters are used:
        # The index for the relative partitioning of root mass produced by new growth to class i.
        global RootNitrogen, RootWeightLoss
        cdef double[3] RootGrowthIndex = [1.0, 0.0, 0.0]
        cdef double rtminc = 0.0000001  # the threshold ratio of root mass capable of growth
        # to soil cell volume (g/cm3); when this threshold is reached, a part of root growth in this cell may be extended to adjoining cells.
        # Assign zero to pavail if this is the day of emergence.
        if self.daynum <= day_emerge:
            self.pavail = 0
        adwr1 = np.zeros((40, 20), dtype=np.float64)  # actual growth rate from roots existing in this soil cell.
        # Assign zero to the arrays of actual root growth rate.
        for l in range(nl):
            for k in range(nk):
                self._[0].soil.cells[l][k].root.actual_growth = 0
        # The amount of carbon allocated for root growth is calculated from CarbonAllocatedForRootGrowth, converted to g dry matter per slab, and added to previously allocated carbon that has not been used for growth (pavail). if there is no potential root growth, this will be stored in pavail. Otherwise, zero is assigned to pavail.
        if sumpdr <= 0:
            self.pavail += CarbonAllocatedForRootGrowth * 0.01 * row_space / per_plant_area
            return
        cdef double actgf  # actual growth factor (ratio of available C to potential growth).
        # The ratio of available C to potential root growth (actgf) is calculated. pavail (if not zero) is used here, and zeroed after being used.
        actgf = (self.pavail + CarbonAllocatedForRootGrowth * 0.01 * row_space / per_plant_area) / sumpdr
        self.pavail = 0

        for l in range(self._[0].soil.number_of_layers_with_root):
            for k in range(nk):
                # adwr1(l,k), is proportional to the potential growth rate of roots in this cell.
                if self._[0].soil.cells[l][k].root.age > 0:
                    adwr1[l][k] = self._[0].soil.cells[l][k].root.potential_growth * actgf
        # If extra carbon is available, it is assumed to be added to the taproot.
        if self.extra_carbon > 0:
            # available carbon for taproot growth, in g dry matter per slab.
            # ExtraCarbon is converted to availt (g dry matter per slab).
            availt = self.extra_carbon * 0.01 * row_space / per_plant_area
            sdl = TapRootLength - DepthLastRootLayer
            # distance from the tip of the taproot, cm.
            tpwt = np.zeros((40, 2))  # proportionality factors for allocating added dry matter among taproot soil cells.
            sumwt = 0  # sum of the tpwt array.
            # Extra Carbon (availt) is added to soil cells with roots in the columns immediately to the left and to the right of the location of the plant row.
            for l in reversed(range(LastTaprootLayer + 1)):
                # The weighting factors for allocating the carbon (tpwt) are proportional to the volume of each soil cell and its distance (sdl) from the tip of the taproot.
                sdl += dl(l)
                tpwt[l][0] = sdl * dl(l) * wk(plant_row_column, row_space)
                tpwt[l][1] = sdl * dl(l) * wk(plant_row_column + 1, row_space)
                sumwt += tpwt[l][0] + tpwt[l][1]
            # The proportional amount of mass is added to the mass of the last (inactive) root class in each soil cell.
            for l in range(LastTaprootLayer + 1):
                self._[0].soil.cells[l][plant_row_column].root.weight[NumRootAgeGroups - 1] += availt * tpwt[l][0] / sumwt
                self._[0].soil.cells[l][plant_row_column + 1].root.weight[NumRootAgeGroups - 1] += availt * tpwt[l][1] / sumwt
        # Check each cell if the ratio of root weight capable of growth to cell volume (rtconc) exceeds the threshold rtminc, and call RedistRootNewGrowth() for this cell.
        # Otherwise, all new growth is contained in the same cell, and the actual growth in this cell, ActualRootGrowth(l,k) will be equal to adwr1(l,k).
        for l in range(nl):
            for k in range(nk):
                if self._[0].soil.cells[l][k].root.age > 0:
                    rtconc = 0  # ratio of root weight capable of growth to cell volume.
                    for i in range(NumRootAgeGroups):
                        rtconc += self._[0].soil.cells[l][k].root.weight[i] * cgind[i]
                    rtconc = rtconc / (dl(l) * wk(k, row_space))
                    if rtconc > rtminc:
                        RedistRootNewGrowth(self._[0], l, k, adwr1[l][k], row_space, plant_row_column)
                    else:
                        self._[0].soil.cells[l][k].root.actual_growth += adwr1[l][k]
        # The new actual growth ActualRootGrowth(l,k) in each cell is partitioned among the root classes in it in proportion to the parameters RootGrowthIndex(i), and the previous values of RootWeight(k,l,i), and added to RootWeight(k,l,i).
        sumind = 0  # sum of the growth index for all classes in a cell.
        for i in range(NumRootAgeGroups):
            sumind += RootGrowthIndex[i]

        for l in range(self._[0].soil.number_of_layers_with_root):
            for k in range(nk):
                if self._[0].soil.cells[l][k].root.age > 0:
                    sumgr = 0  # sum of growth index multiplied by root weight, for all classes in a cell.
                    for i in range(NumRootAgeGroups):
                        sumgr += RootGrowthIndex[i] * self._[0].soil.cells[l][k].root.weight[i]
                    for i in range(NumRootAgeGroups):
                        if sumgr > 0:
                            self._[0].soil.cells[l][k].root.weight[i] += self._[0].soil.cells[l][k].root.actual_growth * RootGrowthIndex[i] * self._[0].soil.cells[l][k].root.weight[i] / sumgr
                        else:
                            self._[0].soil.cells[l][k].root.weight[i] += self._[0].soil.cells[l][k].root.actual_growth * RootGrowthIndex[i] / sumind
        # Call function TapRootGrowth() for taproot elongation, if the taproot has not already reached the bottom of the slab.
        if LastTaprootLayer < nl - 1 or TapRootLength < DepthLastRootLayer:
            TapRootGrowth(self._[0], NumRootAgeGroups, plant_row_column)
        # Call functions for growth of lateral roots
        InitiateLateralRoots()
        for l in range(LastTaprootLayer):
            if LateralRootFlag[l] == 2:
                LateralRootGrowthLeft(self._[0], l, NumRootAgeGroups, plant_row_column, row_space)
                LateralRootGrowthRight(self._[0], l, NumRootAgeGroups, plant_row_column, row_space)
        # Initialize DailyRootLoss (weight of sloughed roots) for this day.
        cdef double DailyRootLoss = 0  # total weight of sloughed roots, g per plant per day.
        for l in range(self._[0].soil.number_of_layers_with_root):
            for k in range(nk):
                # Check RootAge to determine if this soil cell contains roots, and then compute root aging and root death by calling RootAging() and RootDeath() for each soil cell with roots.
                if self._[0].soil.cells[l][k].root.age > 0:
                    RootAging(self._[0].soil.cells[l][k], l, k)
                    DailyRootLoss = RootDeath(self._[0].soil.cells[l][k], l, k, DailyRootLoss)
        # Check if cultivation is executed in this day and call RootCultivation().
        for j in range(5):
            if CultivationDate[j] == self.daynum:
                DailyRootLoss = RootCultivation(self._[0].soil.cells, NumRootAgeGroups, CultivationDepth[j], DailyRootLoss, row_space)
        # Convert DailyRootLoss to g per plant units and add it to RootWeightLoss.
        DailyRootLoss = DailyRootLoss * 100. * per_plant_area / row_space
        RootWeightLoss += DailyRootLoss
        # Adjust RootNitrogen (root N content) for loss by death of roots.
        RootNitrogen -= DailyRootLoss * self.root_nitrogen_concentration
        self.cumulative_nitrogen_loss += DailyRootLoss * self.root_nitrogen_concentration
        # Call function RootSummation().
        RootSummation(self._[0], NumRootAgeGroups, row_space, per_plant_area)

    def leaf_abscission(self, per_plant_area, first_square_date, defoliate_date):
        global ReserveC
        # If there are almost no leaves, this routine is not executed.
        if self.leaf_area_index <= 0.0001:
            return
        # Compute droplf as a function of LeafAreaIndex.
        p0 = 140 if self.version < 0x500 else 100
        p1 = -1
        droplf = p0 + p1 * self.leaf_area_index  # leaf age until its abscission.
        # Call PreFruitLeafAbscission() to simulate the physiological abscission of prefruiting node leaves.
        PreFruitLeafAbscission(self._[0], droplf, self.daynum, date2doy(first_square_date), date2doy(defoliate_date), self.day_inc)
        # Loop for all vegetative branches and fruiting branches, and call MainStemLeafAbscission() for each fruiting branch to simulate the physiological abscission of the other leaves.
        for k in range(self.number_of_vegetative_branches):
            for l in range(self.vegetative_branches[k].number_of_fruiting_branches):
                MainStemLeafAbscission(self._[0], k, l, droplf, self.daynum, date2doy(defoliate_date))
        # Call DefoliationLeafAbscission() to simulate leaf abscission caused by defoliants.
        if date2doy(defoliate_date) > 0 and self.daynum >= date2doy(defoliate_date):
            DefoliationLeafAbscission(self._[0], date2doy(defoliate_date))
        # If the reserves in the leaf are too high, add the lost reserves to AbscisedLeafWeight and adjust ReserveC.
        if ReserveC > 0:
            # maximum possible amount of reserve C in the leaves.
            resmax = 0.2 * self.leaf_weight
            if ReserveC > resmax:
                self.abscised_leaf_weight += ReserveC - resmax
                ReserveC = resmax
        # Compute the resulting LeafAreaIndex but do not let it get too small.
        self.leaf_area_index = max(0.0001, self.leaf_area / per_plant_area)


def compute_incoming_long_wave_radiation(humidity: float, temperature: float, cloud_cov: float, cloud_cor: float) -> float:
    """LONG WAVE RADIATION EMITTED FROM SKY"""
    stefa1 = 1.38e-12  # Stefan-Boltsman constant.
    vp = 0.01 * humidity * VaporPressure(temperature)  # air vapor pressure, KPa.
    ea0 = clearskyemiss(vp, temperature + 273.161)  # sky emissivity from clear portions of the sky.
    # incoming long wave radiation (ly / sec).
    rlzero = (ea0 * (1 - cloud_cov) + cloud_cov) * stefa1 * (temperature + 273.161) ** 4 - cloud_cor / 41880  # CloudTypeCorr converted from W m-2 to ly sec-1.
    return rlzero


def dayrh(tt: float, tdew: float) -> float:
    """Computes the hourly values of relative humidity, using the hourly air and dew point temperatures. It calls function `VaporPressure`

    If the estimated dew point is higher than the actual air temperature, its value is taken as the air temperature (relative humidity 100%).

    The relative humidity is calculated as the percentage ratio of the saturated vapor pressure at dew point temperature and the saturated vapor pressure at actual air temperature.

    Reference:

    Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling diurnal patterns of air temperature, radiation, wind speed and relative humidity by equations from daily characteristics. Agricultural Systems 51:377-393.
    :param tt: air temperature C at this time of day.
    :type tt: float
    :param tdew: dew point temperature C at this time of day.
    :type tdew: float
    :return: relative humidity
    :rtype: float
    """
    td = min(tt, tdew)  # the dew point temperature (C), is assumed to be tt if tt < tdew.
    esvp = VaporPressure(tt)  # the saturated vapor pressure in the air (mbar).
    vpa = VaporPressure(td)  # the actual vapor pressure in the air (mbar).
    relative_humidity = 100 * vpa / esvp  # relative humidity at this time of day, %.
    return min(100, max(1, relative_humidity))


cdef class Simulation:
    cdef cSimulation _sim
    cdef public unsigned int profile_id
    cdef public unsigned int version
    cdef double relative_radiation_received_by_a_soil_column[20]  # the relative radiation received by a soil column, as affected by shading by plant canopy.
    cdef double max_leaf_area_index
    cdef double ptsred  # The effect of moisture stress on the photosynthetic rate
    cdef double DaysTo1stSqare   # number of days from emergence to 1st square
    cdef double defkgh  # amount of defoliant applied, kg per ha
    cdef double tdfkgh  # total cumulative amount of defoliant
    cdef bool_t idsw  # switch indicating if predicted defoliation date was defined.
    cdef public double skip_row_width  # the smaller distance between skip rows, cm
    cdef public double plants_per_meter  # average number of plants pre meter of row.

    def __init__(self, profile_id=0, version=0x0400):
        self.profile_id = profile_id
        self.version = version
        self.max_leaf_area_index = 0.001

    @property
    def year(self):
        return self._sim.year

    @year.setter
    def year(self, year):
        self._sim.year = year

    @property
    def start_date(self):
        return doy2date(self.year, self._sim.day_start)

    @start_date.setter
    def start_date(self, d):
        self._sim.day_start = date2doy(d)

    @property
    def stop_date(self):
        return doy2date(self.year, self._sim.day_finish)

    @stop_date.setter
    def stop_date(self, d):
        self._sim.day_finish = date2doy(d)

    @property
    def emerge_date(self):
        return doy2date(self.year, self._sim.day_emerge)

    @emerge_date.setter
    def emerge_date(self, d):
        self._sim.day_emerge = date2doy(d)

    @property
    def plant_date(self):
        return doy2date(self.year, self._sim.day_plant)

    @plant_date.setter
    def plant_date(self, d):
        self._sim.day_plant = date2doy(d)

    @property
    def topping_date(self):
        if self.version >= 0x500:
            return doy2date(self.year, self._sim.day_topping)
        return None

    @topping_date.setter
    def topping_date(self, d):
        self._sim.day_topping = date2doy(d)

    @property
    def latitude(self):
        return self._sim.latitude

    @latitude.setter
    def latitude(self, value):
        self._sim.latitude = value or 0

    @property
    def longitude(self):
        return self._sim.longitude

    @longitude.setter
    def longitude(self, value):
        self._sim.longitude = value or 0

    @property
    def elevation(self):
        return self._sim.elevation

    @elevation.setter
    def elevation(self, value):
        self._sim.elevation = value or 0

    @property
    def site_parameters(self):
        return SitePar

    @site_parameters.setter
    def site_parameters(self, parameters):
        for i, p in enumerate(parameters):
            SitePar[i + 1] = p

    @property
    def cultivar_parameters(self):
        return self._sim.cultivar_parameters

    @cultivar_parameters.setter
    def cultivar_parameters(self, parameters):
        for i, p in enumerate(parameters):
            self._sim.cultivar_parameters[i + 1] = p

    @property
    def row_space(self):
        return self._sim.row_space

    @row_space.setter
    def row_space(self, value):
        self._sim.row_space = value or 0

    @property
    def first_square_date(self):
        return doy2date(self.year, self._sim.first_square)

    @first_square_date.setter
    def first_square_date(self, value):
        self._sim.first_square = date2doy(value)

    @property
    def defoliate_date(self):
        return doy2date(self.year, self._sim.day_defoliate)

    @defoliate_date.setter
    def defoliate_date(self, value):
        self._sim.day_defoliate = date2doy(value)

    @property
    def states(self):
        return [self.state(i) for i in
                range(self._sim.day_finish - self._sim.day_start + 1)]

    def state(self, i):
        return State.from_ptr(&self._sim.states[i], self.year, self.version)

    @property
    def climate(self):
        return self._sim.climate

    @climate.setter
    def climate(self, climate):
        alias = {
            "radiation": "Rad",
            "max": "Tmax",
            "min": "Tmin",
            "wind": "Wind",
            "rain": "Rain",
            "dewpoint": "Tdew",
        }
        for i, daily_climate in enumerate(climate):
            self._sim.climate[i] = {
                alias[k]: v for k, v in daily_climate.items()
            }

    def run(self):
        self._init_state()
        try:
            self._simulate()
        except RuntimeError:
            pass

    def _init_state(self):
        cdef State state0 = self.state(0)
        state0.date = self.start_date
        state0.lint_yield = 0
        state0.soil.number_of_layers_with_root = 7
        state0.plant_height = 4.0
        state0.plant_weight = 0
        state0.stem_weight = 0.2
        state0.petiole_weight = 0
        state0.square_weight = 0
        state0.green_bolls_weight = 0
        state0.green_bolls_burr_weight = 0
        state0.open_bolls_weight = 0
        state0.open_bolls_burr_weight = 0
        state0.bloom_weight_loss = 0
        state0.abscised_fruit_sites = 0
        state0.abscised_leaf_weight = 0
        state0.cumulative_nitrogen_loss = 0
        state0.cumulative_transpiration = 0
        state0.cumulative_evaporation = 0
        state0.applied_water = 0
        state0.water_stress = 1
        state0.water_stress_stem = 1
        state0.carbon_stress = 1
        state0.extra_carbon = 0
        state0.leaf_area_index = 0.001
        state0.leaf_area = 0
        state0.leaf_weight = 0.20
        state0.leaf_nitrogen = 0.0112
        state0.number_of_vegetative_branches = 1
        state0.number_of_squares = 0
        state0.number_of_green_bolls = 0
        state0.number_of_open_bolls = 0
        state0.nitrogen_stress = 1
        state0.nitrogen_stress_vegetative = 1
        state0.nitrogen_stress_fruiting = 1
        state0.nitrogen_stress_root = 1
        state0.total_required_nitrogen = 0
        state0.leaf_nitrogen_concentration = .056
        state0.petiole_nitrogen_concentration = 0
        state0.seed_nitrogen_concentration = 0
        state0.burr_nitrogen = 0
        state0.seed_nitrogen = 0
        state0.root_nitrogen_concentration = .026
        state0.square_nitrogen_concentration = 0
        state0.square_nitrogen = 0
        state0.stem_nitrogen = 0.0072
        state0.fruit_growth_ratio = 1
        state0.ginning_percent = 0.35
        state0.number_of_pre_fruiting_nodes = 1
        for i in range(9):
            state0.age_of_pre_fruiting_nodes[i] = 0
            state0.leaf_area_pre_fruiting[i] = 0
            state0.leaf_weight_pre_fruiting[i] = 0
        for k in range(3):
            state0._[0].vegetative_branches[k].number_of_fruiting_branches = 0
            for l in range(30):
                state0._[0].vegetative_branches[k].fruiting_branches[
                    l].number_of_fruiting_nodes = 0
                state0._[0].vegetative_branches[k].fruiting_branches[l].delay_for_new_node = 0
                state0._[0].vegetative_branches[k].fruiting_branches[l].main_stem_leaf = dict(
                    leaf_area=0,
                    leaf_weight=0,
                    petiole_weight=0,
                    potential_growth_for_leaf_area=0,
                    potential_growth_for_leaf_weight=0,
                    potential_growth_for_petiole_weight=0,
                )
                for m in range(5):
                    state0._[0].vegetative_branches[k].fruiting_branches[l].nodes[m] = dict(
                        age=0,
                        fraction=0,
                        average_temperature=0,
                        ginning_percent=0.35,
                        stage=Stage.NotYetFormed,
                        leaf=dict(
                            age=0,
                            potential_growth=0,
                            area=0,
                            weight=0,
                        ),
                        square=dict(
                            potential_growth=0,
                            weight=0,
                        ),
                        boll=dict(
                            age=0,
                            potential_growth=0,
                            weight=0,
                        ),
                        burr=dict(
                            potential_growth=0,
                            weight=0,
                        ),
                        petiole=dict(
                            potential_growth=0,
                            weight=0,
                        ),
                    )
        self._sim.states[0] = state0._[0]

    def _simulate(self):
        for i in range(self._sim.day_finish - self._sim.day_start + 1):
            self._simulate_this_day(i)
            if i < self._sim.day_finish - self._sim.day_start:
                CopyState(self._sim, i)

    def _energy_balance(self, u, ihr, k, ess, etp1):
        """
        This function solves the energy balance equations at the soil surface, and at the foliage / atmosphere interface. It computes the resulting temperatures of the soil surface and the plant canopy.

        Units for all energy fluxes are: cal cm-2 sec-1.
        It is called from SoilTemperature(), on each hourly time step and for each soil column.
        It calls functions clearskyemiss(), VaporPressure(), SensibleHeatTransfer(), SoilSurfaceBalance() and CanopyBalance.()

        :param ihr: the time of day in hours.
        :param k: soil column number.
        :param ess: evaporation from surface of a soil column (mm / sec).
        :param etp1: actual transpiration rate (mm / sec).
        :param sf: fraction of shaded soil area
        """
        # Constants used:
        cdef double wndfac = 0.60  # Ratio of wind speed under partial canopy cover.
        cdef double cswint = 0.75  # proportion of short wave radiation (on fully shaded soil surface) intercepted by the canopy.
        # Set initial values
        cdef double sf = 1 - self.relative_radiation_received_by_a_soil_column[k]
        cdef double thet = self._sim.states[u].hours[ihr].temperature + 273.161  # air temperature, K
        cdef double so = SoilTemp[0][k]  # soil surface temperature, K
        cdef double so2 = SoilTemp[1][k]  # 2nd soil layer temperature, K
        cdef double so3 = SoilTemp[2][k]  # 3rd soil layer temperature, K
        # Compute soil surface albedo (based on Horton and Chung, 1991):
        ag = compute_soil_surface_albedo(self._sim.states[u].soil.cells[0][k].water_content, FieldCapacity[0], thad[0], SitePar[15], SitePar[16])

        rzero, rss, rsup = compute_incoming_short_wave_radiation(self._sim.states[u].hours[ihr].radiation, sf * cswint, ag)
        rlzero = compute_incoming_long_wave_radiation(self._sim.states[u].hours[ihr].humidity, self._sim.states[u].hours[ihr].temperature, self._sim.states[u].hours[ihr].cloud_cov, self._sim.states[u].hours[ihr].cloud_cor)

        # Set initial values of canopy temperature and air temperature in canopy.
        cdef double tv  # temperature of plant foliage (K)
        cdef double tafk  # temperature (K) of air inside the canopy.
        if sf < 0.05:  # no vegetation
            tv = thet
            tafk = thet
        # Wind velocity in canopy is converted to cm / s.
        cdef double wndhr  # wind speed in cm /sec
        wndhr = self._sim.states[u].hours[ihr].wind_speed * 100
        cdef double rocp  # air density * specific heat at constant pressure = 0.24 * 2 * 1013 / 5740
        # divided by tafk.
        cdef double c2  # multiplier for sensible heat transfer (at plant surface).
        cdef double rsv  # global radiation absorbed by the vegetation
        if sf >= 0.05:  # a shaded soil column
            tv = FoliageTemp[k]  # vegetation temperature
            # Short wave radiation intercepted by the canopy:
            rsv = (
                    rzero * (1 - self._sim.states[u].hours[ihr].albedo) * sf * cswint  # from above
                    + rsup * (1 - self._sim.states[u].hours[ihr].albedo) * sf * cswint  # reflected from soil surface
            )
            # Air temperature inside canopy is the average of soil, air, and plant temperatures, weighted by 0.1, 0.3, and 0.6, respectively.
            tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv)

            # Call SensibleHeatTransfer() to compute sensible heat transfer coefficient. Factor 2.2 for sensible heat transfer: 2 sides of leaf plus stems and petioles.
            # sensible heat transfer coefficient for soil
            varcc = SensibleHeatTransfer(tv, tafk, self._sim.states[u].plant_height, wndhr)  # canopy to air
            rocp = 0.08471 / tafk
            c2 = 2.2 * sf * rocp * varcc
        cdef double soold = so  # previous value of soil surface temperature
        cdef double tvold = tv  # previous value of vegetation temperature
        # Starting iterations for soil and canopy energy balance
        for menit in range(30):
            soold = so
            wndcanp = (1 - sf * (1 - wndfac)) * wndhr  # estimated wind speed under canopy
            # Call SensibleHeatTransfer() to compute sensible heat transfer for soil surface to air
            tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv)
            # sensible heat transfer coefficientS for soil
            varc = SensibleHeatTransfer(so, tafk, 0, wndcanp)
            rocp = 0.08471 / tafk
            hsg = rocp * varc  # multiplier for computing sensible heat transfer soil to air.
            # Call SoilSurfaceBalance() for energy balance in soil surface / air interface.
            SoilSurfaceBalance(self._sim.states[u], ihr, k, ess, rlzero, rss, sf, hsg, so, so2, so3, thet, tv)

            if sf >= 0.05:
                # This section executed for shaded columns only.
                tvold = tv
                # Compute canopy energy balance for shaded columns
                CanopyBalance(ihr, k, etp1, rlzero, rsv, c2, sf, so, thet, tv, self._sim.day_start + u)
                if menit >= 10:
                    # The following is used to reduce fluctuations.
                    so = (so + soold) / 2
                    tv = (tv + tvold) / 2
            if abs(tv - tvold) <= 0.05 and abs(so - soold) <= 0.05:
                break
        else:
            raise SimulationEnd  # If more than 30 iterations are needed - stop simulation.
        # After convergence - set global variables for the following temperatures:
        if sf >= 0.05:
            FoliageTemp[k] = tv
        SoilTemp[0][k] = so
        SoilTemp[1][k] = so2
        SoilTemp[2][k] = so3

    def _soil_temperature(self, u):
        """
        This is the main part of the soil temperature sub-model.
        It is called daily from self._simulate_this_day.
        It calls the following functions:
        _energy_balance(), PredictEmergence(), SoilHeatFlux(), SoilTemperatureInit().

        References:

        Benjamin, J.G., Ghaffarzadeh, M.R. and Cruse, R.M. 1990. Coupled water and heat transport in ridged soils. Soil Sci. Soc. Am. J. 54:963-969.

        Chen, J. 1984. Uncoupled multi-layer model for the transfer of sensible and latent heat flux densities from vegetation. Boundary-Layer Meteorology 28:213-225.

        Chen, J. 1985. A graphical extrapolation method to determine canopy resistance from measured temperature and humidity profiles above a crop canopy. Agric. For. Meteorol. 37:75-88.

        Clothier, B.E., Clawson, K.L., Pinter, P.J.Jr., Moran, M.S., Reginato, R.J. and Jackson, R.D. 1986. Estimation of soil heat flux from net radiation during the growth of alfalfa. Agric. For. Meteorol. 37:319-329.

        Costello, T.A. and Braud, H.J. Jr. 1989. Thermal diffusivity of soil by nonlinear regression analysis of soil temperature data. Trans. ASAE 32:1281-1286.

        De Vries, D.A. 1963. Thermal properties of soils. In: W.R. Van Wijk (ed) Physics of plant environment, North Holland, Amsterdam, pp 210-235.

        Deardorff, J.W. 1978. Efficient prediction of ground surface temperature and moisture with inclusion of a layer of vegetation. J. Geophys. Res. 83 (C4):1889-1903.

        Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of daily and hourly net radiation. CIMIS Final Report June 1988, pp. 58-79.

        Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling diurnal patterns of air temperature, radiation, wind speed and relative humidity by equations from daily characteristics. Agricultural Systems 51:377-393.

        Hadas, A. 1974. Problem involved in measuring the soil thermal conductivity and diffusivity in a moist soil. Agric. Meteorol. 13:105-113.

        Hadas, A. 1977. Evaluation of theoretically predicted thermal conductivities of soils under field and laboratory conditions. Soil Sci. Soc. Am. J. 41:460-466.

        Hanks, R.J., Austin, D.D. and Ondrechen, W.T. 1971. Soil temperature estimation by a numerical method. Soil Sci. Soc. Am. Proc. 35:665-667.

        Hares, M.A. and Novak, M.D. 1992. Simulation of surface energy balance and soil temperature under strip tillage: I. Model description. Soil Sci. Soc. Am. J. 56:22-29.

        Hares, M.A. and Novak, M.D. 1992. Simulation of surface energy balance and soil temperature under strip tillage: II. Field test. Soil Sci. Soc. Am. J. 56:29-36.

        Horton, E. and Wierenga, P.J. 1983. Estimating the soil heat flux from observations of soil temperature near the surface. Soil Sci. Soc. Am. J. 47:14-20.

        Horton, E., Wierenga, P.J. and Nielsen, D.R. 1983. Evaluation of methods for determining apparent thermal diffusivity of soil near the surface. Soil Sci. Soc. Am. J. 47:25-32.

        Horton, R. 1989. Canopy shading effects on soil heat and water flow. Soil Sci. Soc. Am. J. 53:669-679.

        Horton, R., and Chung, S-O, 1991. Soil Heat Flow. Ch. 17 in: Hanks, J., and Ritchie, J.T., (Eds.) Modeling Plant and Soil Systems. Am. Soc. Agron., Madison, WI, pp 397-438.

        Iqbal, M. 1983. An Introduction to Solar Radiation. Academic Press.

        Kimball, B.A., Jackson, R.D., Reginato, R.J., Nakayama, F.S. and Idso, S.B. 1976. Comparison of field-measured and calculated soil heat fluxes. Soil Sci. Soc. Am. J. 40:18-28.

        Lettau, B. 1971. Determination of the thermal diffusivity in the upper layers of a natural ground cover. Soil Sci. 112:173-177.

        Monin, A.S. 1973. Boundary layers in planetary atmospheres. In: P. Morrel (ed.), Dynamic meteorology, D. Reidel Publishing Company, Boston, pp. 419-458.

        Spitters, C.J.T., Toussaint, H.A.J.M. and Goudriaan, J. 1986. Separating the diffuse and direct component of global radiation and its implications for modeling canopy photosynthesis. Part I. Components of incoming radiation. Agric. For. Meteorol. 38:217-229.

        Wierenga, P.J. and de Wit, C.T. 1970. Simulation of heat flow in soils. Soil Sci. Soc. Am. Proc. 34:845-848.

        Wierenga, P.J., Hagan, R.M. and Nielsen, D.R. 1970. Soil temperature profiles during infiltration and redistribution of cool and warm irrigation water. Water Resour. Res. 6:230-238.

        Wierenga, P.J., Nielsen, D.R. and Hagan, R.M. 1969. Thermal properties of soil based upon field and laboratory measurements. Soil Sci. Soc. Am. Proc. 33:354-360.
        """
        global DeepSoilTemperature
        if u == 0:
            SoilTemperatureInit(self._sim)
        # Compute dts, the daily change in deep soil temperature (C), as a site-dependent function of Daynum.
        cdef double dts = 2 * pi * SitePar[10] / 365 * cos(2 * pi * (self._sim.states[u].daynum - SitePar[11]) / 365)
        # Define iter1 and dlt for hourly time step.
        cdef int iter1 = 24  # number of iterations per day.
        cdef double dlt = 3600  # time (seconds) of one iteration.
        cdef int kk = 1  # number of soil columns for executing computations.
        # If there is no canopy cover, no horizontal heat flux is assumed, kk = 1.
        # Otherwise it is equal to the number of columns in the slab.
        cdef double shadeav = 0  # average shaded area in all shaded soil columns.
        # isw defines the type of soil temperature computation.
        if isw > 1:
            shadetot = 0  # sum of shaded area in all shaded soil columns.
            nshadedcol = 0  # number of at least partially shaded soil columns.
            kk = nk
            for k in range(nk):
                if self.relative_radiation_received_by_a_soil_column[k] <= 0.99:
                    shadetot += 1 - self.relative_radiation_received_by_a_soil_column[k]
                    nshadedcol += 1

            if nshadedcol > 0:
                shadeav = shadetot / nshadedcol
        # Set daily averages of soil temperature to zero.
        for l in range(nl):
            for k in range(nk):
                SoilTempDailyAvrg[l][k] = 0
        # es and ActualSoilEvaporation are computed as the average for the whole soil slab, weighted by column widths.
        cdef double es = 0  # potential evaporation rate, mm day-1
        self._sim.states[u].potential_evaporation = 0
        self._sim.states[u].actual_soil_evaporation = 0
        # Start hourly loop of iterations.
        for ihr in range(iter1):
            # Update the temperature of the last soil layer (lower boundary conditions).
            DeepSoilTemperature += dts * dlt / 86400
            etp0 = 0  # actual transpiration (mm s-1) for this hour
            if self._sim.states[u].evapotranspiration > 0.000001:
                etp0 = self._sim.states[u].actual_transpiration * self._sim.states[u].hours[ihr].ref_et / self._sim.states[u].evapotranspiration / dlt
            # Compute vertical transport for each column
            for k in range(kk):
                #  Set SoilTemp for the lowest soil layer.
                SoilTemp[nl - 1][k] = DeepSoilTemperature
                # Compute transpiration from each column, weighted by its relative shading.
                etp1 = 0  # actual hourly transpiration (mm s-1) for a column.
                if shadeav > 0.000001:
                    etp1 = etp0 * (1 - self.relative_radiation_received_by_a_soil_column[k]) / shadeav
                ess = 0  # evaporation rate from surface of a soil column (mm / sec).
                # The potential evaporation rate (escol1k) from a column is the sum of the radiation component of the Penman equation(es1hour), multiplied by the relative radiation reaching this column, and the wind and vapor deficit component of the Penman equation (es2hour).
                # potential evaporation fron soil surface of a column, mm per hour.
                escol1k = self._sim.states[u].hours[ihr].et1 * self.relative_radiation_received_by_a_soil_column[k] + self._sim.states[u].hours[ihr].et2
                es += escol1k * wk(k, self.row_space)
                # Compute actual evaporation from soil surface. update cell.water_content of the soil soil cell, and add to daily sum of actual evaporation.
                evapmax = 0.9 * (self._sim.states[u].soil.cells[0][k].water_content - thad[0]) * 10 * dl(0)  # maximum possible evaporatio from a soil cell near the surface.
                escol1k = min(evapmax, escol1k)
                self._sim.states[u].soil.cells[0][k].water_content -= 0.1 * escol1k / dl(0)
                self._sim.states[u].actual_soil_evaporation += escol1k * wk(k, self.row_space)
                ess = escol1k / dlt
                # Call self._energy_balance to compute soil surface and canopy temperature.
                self._energy_balance(u, ihr, k, ess, etp1)
            # Compute soil temperature flux in the vertical direction.
            # Assign iv = 1, layer = 0, nn = nl.
            iv = 1  # indicates vertical (=1) or horizontal (=0) flux.
            nn = nl  # number of array members for heat flux.
            layer = 0  # soil layer number
            # Loop over kk columns, and call SoilHeatFlux().
            for k in range(kk):
                SoilHeatFlux(self._sim.states[u], dlt, iv, nn, layer, k, self.row_space)
            # If no horizontal heat flux is assumed, make all array members of SoilTemp equal to the value computed for the first column. Also, do the same for array memebers of cell.water_content.
            if isw <= 1:
                for l in range(nl):
                    for k in range(nk):
                        SoilTemp[l][k] = SoilTemp[l][0]
                        if l == 0:
                            self._sim.states[u].soil.cells[l][k].water_content = self._sim.states[u].soil.cells[l][0].water_content
            # Compute horizontal transport for each layer

            # Compute soil temperature flux in the horizontal direction, when isw = 2.
            # Assign iv = 0 and nn = nk. Start loop for soil layers, and call SoilHeatFlux.
            if isw > 1:
                iv = 0
                nn = nk
                for l in range(nl):
                    layer = l
                    SoilHeatFlux(self._sim.states[u], dlt, iv, nn, layer, l, self.row_space)
            # Compute average temperature of soil layers, in degrees C.
            tsolav = [0] * nl  # hourly average soil temperature C, of a soil layer.
            for l in range(nl):
                for k in range(nk):
                    SoilTempDailyAvrg[l][k] += SoilTemp[l][k]
                    tsolav[l] += SoilTemp[l][k] - 273.161
                tsolav[l] /= nk
            # Compute average temperature of foliage, in degrees C. The average is weighted by the canopy shading of each column, only columns which are shaded 5% or more by canopy are used.
            tfc = 0  # average foliage temperature, weighted by shading in each column
            shading = 0  # sum of shaded area in all shaded columns, used to compute TFC
            for k in range(nk):
                if self.relative_radiation_received_by_a_soil_column[k] <= 0.95:
                    tfc += (FoliageTemp[k] - 273.161) * (1 - self.relative_radiation_received_by_a_soil_column[k])
                    shading += 1 - self.relative_radiation_received_by_a_soil_column[k]
            if shading >= 0.01:
                tfc = tfc / shading
            # If emergence date is to be simulated, call PredictEmergence().
            if isw == 0 and self._sim.states[u].daynum >= self._sim.day_plant:
                PredictEmergence(self._sim, u, ihr)
        # At the end of the day compute actual daily evaporation and its cumulative sum.
        if kk == 1:
            es /= wk(1, self.row_space)
            self._sim.states[u].actual_soil_evaporation /= wk(1, self.row_space)
        else:
            es /= self.row_space
            self._sim.states[u].actual_soil_evaporation /= self.row_space
        self._sim.states[u].cumulative_evaporation += self._sim.states[u].actual_soil_evaporation
        if self._sim.states[u].kday > 0:
            self._sim.states[u].potential_evaporation = es
        # compute daily averages.
        for l in range(nl):
            for k in range(nk):
                SoilTempDailyAvrg[l][k] /= iter1

    def _column_shading(self, u):
        cdef cState state = self._sim.states[u]
        zint = 1.0756 * state.plant_height / self.row_space
        sw = 0
        for k in range(nk):
            if k <= self._sim.plant_row_column:
                j = self._sim.plant_row_column - k
                sw += wk(j, self.row_space)
                sw0 = sw
                sw1 = sw - wk(j, self.row_space) / 2
                k0 = j
            else:
                sw += wk(k, self.row_space)
                sw1 = sw - sw0 - wk(k, self.row_space) / 2
                k0 = k
            shade = 0
            if sw1 < state.plant_height:
                shade = 1 - (sw1 / state.plant_height) * (sw1 / state.plant_height)
                if state.light_interception < zint and state.leaf_area_index < self.max_leaf_area_index:
                    shade *= state.light_interception / zint
            self.relative_radiation_received_by_a_soil_column[k0] = max(0.05, 1 - shade)

    def _daily_climate(self, u):
        cdef double declination  # daily declination angle, in radians
        cdef double sunr  # time of sunrise, hours.
        cdef double suns  # time of sunset, hours.
        cdef double tmpisr  # extraterrestrial radiation, \frac{W}{m^2}
        hour = timedelta(hours=1)
        result = compute_day_length((self.latitude, self.longitude), doy2date(self.year, self._sim.states[u].daynum))
        declination = result["declination"]
        zero = result["sunr"].replace(hour=0, minute=0, second=0, microsecond=0)
        sunr = (result["sunr"] - zero) / hour
        suns = (result["suns"] - zero) / hour
        tmpisr = result["tmpisr"]
        self._sim.states[u].solar_noon = (result["solar_noon"] - zero) / hour
        self._sim.states[u].day_length = result["day_length"] / hour

        cdef double xlat = self.latitude * pi / 180  # latitude converted to radians.
        cdef double cd = cos(xlat) * cos(declination)  # amplitude of the sine of the solar height.
        cdef double sd = sin(xlat) * sin(declination)  # seasonal offset of the sine of the solar height.
        # The computation of the daily integral of global radiation (from sunrise to sunset) is based on Spitters et al. (1986).
        cdef double c11 = 0.4  # constant parameter
        cdef double radsum
        if abs(sd / cd) >= 1:
            radsum = 0
        else:
            # dsbe is the integral of sinb * (1 + c11 * sinb) from sunrise to sunset.
            dsbe = acos(-sd / cd) * 24 / pi * (sd + c11 * sd * sd + 0.5 * c11 * cd * cd) + 12 * (
                    cd * (2 + 3 * c11 * sd)) * sqrt(1 - (sd / cd) * (sd / cd)) / pi
            # The daily radiation integral is computed for later use in function Radiation.
            # Daily radiation intedral is converted from langleys to Watt m - 2, and divided by dsbe.
            # 11.630287 = 1000000 / 3600 / 23.884
            radsum = self._sim.climate[u].Rad * 11.630287 / dsbe
        cdef double rainToday
        rainToday = self._sim.climate[u].Rain  # the amount of rain today, mm
        # Set 'pollination switch' for rainy days (as in GOSSYM).
        self._sim.states[u].pollination_switch = rainToday < 2.5
        # Call SimulateRunoff() only if the daily rainfall is more than 2 mm.
        # Note: this is modified from the original GOSSYM - RRUNOFF routine. It is called here for rainfall only, but it is not activated when irrigation is applied.
        cdef double runoffToday = 0  # amount of runoff today, mm
        if rainToday >= 2.0:
            runoffToday = SimulateRunoff(self._sim, u, SandVolumeFraction[0], ClayVolumeFraction[0], NumIrrigations)
            if runoffToday < rainToday:
                rainToday -= runoffToday
            else:
                rainToday = 0
            self._sim.climate[u].Rain = rainToday
        self._sim.states[u].runoff = runoffToday
        # Parameters for the daily wind function are now computed:
        cdef double t1 = sunr + SitePar[1]  # the hour at which wind begins to blow (SitePar(1) hours after sunrise).
        cdef double t2 = self._sim.states[u].solar_noon + SitePar[
            2]  # the hour at which wind speed is maximum (SitePar(2) hours after solar noon).
        cdef double t3 = suns + SitePar[3]  # the hour at which wind stops to blow (SitePar(3) hours after sunset).
        cdef double wnytf = SitePar[4]  # used for estimating night time wind (from time t3 to time t1 next day).

        for ihr in range(24):
            ti = ihr + 0.5
            sinb = sd + cd * cos(pi * (ti - self._sim.states[u].solar_noon) / 12)
            self._sim.states[u].hours[ihr].radiation = radiation(radsum, sinb, c11)
            self._sim.states[u].hours[ihr].temperature = daytmp(self._sim, u, ti, SitePar[8], LastDayWeatherData, sunr,
                                                                suns)
            self._sim.states[u].hours[ihr].dew_point = tdewhour(self._sim, u, LastDayWeatherData, ti,
                                                                self._sim.states[u].hours[ihr].temperature, sunr,
                                                                self._sim.states[u].solar_noon, SitePar[8], SitePar[12],
                                                                SitePar[13], SitePar[14])
            self._sim.states[u].hours[ihr].humidity = dayrh(self._sim.states[u].hours[ihr].temperature,
                                                            self._sim.states[u].hours[ihr].dew_point)
            self._sim.states[u].hours[ihr].wind_speed = daywnd(ti, self._sim.climate[u].Wind, t1, t2, t3, wnytf)
        # Compute average daily temperature, using function AverageAirTemperatures.
        AverageAirTemperatures(self._sim.states[u].hours, self._sim.states[u].average_temperature, DayTimeTemp,
                               NightTimeTemp)
        # Compute potential evapotranspiration.
        EvapoTranspiration(self._sim.states[u], self.latitude, self.elevation, declination, tmpisr, SitePar[7])

    def _stress(self, u):
        global AverageLwpMin, AverageLwp, LwpMin, LwpMax
        # The following constant parameters are used:
        cdef double[9] vstrs = [-3.0, 3.229, 1.907, 0.321, -0.10, 1.230, 0.340, 0.30, 0.05]
        # Call LeafWaterPotential() to compute leaf water potentials.
        LeafWaterPotential(self._sim.states[u], self.row_space)
        # The running averages, for the last three days, are computed:
        # AverageLwpMin is the average of LwpMin, and AverageLwp of LwpMin + LwpMax.
        AverageLwpMin += (LwpMin - LwpMinX[2]) / 3
        AverageLwp += (LwpMin + LwpMax - LwpX[2]) / 3
        for i in (2, 1):
            LwpMinX[i] = LwpMinX[i - 1]
            LwpX[i] = LwpX[i - 1]
        LwpMinX[0] = LwpMin
        LwpX[0] = LwpMin + LwpMax
        if self._sim.states[u].kday < 5:
            self.ptsred = 1
            self._sim.states[u].water_stress_stem = 1
            return
        # The computation of ptsred, the effect of moisture stress on the photosynthetic rate, is based on the following work:
        # Ephrath, J.E., Marani, A., Bravdo, B.A., 1990. Effects of moisture stress on stomatal resistance and photosynthetic rate in cotton (Gossypium hirsutum) 1. Controlled levels of stress. Field Crops Res.23:117-131.
        # It is a function of AverageLwpMin (average LwpMin for the last three days).
        if AverageLwpMin < vstrs[0]:
            AverageLwpMin = vstrs[0]
        self.ptsred = vstrs[1] + AverageLwpMin * (vstrs[2] + vstrs[3] * AverageLwpMin)
        if self.ptsred > 1:
            self.ptsred = 1
        # The general moisture stress factor (WaterStress) is computed as an empirical function of AverageLwp. psilim, the value of AverageLwp at the maximum value of the function, is used for truncating it.
        # The minimum value of WaterStress is 0.05, and the maximum is 1.
        cdef double psilim  # limiting value of AverageLwp.
        cdef double WaterStress
        psilim = -0.5 * vstrs[5] / vstrs[6]
        if AverageLwp > psilim:
            WaterStress = 1
        else:
            WaterStress = vstrs[4] - AverageLwp * (vstrs[5] + vstrs[6] * AverageLwp)
            if WaterStress > 1:
                WaterStress = 1
            if WaterStress < 0.05:
                WaterStress = 0.05
        # Water stress affecting plant height and stem growth(WaterStressStem) is assumed to be more severe than WaterStress, especially at low WaterStress values.
        self._sim.states[u].water_stress_stem = WaterStress * (1 + vstrs[7] * (2 - WaterStress)) - vstrs[7]
        if self._sim.states[u].water_stress_stem < vstrs[8]:
            self._sim.states[u].water_stress_stem = vstrs[8]
        self._sim.states[u].water_stress = WaterStress

    def _get_net_photosynthesis(self, u):
        """
        References:
        Baker et. al. (1972). Simulation of Growth and Yield in Cotton: I. Gross photosynthesis, respiration and growth. Crop Sci. 12:431-435.
        Harper et. al. (1973) Carbon dioxide and the photosynthesis of field crops.  A metered carbon dioxide release in cotton under field conditions.  Agron. J. 65:7-11.
        Baker (1965)  Effects of certain environmental factors on net assimilation in cotton.  Crop Sci. 5:53-56 (Fig 5).
        """
        global NetPhotosynthesis, CumNetPhotosynth
        cdef cState state = self._sim.states[u]
        # constants
        cdef double gsubr = 0.375  # the growth respiration factor.
        cdef double rsubo = 0.0032  # maintenance respiration factor.
        cdef double[4] vpnet = [1.30, 0.034, 0.010, 0.32]
        # Note: co2parm is for icrease in ambient CO2 concentration changes from 1959 (308 ppm).
        # The first 28 values (up to 1987) are from GOSSYM. The other values (up to 2004) are derived from data of the Carbon Dioxide Information Analysis Center (CDIAC).

        # Exit the function and end simulation if there are no leaves
        if state.leaf_area_index <= 0:
            raise SimulationEnd
        # Get the CO2 correction factor (pnetcor) for photosynthesis, using ambient_co2_factor and a factor that may be variety specific (vpnet[0]).
        cdef double pnetcor = ambient_co2_factor(self.year) * vpnet[0]  # correction factor for gross photosynthesis
        # Compute ptnfac, the effect of leaf N concentration on photosynthesis, using an empirical relationship.
        cdef double ptnfac = vpnet[3] + (state.leaf_nitrogen_concentration - vpnet[2]) * (1 - vpnet[3]) / (
                vpnet[1] - vpnet[2])  # correction factor for low nitrogen content in leaves.
        if ptnfac > 1:
            ptnfac = 1
        if ptnfac < vpnet[3]:
            ptnfac = vpnet[3]

        # Convert the average daily short wave radiation from langley per day, to Watts per square meter (wattsm).
        cdef double wattsm = self._sim.climate[u].Rad * 697.45 / (
                state.day_length * 60)  # average daily global radiation, W m^{-2}
        # Compute pstand as an empirical function of wattsm (based on Baker et al., 1972)
        cdef double pstand  # gross photosynthesis for a non-stressed full canopy
        pstand = 2.3908 + wattsm * (1.37379 - wattsm * 0.00054136)
        # Convert it to gross photosynthesis per plant (pplant), using per_plant_area and corrections for light interception by canopy, ambient CO2 concentration, water stress and low N in the leaves.
        cdef double pplant  # actual gross photosynthetic rate, g per plant per day.
        pplant = 0.001 * pstand * state.light_interception * self._sim.per_plant_area * self.ptsred * pnetcor * ptnfac
        # Compute the photorespiration factor (rsubl) as a linear function of average day time temperature.
        cdef double rsubl = 0.0032125 + 0.0066875 * DayTimeTemp  # photorespiration factor
        # Photorespiration (lytres) is computed as a proportion of gross photosynthetic rate.
        cdef double lytres  # rate of photorespiration, g per plant per day.
        lytres = rsubl * pplant
        # Old stems are those more than voldstm = 32 calendar days old.
        # Maintenance respiration is computed on the basis of plant dry weight, minus the old stems and the dry tissue of opened bolls.
        cdef oldstmwt  # weight of old stems
        cdef int voldstm = 32
        cdef int kkday = state.kday - voldstm  # day of least recent actively growing stems.
        if kkday < 1:
            oldstmwt = 0
        else:
            oldstmwt = StemWeight[kkday]
        cdef double bmain  # maintenance respiration, g per plant per day.
        bmain = (state.plant_height - state.open_bolls_weight - state.open_bolls_burr_weight - oldstmwt) * rsubo
        # Net photosynthesis is computed by subtracting photo-respiration and maintenance respiration from the gross rate of photosynthesis.
        # To avoid computational problem, make sure that pts is positive and non-zero.
        cdef double pts  # intermediate computation of NetPhotosynthesis.
        pts = pplant - lytres - bmain
        if pts < 0.00001:
            pts = 0.00001
        # the growth respiration (gsubr) supplies energy for converting the supplied carbohydrates to plant tissue dry matter.
        # 0.68182 converts CO2 to CH2O. NetPhotosynthesis is the computed net photosynthesis, in g per plant per day.
        NetPhotosynthesis = pts / (1 + gsubr) * 0.68182
        # CumNetPhotosynth is the cumulative value of NetPhotosynthesis, from day of emergence.
        CumNetPhotosynth += NetPhotosynthesis

    def _potential_fruit_growth(self, u):
        """
        This function simulates the potential growth of fruiting sites of cotton plants. It is called from PlantGrowth(). It calls TemperatureOnFruitGrowthRate()

        The following global variables are set here:
            PotGroAllBolls, PotGroAllBurrs, PotGroAllSquares.

        References:
        Marani, A. 1979. Growth rate of cotton bolls and their components. Field Crops Res. 2:169-175.
        Marani, A., Phene, C.J. and Cardon, G.E. 1992. CALGOS, a version of GOSSYM adapted for irrigated cotton.  III. leaf and boll growth routines. Beltwide Cotton Grow, Res. Conf. 1992:1361-1363.
        """
        global PotGroAllSquares, PotGroAllBolls, PotGroAllBurrs
        # The constant parameters used:
        cdef double[5] vpotfrt = [0.72, 0.30, 3.875, 0.125, 0.17]
        # Compute tfrt for the effect of temperature on boll and burr growth rates. Function TemperatureOnFruitGrowthRate() is used (with parameters derived from GOSSYM), for day time and night time temperatures, weighted by day and night lengths.
        cdef double tfrt  # the effect of temperature on rate of boll, burr or square growth.
        tfrt = (self._sim.states[u].day_length * TemperatureOnFruitGrowthRate(DayTimeTemp) + (24 - self._sim.states[u].day_length) * TemperatureOnFruitGrowthRate(NightTimeTemp)) / 24
        # Assign zero to sums of potential growth of squares, bolls and burrs.
        PotGroAllSquares = 0
        PotGroAllBolls = 0
        PotGroAllBurrs = 0
        # Assign values for the boll growth equation parameters. These are cultivar - specific.
        cdef double agemax = self.cultivar_parameters[9]  # maximum boll growth period (physiological days).
        cdef double rbmax = self.cultivar_parameters[10]  # maximum rate of boll (seed and lint) growth, g per boll per physiological day.
        cdef double wbmax = self.cultivar_parameters[11]  # maximum possible boll (seed and lint) weight, g per boll.
        # Loop for all vegetative stems.
        for k in range(self._sim.states[u].number_of_vegetative_branches):  # loop of vegetative stems
            for l in range(self._sim.states[u].vegetative_branches[k].number_of_fruiting_branches):  # loop of fruiting branches
                for m in range(self._sim.states[u].vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes):  # loop for nodes on a fruiting branch
                    # Calculate potential square growth for node (k,l,m).
                    # Sum potential growth rates of squares as PotGroAllSquares.
                    if self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].stage == Stage.Square:
                        # ratesqr is the rate of square growth, g per square per day.
                        # The routine for this is derived from GOSSYM, and so are the parameters used.
                        ratesqr = tfrt * vpotfrt[3] * exp(-vpotfrt[2] + vpotfrt[3] * self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].age)
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].square.potential_growth = ratesqr * self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].fraction
                        PotGroAllSquares += self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].square.potential_growth
                    # Growth of seedcotton is simulated separately from the growth of burrs. The logistic function is used to simulate growth of seedcotton. The constants of this function for cultivar 'Acala-SJ2', are based on the data of Marani (1979); they are derived from calibration for other cultivars
                    # agemax is the age of the boll (in physiological days after bloom) at the time when the boll growth rate is maximal.
                    # rbmax is the potential maximum rate of boll growth (g seeds plus lint dry weight per physiological day) at this age.
                    # wbmax is the maximum potential weight of seed plus lint (g dry weight per boll).
                    # The auxiliary variable pex is computed as
                    #    pex = exp(-4 * rbmax * (t - agemax) / wbmax)
                    # where t is the physiological age of the boll after bloom (= agebol).
                    # Boll weight (seed plus lint) at age T, according to the logistic function is:
                    #    wbol = wbmax / (1 + pex)
                    # and the potential boll growth rate at this age will be the derivative of this function:
                    #    ratebol = 4 * rbmax * pex / (1. + pex)**2
                    elif self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].stage == Stage.YoungGreenBoll or self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].stage == Stage.GreenBoll:
                        # pex is an intermediate variable to compute boll growth.
                        pex = exp(-4 * rbmax * (self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].boll.age - agemax) / wbmax)
                        # ratebol is the rate of boll (seed and lint) growth, g per boll per day.
                        ratebol = 4 * tfrt * rbmax * pex / (1 + pex) ** 2
                        # Potential growth rate of the burrs is assumed to be constant (vpotfrt[4] g dry weight per day) until the boll reaches its final volume. This occurs at the age of 22 physiological days in 'Acala-SJ2'. Both ratebol and ratebur are modified by temperature (tfrt) and ratebur is also affected by water stress (wfdb).
                        # Compute wfdb for the effect of water stress on burr growth rate. wfdb is the effect of water stress on rate of burr growth.
                        wfdb = vpotfrt[0] + vpotfrt[1] * self._sim.states[u].water_stress
                        if wfdb < 0:
                            wfdb = 0
                        if wfdb > 1:
                            wfdb = 1
                        ratebur = None  # rate of burr growth, g per boll per day.
                        if self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].boll.age >= 22:
                            ratebur = 0
                        else:
                            ratebur = vpotfrt[4] * tfrt * wfdb
                        # Potential boll (seeds and lint) growth rate (ratebol) and potential burr growth rate (ratebur) are multiplied by FruitFraction to compute PotGroBolls and PotGroBurrs for node (k,l,m).
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].boll.potential_growth = ratebol * self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].fraction
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].burr.potential_growth = ratebur * self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].fraction
                        # Sum potential growth rates of bolls and burrs as PotGroAllBolls and PotGroAllBurrs, respectively.
                        PotGroAllBolls += self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].boll.potential_growth
                        PotGroAllBurrs += self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].burr.potential_growth

                    # If these are not green bolls, their potential growth is 0. End loop.
                    else:
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].boll.potential_growth = 0
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].burr.potential_growth = 0

    def _potential_leaf_growth(self, u):
        """
        This function simulates the potential growth of leaves of cotton plants. It is called from self._growth(). It calls function temperature_on_leaf_growth_rate().

        The following monomolecular growth function is used :
            leaf area = smax * (1 - exp(-c * pow(t,p)))
        where    smax = maximum leaf area.
            t = time (leaf age).
            c, p = constant parameters.
        Note: p is constant for all leaves, whereas smax and c depend on leaf position.
        The rate per day (the derivative of this function) is :
            r = smax * c * p * exp(-c * pow(t,p)) * pow(t, (p-1))
        """
        global PotGroAllLeaves, PotGroAllPetioles
        # The following constant parameters. are used in this function:
        cdef double p = 1.6  # parameter of the leaf growth rate equation.
        cdef double[14] vpotlf = [3.0, 0.95, 1.2, 13.5, -0.62143, 0.109365, 0.00137566, 0.025, 0.00005, 30., 0.02, 0.001, 2.50, 0.18]
        # Calculate water stress reduction factor for leaf growth rate (wstrlf). This has been empirically calibrated in COTTON2K.
        cdef double wstrlf = self._sim.states[u].water_stress * (1 + vpotlf[0] * (2 - self._sim.states[u].water_stress)) - vpotlf[0]
        if wstrlf < 0.05:
            wstrlf = 0.05
        # Calculate wtfstrs, the effect of leaf water stress on state.leaf_weight_area_ratio (the ratio of leaf dry weight to leaf area). This has also been empirically calibrated in COTTON2K.
        cdef double wtfstrs = vpotlf[1] + vpotlf[2] * (1 - wstrlf)
        # Compute the ratio of leaf dry weight increment to leaf area increment (g per dm2), as a function of average daily temperature and water stress. Parameters for the effect of temperature are adapted from GOSSYM.
        cdef double tdday = self._sim.states[u].average_temperature  # limited value of today's average temperature.
        if tdday < vpotlf[3]:
            tdday = vpotlf[3]
        self._sim.states[u].leaf_weight_area_ratio = wtfstrs / (vpotlf[4] + tdday * (vpotlf[5] - tdday * vpotlf[6]))
        # Assign zero to total potential growth of leaf and petiole.
        PotGroAllLeaves = 0
        PotGroAllPetioles = 0
        cdef double c = 0  # parameter of the leaf growth rate equation.
        cdef double smax = 0  # maximum possible leaf area, a parameter of the leaf growth rate equation.
        cdef double rate  # growth rate of area of a leaf.
        # Compute the potential growth rate of prefruiting leaves. smax and c are functions of prefruiting node number.
        for j in range(self._sim.states[u].number_of_pre_fruiting_nodes):
            if self._sim.states[u].leaf_area_pre_fruiting[j] <= 0:
                PotGroLeafAreaPreFru[j] = 0
                PotGroLeafWeightPreFru[j] = 0
                PotGroPetioleWeightPreFru[j] = 0
            else:
                jp1 = j + 1
                smax = jp1 * (self.cultivar_parameters[2] - self.cultivar_parameters[3] * jp1)
                if smax < self.cultivar_parameters[4]:
                    smax = self.cultivar_parameters[4]
                c = vpotlf[7] + vpotlf[8] * jp1 * (jp1 - vpotlf[9])
                rate = smax * c * p * exp(-c * pow(self._sim.states[u].age_of_pre_fruiting_nodes[j], p)) * pow(self._sim.states[u].age_of_pre_fruiting_nodes[j], (p - 1))
                # Growth rate is modified by water stress and a function of average temperature.
                # Compute potential growth of leaf area, leaf weight and petiole weight for leaf on node j. Add leaf weight potential growth to PotGroAllLeaves.
                # Add potential growth of petiole weight to PotGroAllPetioles.
                if rate >= 1e-12:
                    PotGroLeafAreaPreFru[j] = rate * wstrlf * temperature_on_leaf_growth_rate(self._sim.states[u].average_temperature)
                    PotGroLeafWeightPreFru[j] = PotGroLeafAreaPreFru[j] * self._sim.states[u].leaf_weight_area_ratio
                    PotGroPetioleWeightPreFru[j] = PotGroLeafAreaPreFru[j] * self._sim.states[u].leaf_weight_area_ratio * vpotlf[13]
                    PotGroAllLeaves += PotGroLeafWeightPreFru[j]
                    PotGroAllPetioles += PotGroPetioleWeightPreFru[j]
        # denfac is the effect of plant density on leaf growth rate.
        cdef double denfac = 1 - vpotlf[12] * (1 - self._sim.density_factor)
        for k in range(self._sim.states[u].number_of_vegetative_branches):
            for l in range(self._sim.states[u].vegetative_branches[k].number_of_fruiting_branches):
                # smax and c are  functions of fruiting branch number.
                # smax is modified by plant density, using the density factor denfac.
                # Compute potential main stem leaf growth, assuming that the main stem leaf is initiated at the same time as leaf (k,l,0).
                if self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.leaf_area <= 0:
                    self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_leaf_area = 0
                    self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_leaf_weight = 0
                    self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_petiole_weight = 0
                else:
                    lp1 = l + 1
                    smax = self.cultivar_parameters[5] + self.cultivar_parameters[6] * lp1 * (self.cultivar_parameters[7] - lp1)
                    smax = smax * denfac
                    if smax < self.cultivar_parameters[4]:
                        smax = self.cultivar_parameters[4]
                    c = vpotlf[10] + lp1 * vpotlf[11]
                    if self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[0].leaf.age > 70:
                        rate = 0
                    else:
                        rate = smax * c * p * exp(-c * pow(self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[0].leaf.age, p)) * pow(self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[0].leaf.age, (p - 1))
                    # Add leaf and petiole weight potential growth to SPDWL and SPDWP.
                    if rate >= 1e-12:
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_leaf_area = rate * wstrlf * temperature_on_leaf_growth_rate(self._sim.states[u].average_temperature)
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_leaf_weight = self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_leaf_area * self._sim.states[u].leaf_weight_area_ratio
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_petiole_weight = self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_leaf_area * self._sim.states[u].leaf_weight_area_ratio * vpotlf[13]
                        PotGroAllLeaves += self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_leaf_weight
                        PotGroAllPetioles += self._sim.states[u].vegetative_branches[k].fruiting_branches[l].main_stem_leaf.potential_growth_for_petiole_weight
                # Assign smax value of this main stem leaf to smaxx, c to cc.
                # Loop over the nodes of this fruiting branch.
                smaxx = smax  # value of smax for the corresponding main stem leaf.
                cc = c  # value of c for the corresponding main stem leaf.
                for m in range(self._sim.states[u].vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes):
                    if self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.area <= 0:
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.potential_growth = 0
                        self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].petiole.potential_growth = 0
                    # Compute potential growth of leaf area and leaf weight for leaf on fruiting branch node (k,l,m).
                    # Add leaf and petiole weight potential growth to spdwl and spdwp.
                    else:
                        mp1 = m + 1
                        # smax and c are reduced as a function of node number on this fruiting branch.
                        smax = smaxx * (1 - self.cultivar_parameters[8] * mp1)
                        c = cc * (1 - self.cultivar_parameters[8] * mp1)
                        # Compute potential growth for the leaves on fruiting branches.
                        if self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.age > 70:
                            rate = 0
                        else:
                            rate = smax * c * p * exp(-c * pow(self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.age, p)) * pow(self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.age, (p - 1))
                        if rate >= 1e-12:
                            # Growth rate is modified by water stress. Potential growth is computed as a function of average temperature.
                            self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.potential_growth = rate * wstrlf * temperature_on_leaf_growth_rate(self._sim.states[u].average_temperature)
                            self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].petiole.potential_growth = self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.potential_growth * self._sim.states[u].leaf_weight_area_ratio * vpotlf[13]
                            PotGroAllLeaves += self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.potential_growth * self._sim.states[u].leaf_weight_area_ratio
                            PotGroAllPetioles += self._sim.states[u].vegetative_branches[k].fruiting_branches[l].nodes[m].petiole.potential_growth

    def _growth(self, u):
        state = self.state(u)
        global PotGroStem, PotGroAllRoots
        # Call self._potential_leaf_growth(u) to compute potential growth rate of leaves.
        self._potential_leaf_growth(u)
        # If it is after first square, call self._potential_fruit_growth(u) to compute potential growth rate of squares and bolls.
        if self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].stage != Stage.NotYetFormed:
            self._potential_fruit_growth(u)
        # Active stem tissue(stemnew) is the difference between state.stem_weight and the value of StemWeight(kkday).
        cdef int voldstm = 32  # constant parameter(days for stem tissue to become "old")
        cdef int kkday = self._sim.states[u].kday - voldstm  # age of young stem tissue
        if kkday < 1:
            kkday = 1
        cdef double stemnew = self._sim.states[u].stem_weight - StemWeight[kkday]  # dry weight of active stem tissue.
        # Call PotentialStemGrowth() to compute PotGroStem, potential growth rate of stems.
        # The effect of temperature is introduced, by multiplying potential growth rate by DayInc.
        # Stem growth is also affected by water stress(WaterStressStem).PotGroStem is limited by (maxstmgr * per_plant_area) g per plant per day.
        PotGroStem = PotentialStemGrowth(stemnew, self._sim.states[u].kday,
                                     self._sim.states[u].vegetative_branches[0].fruiting_branches[2].nodes[0].stage, self._sim.density_factor,
                                     self.cultivar_parameters[12], self.cultivar_parameters[13], self.cultivar_parameters[14],
                                     self.cultivar_parameters[15], self.cultivar_parameters[16], self.cultivar_parameters[17],
                                     self.cultivar_parameters[18]) * self._sim.states[u].day_inc * self._sim.states[u].water_stress_stem
        cdef double maxstmgr = 0.067  # maximum posible potential stem growth, g dm - 2 day - 1.
        if PotGroStem > maxstmgr * self._sim.per_plant_area:
            PotGroStem = maxstmgr * self._sim.per_plant_area
        # Call PotentialRootGrowth() to compute potential growth rate on roots.
        cdef double sumpdr  # total potential growth rate of roots in g per slab.this is computed in PotentialRootGrowth() and used in ActualRootGrowth().
        sumpdr = State.from_ptr(&self._sim.states[u], self.year, self.version).potential_root_growth(3, self._sim.states[u].soil.number_of_layers_with_root, self._sim.per_plant_area)
        # Total potential growth rate of roots is converted from g per slab(sumpdr) to g per plant(PotGroAllRoots).
        PotGroAllRoots = sumpdr * 100 * self._sim.per_plant_area / self.row_space
        # Limit PotGroAllRoots to(maxrtgr * per_plant_area) g per plant per day.
        cdef double maxrtgr = 0.045  # maximum possible potential root growth, g dm - 2 day - 1.
        if PotGroAllRoots > maxrtgr * self._sim.per_plant_area:
            PotGroAllRoots = maxrtgr * self._sim.per_plant_area
        # Call DryMatterBalance() to compute carbon balance, allocation of carbon to plant parts, and carbon stress.DryMatterBalance() also computes and returns the values of the following arguments:
        # cdleaf is carbohydrate requirement for leaf growth, g per plant per day.
        # cdpet is carbohydrate requirement for petiole growth, g per plant per day.
        # cdroot is carbohydrate requirement for root growth, g per plant per day.
        # cdstem is carbohydrate requirement for stem growth, g per plant per day.
        cdef double cdstem, cdleaf, cdpet, cdroot
        DryMatterBalance(self._sim.states[u], cdstem, cdleaf, cdpet, cdroot, self._sim.per_plant_area)
        # If it is after first square, call ActualFruitGrowth() to compute actual
        # growth rate of squares and bolls.
        if self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].stage != Stage.NotYetFormed:
            ActualFruitGrowth(self._sim.states[u])
        # Initialize state.leaf_weight.It is assumed that cotyledons fall off at time of first square.Also initialize state.leaf_area and state.petiole_weight.
        if self._sim.first_square > 0:
            self._sim.states[u].leaf_weight = 0
            self._sim.states[u].leaf_area = 0
        else:
            cotylwt = 0.20  # weight of cotyledons dry matter.
            self._sim.states[u].leaf_weight = cotylwt
            self._sim.states[u].leaf_area = 0.6 * cotylwt
        self._sim.states[u].petiole_weight = 0
        # Call ActualLeafGrowth to compute actual growth rate of leaves and compute leaf area index.
        ActualLeafGrowth(self._sim.states[u])
        self._sim.states[u].leaf_area_index = self._sim.states[u].leaf_area / self._sim.per_plant_area
        # Add ActualStemGrowth to state.stem_weight, and define StemWeight(Kday) for this day.
        self._sim.states[u].stem_weight += ActualStemGrowth
        StemWeight[self._sim.states[u].kday] = self._sim.states[u].stem_weight
        # Plant density affects growth in height of tall plants.
        cdef double htdenf = 55  # minimum plant height for plant density affecting growth in height.
        cdef double z1  # intermediate variable to compute denf2.
        z1 = (self._sim.states[u].plant_height - htdenf) / htdenf
        if z1 < 0:
            z1 = 0
        if z1 > 1:
            z1 = 1
        cdef double denf2  # effect of plant density on plant growth in height.
        denf2 = 1 + z1 * (self._sim.density_factor - 1)
        # Call AddPlantHeight to compute PlantHeight.
        cdef int l, l1, l2  # node numbers of top three nodes.
        l = self._sim.states[u].vegetative_branches[0].number_of_fruiting_branches - 1
        l1 = l - 1
        if l < 1:
            l1 = 0
        l2 = l - 2
        if l < 2:
            l2 = 0
        cdef double agetop  # average physiological age of top three nodes.
        agetop = (self._sim.states[u].vegetative_branches[0].fruiting_branches[l].nodes[0].age +
                  self._sim.states[u].vegetative_branches[0].fruiting_branches[l1].nodes[0].age +
                  self._sim.states[u].vegetative_branches[0].fruiting_branches[l2].nodes[0].age) / 3
        if self.version < 0x500 or self._sim.day_topping <= 0 or self._sim.states[u].daynum < self._sim.day_topping:
            self._sim.states[u].plant_height += AddPlantHeight(denf2, self._sim.states[u].day_inc, self._sim.states[u].number_of_pre_fruiting_nodes,
                                                 self._sim.states[u].vegetative_branches[0].fruiting_branches[1].nodes[0].stage,
                                                 self._sim.states[u].age_of_pre_fruiting_nodes[
                                                     self._sim.states[u].number_of_pre_fruiting_nodes - 1],
                                                 self._sim.states[u].age_of_pre_fruiting_nodes[
                                                     self._sim.states[u].number_of_pre_fruiting_nodes - 2], agetop,
                                                 self._sim.states[u].water_stress_stem, self._sim.states[u].carbon_stress,
                                                 self._sim.states[u].nitrogen_stress_vegetative, self.cultivar_parameters[19],
                                                 self.cultivar_parameters[20], self.cultivar_parameters[21],
                                                 self.cultivar_parameters[22], self.cultivar_parameters[23],
                                                 self.cultivar_parameters[24], self.cultivar_parameters[25],
                                                 self.cultivar_parameters[26])
        # Call ActualRootGrowth() to compute actual root growth.
        state.compute_actual_root_growth(sumpdr, self.row_space, self._sim.per_plant_area, 3, self._sim.day_emerge, self._sim.plant_row_column)

    def _add_vegetative_branch(self, u, delayVegByCStress, stemNRatio, DaysTo1stSqare, PhenDelayByNStress):
        """
        This function decides whether a new vegetative branch is to be added, and then forms it. It is called from CottonPhenology().
        """
        if self._sim.states[u].number_of_vegetative_branches == 3:
            return
        # The following constant parameters are used:
        cdef double[3] vpvegb = [13.39, -0.696, 0.012]
        # TimeToNextVegBranch is computed as a function of this average temperature.
        cdef double TimeToNextVegBranch  # time, in physiological days, for the next vegetative branch to be formed.
        TimeToNextVegBranch = vpvegb[0] + self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].average_temperature * (vpvegb[1] + self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].average_temperature * vpvegb[2])
        # Compare the age of the first fruiting site of the last formed vegetative branch with TimeToNextVegBranch plus DaysTo1stSqare and the delays caused by stresses, in order to decide if a new vegetative branch is to be formed.
        if self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].age < TimeToNextVegBranch + delayVegByCStress + PhenDelayByNStress + DaysTo1stSqare:
            return
        # When a new vegetative branch is formed, increase NumVegBranches by 1.
        self._sim.states[u].number_of_vegetative_branches += 1
        # Assign 1 to FruitFraction and FruitingCode of the first site of this branch.
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].fraction = 1
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].stage = Stage.Square
        # Add a new leaf to the first site of this branch.
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].leaf.area = self.cultivar_parameters[34]
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].leaf.weight = self.cultivar_parameters[34] * self._sim.states[u].leaf_weight_area_ratio
        # Add a new mainstem leaf to the first node of this branch.
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].main_stem_leaf.leaf_area = self.cultivar_parameters[34]
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].main_stem_leaf.leaf_weight = self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].main_stem_leaf.leaf_area * self._sim.states[u].leaf_weight_area_ratio
        # The initial mass and nitrogen in the new leaves are substracted from the stem.
        self._sim.states[u].stem_weight -= self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].leaf.weight + self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].main_stem_leaf.leaf_weight
        self._sim.states[u].leaf_weight += self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].leaf.weight + self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].main_stem_leaf.leaf_weight
        cdef double addlfn  # nitrogen moved to new leaves from stem.
        addlfn = (self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].leaf.weight + self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].main_stem_leaf.leaf_weight) * stemNRatio
        self._sim.states[u].leaf_nitrogen += addlfn
        self._sim.states[u].stem_nitrogen -= addlfn
        # Assign the initial value of the average temperature of the first site.
        # Define initial NumFruitBranches and NumNodes for the new vegetative branch.
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].nodes[0].average_temperature = self._sim.states[u].average_temperature
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].number_of_fruiting_branches = 1
        self._sim.states[u].vegetative_branches[self._sim.states[u].number_of_vegetative_branches - 1].fruiting_branches[0].number_of_fruiting_nodes = 1

    def _phenology(self, u):
        """
        This is is the main function for simulating events of phenology and abscission in the cotton plant. It is called each day from DailySimulation().

        CottonPhenology() calls PreFruitingNode(), days_to_first_square(), create_first_square(), _add_vegetative_branch(), _add_fruiting_branch(), AddFruitingNode(), SimulateFruitingSite(), LeafAbscission(), FruitingSitesAbscission().
        """
        state = self.state(u)
        # The following constant parameters are used:
        cdef double[8] vpheno = [0.65, -0.83, -1.67, -0.25, -0.75, 10.0, 15.0, 7.10]

        cdef int nwfl = 0  # the node of the most recent white flower. Note: this variable is not used. It is kept for compatibility with previous versions, and may be use in future versions.
        self._sim.states[u].number_of_fruiting_sites = 0
        cdef double stemNRatio  # the ratio of N to dry matter in the stems.
        stemNRatio = self._sim.states[u].stem_nitrogen / self._sim.states[u].stem_weight
        # Compute the phenological delays:
        # PhenDelayByNStress, the delay caused by nitrogen stress, is assumed to be a function of the vegetative nitrogen stress.
        cdef double PhenDelayByNStress = vpheno[0] * (1 - self._sim.states[u].nitrogen_stress_vegetative)  # phenological delay caused by vegetative nitrogen stress.
        if PhenDelayByNStress > 1:
            PhenDelayByNStress = 1
        elif PhenDelayByNStress < 0:
            PhenDelayByNStress = 0

        cdef double delayVegByCStress  # delay in formation of new fruiting branches caused by carbon stress.
        delayVegByCStress = self.cultivar_parameters[27] + self._sim.states[u].carbon_stress * (vpheno[3] + vpheno[4] * self._sim.states[u].carbon_stress)
        if delayVegByCStress > 1:
            delayVegByCStress = 1
        elif delayVegByCStress < 0:
            delayVegByCStress = 0

        cdef double delayFrtByCStress  # delay in formation of new fruiting sites caused by carbon stress.
        delayFrtByCStress = self.cultivar_parameters[28] + self._sim.states[u].carbon_stress * (vpheno[1] + vpheno[2] * self._sim.states[u].carbon_stress)
        if delayFrtByCStress > self.cultivar_parameters[29]:
            delayFrtByCStress = self.cultivar_parameters[29]
        if delayFrtByCStress < 0:
            delayFrtByCStress = 0
        # The following section is executed if the first square has not yet been formed. Function days_to_first_square() is called to compute the number of days to 1st square, and function PreFruitingNode() is called to simulate the formation of prefruiting nodes.
        if self._sim.first_square <= 0:
            self.DaysTo1stSqare = days_to_first_square(state.average_temperature, state.water_stress, state.nitrogen_stress_vegetative, self.cultivar_parameters[30])
            PreFruitingNode(self._sim.states[u], stemNRatio, self._sim.cultivar_parameters)
            # When first square is formed, FirstSquare is assigned the day of year.
            # Function create_first_square() is called for formation of first square.
            if state.kday >= <int>self.DaysTo1stSqare:
                self._sim.first_square = self._sim.states[u].daynum
                state.create_first_square(stemNRatio, self.cultivar_parameters[34])
            # if a first square has not been formed, call LeafAbscission() and exit.
            else:
                state.leaf_abscission(self._sim.per_plant_area, self.first_square_date, self.defoliate_date)
                return
        # The following is executed after the appearance of the first square.
        # If there are only one or two vegetative branches, and if plant population allows it, call _add_vegetative_branch() to decide if a new vegetative branch is to be added. Note that dense plant populations (large per_plant_area) prevent new vegetative branch formation.
        if self._sim.states[u].number_of_vegetative_branches == 1 and self._sim.per_plant_area >= vpheno[5]:
            self._add_vegetative_branch(u, delayVegByCStress, stemNRatio, self.DaysTo1stSqare, PhenDelayByNStress)
        if self._sim.states[u].number_of_vegetative_branches == 2 and self._sim.per_plant_area >= vpheno[6]:
            self._add_vegetative_branch(u, delayVegByCStress, stemNRatio, self.DaysTo1stSqare, PhenDelayByNStress)
        # The maximum number of nodes per fruiting branch (nidmax) is affected by plant density. It is computed as a function of density_factor.
        cdef int nidmax  # maximum number of nodes per fruiting branch.
        nidmax = <int>(vpheno[7] * self._sim.density_factor + 0.5)
        if nidmax > 5:
            nidmax = 5
        # Start loop over all existing vegetative branches.
        # Call AddFruitingBranch() to decide if a new node (and a new fruiting branch) is to be added on this stem.
        for k in range(self._sim.states[u].number_of_vegetative_branches):
            if self._sim.states[u].vegetative_branches[k].number_of_fruiting_branches < 30:
                self.state(u).add_fruiting_branch(k, self._sim.density_factor, delayVegByCStress, delayFrtByCStress, stemNRatio, PhenDelayByNStress, self.cultivar_parameters[35], self.cultivar_parameters[34], self.topping_date)
            # Loop over all existing fruiting branches, and call AddFruitingNode() to decide if a new node on this fruiting branch is to be added.
            for l in range(self._sim.states[u].vegetative_branches[k].number_of_fruiting_branches):
                if self._sim.states[u].vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes < nidmax:
                    AddFruitingNode(self._sim.states[u], k, l, delayFrtByCStress, stemNRatio, self._sim.density_factor, self._sim.cultivar_parameters, PhenDelayByNStress)
                # Loop over all existing fruiting nodes, and call SimulateFruitingSite() to simulate the condition of each fruiting node.
                for m in range(self._sim.states[u].vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes):
                    SimulateFruitingSite(self._sim, u, k, l, m, nwfl, self._sim.states[u].water_stress)
        # Call FruitingSitesAbscission() to simulate the abscission of fruiting parts.
        FruitingSitesAbscission(self._sim, u)
        # Call LeafAbscission() to simulate the abscission of leaves.
        state.leaf_abscission(self._sim.per_plant_area, self.first_square_date, self.defoliate_date)

    def _defoliate(self, u):
        """This function simulates the effects of defoliating chemicals applied on the cotton. It is called from SimulateThisDay()."""
        global PercentDefoliation
        cdef cState state = self._sim.states[u]
        # constant parameters:
        cdef double p1 = -50.0
        cdef double p2 = 0.525
        cdef double p3 = 7.06
        cdef double p4 = 0.85
        cdef double p5 = 2.48
        cdef double p6 = 0.0374
        cdef double p7 = 0.0020

        # If this is first day set initial values of tdfkgh, defkgh to 0.
        if state.daynum <= self._sim.day_emerge:
            self.tdfkgh = 0
            self.defkgh = 0
            self.idsw = False
        # Start a loop for five possible defoliant applications.
        for i in range(5):
            # If there are open bolls and defoliation prediction has been set, execute the following.
            if state.number_of_open_bolls > 0 and DefoliantAppRate[i] <= -99.9:
                # percentage of open bolls in total boll number
                OpenRatio = <int>(100 * state.number_of_open_bolls / (state.number_of_open_bolls + state.number_of_green_bolls))
                if i == 0 and not self.idsw:
                    # If this is first defoliation - check the percentage of boll opening.
                    # If it is after the defined date, or the percent boll opening is greater than the defined threshold - set defoliation date as this day and set a second prediction.
                    if (state.daynum >= DefoliationDate[0] and DefoliationDate[0] > 0) or OpenRatio > DefoliationMethod[i]:
                        self.idsw = True
                        DefoliationDate[0] = state.daynum
                        DefoliantAppRate[1] = -99.9
                        if state.daynum < self._sim.day_defoliate or self._sim.day_defoliate <= 0:
                            self._sim.day_defoliate = state.daynum
                        DefoliationMethod[0] = 0
                # If 10 days have passed since the last defoliation, and the leaf area index is still greater than 0.2, set another defoliation.
                if i >= 1:
                    if state.daynum == (DefoliationDate[i - 1] + 10) and state.leaf_area_index >= 0.2:
                        DefoliationDate[i] = state.daynum
                        if i < 4:
                            DefoliantAppRate[i + 1] = -99.9
                        DefoliationMethod[i] = 0
            if state.daynum == DefoliationDate[i]:
                # If it is a predicted defoliation, assign tdfkgh as 2.5 .
                # Else, compute the amount intercepted by the plants in kg per ha (defkgh), and add it to tdfkgh.
                if DefoliantAppRate[i] < -99:
                    self.tdfkgh = 2.5
                else:
                    if DefoliationMethod[i] == 0:
                        self.defkgh += DefoliantAppRate[i] * 0.95 * 1.12085 * 0.75
                    else:
                        self.defkgh += DefoliantAppRate[i] * state.light_interception * 1.12085 * 0.75
                    self.tdfkgh += self.defkgh
            # If this is after the first day of defoliant application, compute the percent of leaves to be defoliated (PercentDefoliation), as a function of average daily temperature, leaf water potential, days after first defoliation application, and tdfkgh. The regression equation is modified from the equation suggested in GOSSYM.
            if DefoliationDate[i] > 0 and state.daynum > self._sim.day_defoliate:
                dum = -LwpMin * 10  # value of LwpMin in bars.
                PercentDefoliation = p1 + p2 * state.average_temperature + p3 * self.tdfkgh + p4 * (state.daynum - self._sim.day_defoliate) + p5 * dum - p6 * dum * dum + p7 * state.average_temperature * self.tdfkgh * (state.daynum - self._sim.day_defoliate) * dum
                if PercentDefoliation < 0:
                    PercentDefoliation = 0
                perdmax = 40  # maximum possible percent of defoliation.
                if PercentDefoliation > perdmax:
                    PercentDefoliation = perdmax


    def _simulate_this_day(self, u):
        global isw
        state = self.state(u)
        if 0 < self._sim.day_emerge <= self._sim.day_start + u:
            state.kday = (self._sim.day_start + u) - self._sim.day_emerge + 1
            if state.leaf_area_index > self.max_leaf_area_index:
                self.max_leaf_area_index = state.leaf_area_index
            state.light_interception = compute_light_interception(state.leaf_area_index, self.max_leaf_area_index, state.plant_height, self.row_space, version=self.version)
            self._column_shading(u)
        else:
            state.kday = 0
            state.light_interception = 0
            self.relative_radiation_received_by_a_soil_column = [1] * 20
        # The following functions are executed each day (also before emergence).
        self._daily_climate(u)  # computes climate variables for today.
        self._soil_temperature(u)  # executes all modules of soil and canopy temperature.
        SoilProcedures(self._sim, u)  # executes all other soil processes.
        SoilNitrogen(self._sim, u)  # computes nitrogen transformations in the soil.
        SoilSum(self._sim.states[u], self._sim.row_space)  # computes totals of water and N in the soil.
        # The following is executed each day after plant emergence:
        if state.daynum >= self._sim.day_emerge and isw > 0:
            # If this day is after emergence, assign to isw the value of 2.
            isw = 2
            state.day_inc = PhysiologicalAge(self._sim.states[u].hours)  # physiological days increment for this day. computes physiological age
            self._defoliate(u)  # effects of defoliants applied.
            self._stress(u)  # computes water stress factors.
            self._get_net_photosynthesis(u)  # computes net photosynthesis.
            self._growth(u)  # executes all modules of plant growth.
            self._phenology(u)  # executes all modules of plant phenology.
            PlantNitrogen(self._sim, u)  # computes plant nitrogen allocation.
            CheckDryMatterBal(self._sim.states[u])  # checks plant dry matter balance.
        # Check if the date to stop simulation has been reached, or if this is the last day with available weather data. Simulation will also stop when no leaves remain on the plant.
        if state.daynum >= LastDayWeatherData or (state.kday > 10 and state.leaf_area_index < 0.0002):
            raise SimulationEnd

    def _init_grid(self):
        """
        This function initializes the soil grid variables. It is executed once at the beginning of the simulation. It is called from ReadInput().

        The following global or file-scope variables are set here:
        nk, nl.

        The following global variables are referenced here:
        maxk, maxl.
        """
        # PlantRowLocation is the distance from edge of slab, cm, of the plant row.
        global PlantRowLocation, nl, nk
        PlantRowLocation = 0.5 * self.row_space
        if self.skip_row_width > 1:
            # If there is a skiprow arrangement, RowSpace and PlantRowLocation are redefined.
            self.row_space = 0.5 * (
                    self.row_space + self.skip_row_width)  # actual width of the soil slab (cm)
            PlantRowLocation = 0.5 * self.skip_row_width
        # Compute sim.plant_population - number of plants per hectar, and per_plant_area - the average surface area per plant, in dm2, and the empirical plant density factor (density_factor). This factor will be used to express the effect of plant density on some plant growth rate functions.
        # NOTE: density_factor = 1 for 5 plants per sq m (or 50000 per ha).
        self._sim.plant_population = self.plants_per_meter / self.row_space * 1000000
        self._sim.per_plant_area = 1000000 / self._sim.plant_population
        self._sim.density_factor = exp(self.cultivar_parameters[1] * (5 - self._sim.plant_population / 10000))
        # Define the numbers of rows and columns in the soil slab (nl, nk).
        # Define the depth, in cm, of consecutive nl layers.
        # NOTE: maxl and maxk are defined as constants in file "global.h".
        nl = maxl
        nk = maxk
        # The width of the slab columns is computed by dividing the row spacing by the number of columns. It is assumed that slab width is equal to the average row spacing, and column widths are uniform.
        # NOTE: wk is an array - to enable the option of non-uniform column widths in the future.
        # PlantRowColumn (the column including the plant row) is now computed from PlantRowLocation (the distance of the plant row from the edge of the slab).
        cdef double sumwk = 0  # sum of column widths
        self._sim.plant_row_column = 0
        for k in range(nk):
            sumwk = sumwk + wk(k, self.row_space)
            if self._sim.plant_row_column == 0 and sumwk > PlantRowLocation:
                if (sumwk - PlantRowLocation) > (0.5 * wk(k, self.row_space)):
                    self._sim.plant_row_column = k - 1
                else:
                    self._sim.plant_row_column = k

    def read_input(self, lyrsol, **kwargs):
        """This is the main function for reading input."""
        InitializeGlobal()
        initialize_switch(self._sim)
        self._init_grid()
        read_agricultural_input(self._sim, kwargs.get("agricultural_inputs", []))
        InitializeSoilData(self._sim, lyrsol)
        InitializeSoilTemperature()
        InitializeRootData(self._sim)
