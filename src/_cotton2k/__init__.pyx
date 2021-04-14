# distutils: language=c++
# cython: language_level=3
from libc.stdlib cimport malloc
from libc.math cimport exp

from datetime import datetime, date

from _cotton2k._global cimport *
from _cotton2k._structs cimport *
from _cotton2k._io cimport *

cdef extern from "Cottonmodel.h":

    cdef cppclass C2KApp:
        C2KApp() except +
        void DailySimulation(cSimulation &)


def _date2doy(d):
    if isinstance(d, str):
        d = datetime.strptime(d, "%Y-%m-%d")
    if isinstance(d, date):
        return d.timetuple().tm_yday
    elif isinstance(d, int) and d > 0:
        return d
    else:
        return 0

cdef void initialize_switch(cSimulation &sim):
    global isw, Kday
    # If the date of emergence has not been given, emergence will be simulated
    # by the model. In this case, isw = 0, and a check is performed to make
    # sure that the date of planting has been given.
    if sim.day_emerge <= 0:
        if sim.day_plant <= 0:
            raise Exception(" planting date or emergence date must be given in the profile file !!")
        isw = 0
    # If the date of emergence has been given in the input: isw = 1 if
    # simulation starts before emergence, or isw = 2 if simulation starts at emergence.
    elif sim.day_emerge > sim.day_start:
        isw = 1
    else:
        isw = 2
        Kday = 1

cdef double SkipRowWidth  # the smaller distance between skip rows, cm
cdef double PlantsPerM  # average number of plants pre meter of row.

cdef void InitializeGrid(cSimulation &sim):
    """
    This function initializes the soil grid variables. It is executed once at the beginning of the simulation. It is called from ReadInput().

    The following global or file-scope variables are set here:
    DensityFactor, dl, nk, nl, PerPlantArea, PlantPopulation,
    PlantRowColumn, PlantRowLocation, RowSpace, wk.

    The following global variables are referenced here:
    PlantsPerM, SkipRowWidth, VarPar, maxk, maxl."""
    # PlantRowLocation is the distance from edge of slab, cm, of the plant row.
    global PlantRowLocation, PlantPopulation, PerPlantArea, DensityFactor, nl, nk, SkipRowWidth, PlantsPerM
    PlantRowLocation = 0.5 * sim.row_space
    if (SkipRowWidth > 1):
        # If there is a skiprow arrangement, RowSpace and PlantRowLocation are redefined.
        sim.row_space = 0.5 * (sim.row_space + SkipRowWidth)  # actual width of the soil slab (cm)
        PlantRowLocation = 0.5 * SkipRowWidth
    # Compute PlantPopulation - number of plants per hectar, and PerPlantArea - the average surface area per plant, in dm2, and the empirical plant density factor (DensityFactor). This factor will be used to express the effect of plant density on some plant growth rate functions.
    # NOTE: DensityFactor = 1 for 5 plants per sq m (or 50000 per ha).
    PlantPopulation = PlantsPerM / sim.row_space * 1000000
    PerPlantArea = 1000000 / PlantPopulation
    DensityFactor = exp(VarPar[1] * (5 - PlantPopulation / 10000))
    # Define the numbers of rows and columns in the soil slab (nl, nk).
    # Define the depth, in cm, of consecutive nl layers.
    # NOTE: maxl and maxk are defined as constants in file "global.h".
    nl = maxl
    nk = maxk
    # The width of the slab columns is computed by dividing the row spacing by the number of columns. It is assumed that slab width is equal to the average row spacing, and column widths are uniform.
    # NOTE: wk is an array - to enable the option of non-uniform column widths in the future.
    # PlantRowColumn (the column including the plant row) is now computed from PlantRowLocation (the distance of the plant row from the edge of the slab).
    cdef double sumwk = 0  # sum of column widths
    sim.plant_row_column = 0
    for k in range(nk):
        sumwk = sumwk + wk(k, sim.row_space)
        if sim.plant_row_column == 0 and sumwk > PlantRowLocation:
            if (sumwk - PlantRowLocation) > (0.5 * wk(k, sim.row_space)):
                sim.plant_row_column = k - 1
            else:
                sim.plant_row_column = k


cdef class Soil:
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



cdef class SoilImpedance:

    @property
    def curves(self):
        return {gh2oc[i]: {tstbd[j][i]: impede[j][i] for j in range(inrim)} for i in range(ncurve)}

    @curves.setter
    def curves(self, impedance_table):
        global ncurve, inrim
        ncurve = len(impedance_table)
        inrim = len(impedance_table[0])
        for i, row in enumerate(impedance_table):
            gh2oc[i] = row.pop("water")
            for j, pair in enumerate(sorted(row.items())):
                tstbd[j][i], impede[j][i] = pair


cdef class Climate:
    cdef ClimateStruct *climate
    cdef unsigned int start_day
    cdef unsigned int days
    cdef unsigned int current

    def __init__(self, start_date, climate):
        global LastDayWeatherData
        self.start_day = _date2doy(start_date)
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
            self.climate[i].Tdew = daily_climate.get("dewpoint", tdewest(daily_climate["max"], SitePar[5], SitePar[6]))

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start or 0
            stop = key.stop or self.days
            step = key.step or 1
            if isinstance(start, date):
                start = _date2doy(start) - self.start_day
            if isinstance(stop, date):
                stop = _date2doy(stop) - self.start_day
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
                key = _date2doy(key) - self.start_day
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



cdef read_agricultural_input(cSimulation &sim, inputs):
    global NumNitApps, NumIrrigations
    NumNitApps = 0
    idef = 0
    cdef Irrigation irrigation
    cdef NitrogenFertilizer nf
    for i in inputs:
        if i["type"] == "irrigation":
            irrigation.day = _date2doy(i["date"])  # day of year of this irrigation
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
            nf.day = _date2doy(i["date"])
            nf.amtamm = i.get("ammonium", 0)
            nf.amtnit = i.get("nitrate", 0)
            nf.amtura = i.get("urea", 0)
            nf.mthfrt = i.get("method", 0)
            isdhrz = i.get("drip_horizontal_place", 0)  # horizontal placement of DRIP, cm from left edge of soil slab.
            isddph = i.get("drip_depth", 0)  # vertical placement of DRIP, cm from soil surface.
            if nf.mthfrt == 1 or nf.mthfrt == 3:
                nf.ksdr = SlabLoc(isdhrz, sim.row_space)
                nf.lsdr = SlabLoc(isddph, 0)
            else:
                nf.ksdr = 0
                nf.lsdr = 0
            NFertilizer[NumNitApps] = nf
            NumNitApps += 1
        elif i["type"] == "defoliation prediction":
            DefoliationDate[idef] = _date2doy(i["date"])
            DefoliantAppRate[idef] = -99.9
            if idef == 0:
                DayFirstDef = DefoliationDate[0]
            DefoliationMethod[idef] = i.get("method", 0)
            idef += 1

cdef class FruitingBranch:
    cdef cFruitingBranch _branch
    __slots__ = ("delay_for_new_node", "main_stem_leaf", "nodes")

    def __init__(self, _branch):
        self._branch = _branch

    @property
    def delay_for_new_node(self):
        return self._branch.delay_for_new_node

    @property
    def main_stem_leaf(self):
        return self._branch.main_stem_leaf

    @property
    def nodes(self):
        return [self._branch.nodes[i] for i in range(self._branch.number_of_fruiting_nodes)]

    def __iter__(self):
        for attr in self.__slots__:
            yield attr, getattr(self, attr)


cdef class VegetativeBranch:
    cdef cVegetativeBranch _branch

    def __init__(self, _branch):
        self._branch = _branch

    @property
    def fruiting_branches(self):
        return [FruitingBranch(self._branch.fruiting_branches[i]) for i in range(self._branch.number_of_fruiting_branches)]

    def __iter__(self):
        yield "fruiting_branches", self.fruiting_branches


cdef class State:
    cdef cState _state
    __slots__ = (
        "plant_height",
        "plant_weight",
        "lint_yield",
        "number_of_squares",
        "number_of_green_bolls",
        "number_of_open_bolls",
        "leaf_area_index",
        "ginning_percent",
        "vegetative_branches",
        "hours",
        "soil",
    )

    def __init__(self, _state):
        self._state = _state

    @property
    def plant_height(self):
        return self._state.plant_height

    @plant_height.setter
    def plant_height(self, value):
        self._state.plant_height = value

    @property
    def plant_weight(self):
        return self._state.plant_weight

    @plant_weight.setter
    def plant_weight(self, value):
        self._state.plant_weight = value

    @property
    def lint_yield(self):
        return self._state.lint_yield

    @lint_yield.setter
    def lint_yield(self, value):
        self._state.lint_yield = value

    @property
    def ginning_percent(self):
        return self._state.ginning_percent

    @ginning_percent.setter
    def ginning_percent(self, value):
        self._state.ginning_percent = value

    @property
    def number_of_squares(self):
        return self._state.number_of_squares

    @number_of_squares.setter
    def number_of_squares(self, value):
        self._state.number_of_squares = value

    @property
    def number_of_green_bolls(self):
        return self._state.number_of_green_bolls

    @number_of_green_bolls.setter
    def number_of_green_bolls(self, value):
        self._state.number_of_green_bolls = value

    @property
    def number_of_open_bolls(self):
        return self._state.number_of_open_bolls

    @number_of_open_bolls.setter
    def number_of_open_bolls(self, value):
        self._state.number_of_open_bolls = value

    @property
    def leaf_area_index(self):
        return self._state.leaf_area_index

    @leaf_area_index.setter
    def leaf_area_index(self, value):
        self._state.leaf_area_index = value

    @property
    def vegetative_branches(self):
        return [VegetativeBranch(self._state.vegetative_branches[i]) for i in range(self._state.number_of_vegetative_branches)]

    @property
    def hours(self):
        return self._state.hours

    @property
    def soil(self):
        return self._state.soil

    def __iter__(self):
        for attr in self.__slots__:
            value = getattr(self, attr)
            if value == 0 and attr.startswith("number_of_"):
                continue
            yield attr, value

cdef class Simulation:
    cdef cSimulation _sim

    def _doy2date(self, j):
        try:
            return datetime.strptime(f"{self.year} {j}", "%Y %j").date()
        except:
            return

    @property
    def year(self):
        return self._sim.year

    @year.setter
    def year(self, year):
        self._sim.year = year

    @property
    def start_date(self):
        return self._doy2date(self._sim.day_start)

    @start_date.setter
    def start_date(self, d):
        self._sim.day_start = _date2doy(d)

    @property
    def stop_date(self):
        return self._doy2date(self._sim.day_finish)

    @stop_date.setter
    def stop_date(self, d):
        self._sim.day_finish = _date2doy(d)

    @property
    def emerge_date(self):
        return self._doy2date(self._sim.day_emerge)

    @emerge_date.setter
    def emerge_date(self, d):
        self._sim.day_emerge = _date2doy(d)

    @property
    def plant_date(self):
        return self._doy2date(self._sim.day_plant)

    @plant_date.setter
    def plant_date(self, d):
        self._sim.day_plant = _date2doy(d)

    @property
    def topping_date(self):
        return self._doy2date(self._sim.day_topping)

    @topping_date.setter
    def topping_date(self, d):
        self._sim.day_topping = _date2doy(d)

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
        return VarPar

    @cultivar_parameters.setter
    def cultivar_parameters(self, parameters):
        for i, p in enumerate(parameters):
            VarPar[i + 1] = p

    @property
    def row_space(self):
        return self._sim.row_space

    @row_space.setter
    def row_space(self, value):
        self._sim.row_space = value or 0

    @property
    def skip_row_width(self):
        return SkipRowWidth

    @skip_row_width.setter
    def skip_row_width(self, value):
        global SkipRowWidth
        SkipRowWidth = value or 0

    @property
    def plants_per_meter(self):
        return PlantsPerM

    @plants_per_meter.setter
    def plants_per_meter(self, value):
        global PlantsPerM
        PlantsPerM = value or 0

    @property
    def states(self):
        return [State(self._sim.states[i]) for i in range(self._sim.day_finish - self._sim.day_start + 1)]

    @property
    def climate(self):
        return self._sim.climate

    @climate.setter
    def climate(self, climate):
        for i, daily_climate in enumerate(climate):
            self._sim.climate[i].Rad = daily_climate["radiation"]
            self._sim.climate[i].Tmax = daily_climate["max"]
            self._sim.climate[i].Tmin = daily_climate["min"]
            self._sim.climate[i].Wind = daily_climate["wind"]
            self._sim.climate[i].Rain = daily_climate["rain"]
            self._sim.climate[i].Tdew = daily_climate["dewpoint"]

    def run(self):
        app = new C2KApp()
        app.DailySimulation(self._sim)

    def read_input(self, lyrsol, **kwargs):
        """This is the main function for reading input."""
        InitializeGlobal()
        initialize_switch(self._sim)
        self._sim.states = <cState *> malloc(sizeof(cState) * (self._sim.day_finish - self._sim.day_start + 1))
        InitializeGrid(self._sim)
        read_agricultural_input(self._sim, kwargs.get("agricultural_inputs", []))
        InitializeSoilData(self._sim, lyrsol)
        InitializeSoilTemperature()
        InitializeRootData(self._sim)
