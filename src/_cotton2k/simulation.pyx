# distutils: language=c++
# cython: language_level=3
from libc.stdlib cimport malloc
from libc.math cimport exp
from datetime import datetime, date

from _cotton2k.cxx cimport cSimulation
from _cotton2k.climate cimport ClimateStruct
from _cotton2k.irrigation cimport Irrigation
from _cotton2k.rs cimport SlabLoc, tdewest, wk
from _cotton2k.utils import date2doy
from _cotton2k.state import State
from _cotton2k.state cimport cState


class SimulationEnd(RuntimeError):
    pass


cdef void SimulateThisDay(cSimulation &sim, uint32_t u):
    global isw
    cdef double rracol[20]  # the relative radiation received by a soil column, as affected by shading by plant canopy.
    if sim.day_emerge > 0 and sim.day_start + u >= sim.day_emerge:
        sim.states[u].kday = sim.day_start - sim.day_emerge + u + 1
    else:
        sim.states[u].kday = 0
    # The following functions are executed each day (also before emergence).
    ColumnShading(sim.states[u], rracol, sim.day_emerge, sim.row_space, sim.plant_row_column)  # computes light interception and soil shading.
    DayClim(sim, u)  # computes climate variables for today.
    SoilTemperature(sim, u, rracol)  # executes all modules of soil and canopy temperature.
    SoilProcedures(sim, u)  # executes all other soil processes.
    SoilNitrogen(sim, u)  # computes nitrogen transformations in the soil.
    SoilSum(sim.states[u], sim.row_space)  # computes totals of water and N in the soil.
    # The following is executed each day after plant emergence:
    if sim.states[u].daynum >= sim.day_emerge and isw > 0:
        # If this day is after emergence, assign to isw the value of 2.
        isw = 2
        sim.states[u].day_inc = PhysiologicalAge(sim.states[u].hours)  # physiological days increment for this day. computes physiological age
        Defoliate(sim, u)  # effects of defoliants applied.
        Stress(sim, u)  # computes water stress factors.
        GetNetPhotosynthesis(sim, u, sim.states[u].day_length)  # computes net photosynthesis.
        PlantGrowth(sim, u, 3, sim.states[u].day_length)  # executes all modules of plant growth.
        CottonPhenology(sim, u)  # executes all modules of plant phenology.
        PlantNitrogen(sim, u)  # computes plant nitrogen allocation.
        CheckDryMatterBal(sim.states[u]) # checks plant dry matter balance.
    # Check if the date to stop simulation has been reached, or if this is the last day with available weather data. Simulation will also stop when no leaves remain on the plant.
    if sim.states[u].daynum >= LastDayWeatherData or (sim.states[u].kday > 10 and sim.states[u].leaf_area_index < 0.0002):
        raise SimulationEnd


cdef CopyState(cSimulation &sim, uint32_t i):
    cdef cState state = sim.states[i]
    state.daynum = sim.day_start + i + 1
    sim.states[i + 1] = state


cdef DailySimulation(cSimulation &sim):
    cdef cState state0 = sim.states[0]
    state0.daynum = sim.day_start
    state0.lint_yield = 0
    state0.soil.number_of_layers_with_root = 7
    state0.plant_height = 4.0
    state0.plant_weight = 0
    state0.stem_weight = 0.2
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
        state0.vegetative_branches[k].number_of_fruiting_branches = 0
        for l in range(30):
            state0.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes = 0
            state0.vegetative_branches[k].fruiting_branches[l].delay_for_new_node = 0
            state0.vegetative_branches[k].fruiting_branches[l].main_stem_leaf = dict(
                leaf_area=0,
                leaf_weight=0,
                petiole_weight=0,
                potential_growth_for_leaf_area=0,
                potential_growth_for_leaf_weight=0,
                potential_growth_for_petiole_weight=0,
            )
            for m in range(5):
                state0.vegetative_branches[k].fruiting_branches[l].nodes[m] = dict(
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
    sim.states[0] = state0
    for i in range(sim.day_finish - sim.day_start + 1):
        SimulateThisDay(sim, i)
        if i < sim.day_finish - sim.day_start:
            CopyState(sim, i)


cdef void initialize_switch(cSimulation &sim):
    global isw
    cdef cState state0 = sim.states[0]
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
        sim.states[0].kday = 1

cdef double SkipRowWidth  # the smaller distance between skip rows, cm
cdef double PlantsPerM  # average number of plants pre meter of row.

cdef void InitializeGrid(cSimulation &sim):
    """
    This function initializes the soil grid variables. It is executed once at the beginning of the simulation. It is called from ReadInput().

    The following global or file-scope variables are set here:
    nk, nl.

    The following global variables are referenced here:
    maxk, maxl.
    """
    # PlantRowLocation is the distance from edge of slab, cm, of the plant row.
    global PlantRowLocation, nl, nk, SkipRowWidth, PlantsPerM
    PlantRowLocation = 0.5 * sim.row_space
    if (SkipRowWidth > 1):
        # If there is a skiprow arrangement, RowSpace and PlantRowLocation are redefined.
        sim.row_space = 0.5 * (sim.row_space + SkipRowWidth)  # actual width of the soil slab (cm)
        PlantRowLocation = 0.5 * SkipRowWidth
    # Compute sim.plant_population - number of plants per hectar, and per_plant_area - the average surface area per plant, in dm2, and the empirical plant density factor (density_factor). This factor will be used to express the effect of plant density on some plant growth rate functions.
    # NOTE: density_factor = 1 for 5 plants per sq m (or 50000 per ha).
    sim.plant_population = PlantsPerM / sim.row_space * 1000000
    sim.per_plant_area = 1000000 / sim.plant_population
    sim.density_factor = exp(sim.cultivar_parameters[1] * (5 - sim.plant_population / 10000))
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
            self.climate[i].Tdew = daily_climate.get("dewpoint", tdewest(daily_climate["max"], SitePar[5], SitePar[6]))

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


cdef read_agricultural_input(cSimulation &sim, inputs):
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
            DefoliationDate[idef] = date2doy(i["date"])
            DefoliantAppRate[idef] = -99.9
            if idef == 0:
                sim.day_defoliate = DefoliationDate[0]
            DefoliationMethod[idef] = i.get("method", 0)
            idef += 1


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
        self._sim.day_start = date2doy(d)

    @property
    def stop_date(self):
        return self._doy2date(self._sim.day_finish)

    @stop_date.setter
    def stop_date(self, d):
        self._sim.day_finish = date2doy(d)

    @property
    def emerge_date(self):
        return self._doy2date(self._sim.day_emerge)

    @emerge_date.setter
    def emerge_date(self, d):
        self._sim.day_emerge = date2doy(d)

    @property
    def plant_date(self):
        return self._doy2date(self._sim.day_plant)

    @plant_date.setter
    def plant_date(self, d):
        self._sim.day_plant = date2doy(d)

    @property
    def topping_date(self):
        return self._doy2date(self._sim.day_topping)

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
        return [State(state) for state in self._sim.states[:self._sim.day_finish - self._sim.day_start + 1]]

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
        DailySimulation(self._sim)

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
