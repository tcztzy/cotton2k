# distutils: language=c++
# cython: language_level=3
from datetime import datetime, date, timedelta
from math import sin, cos, acos, sqrt, pi

from libc.math cimport exp
from libc.stdlib cimport malloc
from libcpp cimport bool as bool_t

from _cotton2k.climate import compute_day_length
from _cotton2k.leaf import temperature_on_leaf_growth_rate
from _cotton2k.photosynthesis import ambient_co2_factor
from _cotton2k.utils import date2doy
from .climate cimport ClimateStruct
from .cxx cimport (
    cSimulation,
    SandVolumeFraction,
    ClayVolumeFraction,
    DayTimeTemp,
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
    TotalPetioleWeight,
    PotGroAllSquares,
    PotGroAllBolls,
    PotGroAllBurrs,
    PotGroAllLeaves,
    PotGroAllPetioles,
    PotGroLeafAreaPreFru,
    PotGroLeafWeightPreFru,
    PotGroPetioleWeightPreFru,
)
from .irrigation cimport Irrigation
from .rs cimport SlabLoc, tdewest, wk, dayrad, dayrh, daywnd, PotentialStemGrowth, AddPlantHeight, TemperatureOnFruitGrowthRate
from .state cimport cState, cVegetativeBranch, cFruitingBranch


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
        return [self._branch.nodes[i] for i in
                range(self._branch.number_of_fruiting_nodes)]

    def __iter__(self):
        for attr in self.__slots__:
            yield attr, getattr(self, attr)

cdef class VegetativeBranch:
    cdef cVegetativeBranch _branch

    def __init__(self, _branch):
        self._branch = _branch

    @property
    def fruiting_branches(self):
        return [FruitingBranch(self._branch.fruiting_branches[i]) for i in
                range(self._branch.number_of_fruiting_branches)]

    def __iter__(self):
        yield "fruiting_branches", self.fruiting_branches

cdef class State:
    cdef cState *state

    def keys(self):
        return dict(self.state[0]).keys()

    def __getattr__(self, item):
        state = dict(self.state[0])
        return state[item]

    def __setattr__(self, key, value):
        state = dict(self.state[0])
        state[key] = value
        self.state[0] = state

    def __getitem__(self, item):
        return self.__getattr__(item)

    @property
    def average_temperature(self):
        return sum(hour["temperature"] for hour in self.hours) / 24

    @staticmethod
    cdef State from_ptr(cState *_ptr):
        """Factory function to create WrapperClass objects from
        given my_c_struct pointer.

        Setting ``owner`` flag to ``True`` causes
        the extension type to ``free`` the structure pointed to by ``_ptr``
        when the wrapper object is deallocated."""
        # Call to __new__ bypasses __init__ constructor
        cdef State state = State.__new__(State)
        state.state = _ptr
        return state

cdef class Simulation:
    cdef cSimulation _sim
    cdef public unsigned int version
    cdef double relative_radiation_received_by_a_soil_column[20]  # the relative radiation received by a soil column, as affected by shading by plant canopy.
    cdef double max_leaf_area_index
    cdef double ptsred  # The effect of moisture stress on the photosynthetic rate
    cdef double avtemp  # average temperature from day of emergence.
    cdef double sumstrs  # cumulative effect of water and N stresses on date of first square.
    cdef double DaysTo1stSqare   # number of days from emergence to 1st square
    cdef double defkgh  # amount of defoliant applied, kg per ha
    cdef double tdfkgh  # total cumulative amount of defoliant
    cdef bool_t idsw  # switch indicating if predicted defoliation date was defined.
    cdef public double skip_row_width  # the smaller distance between skip rows, cm
    cdef public double plants_per_meter  # average number of plants pre meter of row.

    def _doy2date(self, j):
        try:
            return datetime.strptime(f"{self.year} {j}", "%Y %j").date()
        except:
            return

    def __init__(self, version=0x0400):
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
    def states(self):
        return [self.state(i) for i in
                range(self._sim.day_finish - self._sim.day_start + 1)]

    def state(self, i):
        return State.from_ptr(&self._sim.states[i])

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
        cdef cState state0 = self.state(0)
        state0.daynum = self._sim.day_start
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
                state0.vegetative_branches[k].fruiting_branches[
                    l].number_of_fruiting_nodes = 0
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
        self._sim.states[0] = state0

    def _simulate(self):
        for i in range(self._sim.day_finish - self._sim.day_start + 1):
            self._simulate_this_day(i)
            if i < self._sim.day_finish - self._sim.day_start:
                CopyState(self._sim, i)

    def _calc_light_interception(self, u):
        cdef cState state = self._sim.states[u]
        if self.version < 0x0500:
            if state.leaf_area_index > self.max_leaf_area_index:
                self.max_leaf_area_index = state.leaf_area_index
            zint = 1.0756 * state.plant_height / self.row_space
            lfint = 0.80 * state.leaf_area_index if state.leaf_area_index <= 0.5 else 1 - exp(
                0.07 - 1.16 * state.leaf_area_index)
            if lfint > zint:
                light_interception = (zint + lfint) / 2
            elif state.leaf_area_index < self.max_leaf_area_index:
                light_interception = lfint
            else:
                light_interception = zint
            if light_interception > 1:
                light_interception = 1
            if light_interception < 0:
                light_interception = 0
            self._sim.states[u].light_interception = light_interception
        else:
            self._sim.states[u].light_interception = 1 - exp(-1.16 * state.leaf_area_index)

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
            self.relative_radiation_received_by_a_soil_column[k0] = 1 - shade
            if self.relative_radiation_received_by_a_soil_column[k0] < 0.05:
                self.relative_radiation_received_by_a_soil_column[k0] = 0.05

    def _daily_climate(self, u):
        cdef double declination  # daily declination angle, in radians
        cdef double sunr  # time of sunrise, hours.
        cdef double suns  # time of sunset, hours.
        cdef double tmpisr  # extraterrestrial radiation, \frac{W}{m^2}
        hour = timedelta(hours=1)
        result = compute_day_length((self.latitude, self.longitude), self._doy2date(self._sim.states[u].daynum))
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
            self._sim.states[u].hours[ihr].radiation = dayrad(ti, radsum, sinb, c11)
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
        tfrt = (self._sim.states[u].day_length * TemperatureOnFruitGrowthRate(DayTimeTemp) + (24 - self._sim.states[u].day_length) * TemperatureOnFruitGrowthRate(NightTimeTemp)) / 24;
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
                        ratebol = 4 * tfrt * rbmax * pex / (1 + pex) ** 2;
                        # Potential growth rate of the burrs is assumed to be constant (vpotfrt[4] g dry weight per day) until the boll reaches its final volume. This occurs at the age of 22 physiological days in 'Acala-SJ2'. Both ratebol and ratebur are modified by temperature (tfrt) and ratebur is also affected by water stress (wfdb).
                        # Compute wfdb for the effect of water stress on burr growth rate. wfdb is the effect of water stress on rate of burr growth.
                        wfdb = vpotfrt[0] + vpotfrt[1] * self._sim.states[u].water_stress;
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
        global PotGroStem, PotGroAllRoots, TotalPetioleWeight
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
        sumpdr = PotentialRootGrowth(self._sim.states[u].soil.cells, 3, self._sim.states[u].soil.number_of_layers_with_root,
                                 self._sim.per_plant_area)
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
        cdef double cdstem, cdleaf, cdpet, cdroot;
        DryMatterBalance(self._sim.states[u], cdstem, cdleaf, cdpet, cdroot, self._sim.per_plant_area)
        # If it is after first square, call ActualFruitGrowth() to compute actual
        # growth rate of squares and bolls.
        if self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].stage != Stage.NotYetFormed:
            ActualFruitGrowth(self._sim.states[u])
        # Initialize state.leaf_weight.It is assumed that cotyledons fall off at time of first square.Also initialize state.leaf_area and TotalPetioleWeight.
        if self._sim.first_square > 0:
            self._sim.states[u].leaf_weight = 0
            self._sim.states[u].leaf_area = 0
        else:
            cotylwt = 0.20  # weight of cotyledons dry matter.
            self._sim.states[u].leaf_weight = cotylwt
            self._sim.states[u].leaf_area = 0.6 * cotylwt
        TotalPetioleWeight = 0
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
        if self._sim.day_topping <= 0 or self._sim.states[u].daynum < self._sim.day_topping:
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
        ComputeActualRootGrowth(self._sim.states[u], sumpdr, self.row_space, self._sim.per_plant_area, 3, self._sim.day_emerge,
                                self._sim.plant_row_column)

    def _days_to_first_square(self, u):
        state = self.state(u)
        if state.daynum <= self._sim.day_emerge:
            self.avtemp = state.average_temperature
            self.sumstrs = 0
        self.avtemp = ((state.daynum - self._sim.day_emerge) * self.avtemp + state.average_temperature) / (state.daynum - self._sim.day_emerge + 1)
        if self.avtemp > 34:
            self.avtemp = 34
        self.sumstrs += 0.08 * (1 - state.water_stress) * 0.3 * (1 - state.nitrogen_stress_vegetative)
        return (132.2 + self.avtemp * (-7 + self.avtemp * 0.125)) * self.cultivar_parameters[30] - self.sumstrs

    def _create_first_square(self, u, stemNRatio):
        """
        This function initiates the first square. It is called from function CottonPhenology().
        """
        # FruitFraction and FruitingCode are assigned 1 for the first fruiting site.
        self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].stage = Stage.Square
        self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].fraction = 1
        # Initialize a new leaf at this position. define its initial weight and area. VarPar[34] is the initial area of a new leaf. The mass and nitrogen of the new leaf are substacted from the stem.
        self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].leaf.area = self.cultivar_parameters[34]
        self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].leaf.weight = self.cultivar_parameters[34] * self._sim.states[u].leaf_weight_area_ratio
        self._sim.states[u].stem_weight -= self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].leaf.weight
        self._sim.states[u].leaf_weight += self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].leaf.weight
        self._sim.states[u].leaf_nitrogen += self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].leaf.weight * stemNRatio
        self._sim.states[u].stem_nitrogen -= self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].leaf.weight * stemNRatio
        # Define the initial values of NumFruitBranches, NumNodes, state.fruit_growth_ratio, and AvrgNodeTemper.
        self._sim.states[u].vegetative_branches[0].number_of_fruiting_branches = 1
        self._sim.states[u].vegetative_branches[0].fruiting_branches[0].number_of_fruiting_nodes = 1
        self._sim.states[u].fruit_growth_ratio = 1
        self._sim.states[u].vegetative_branches[0].fruiting_branches[0].nodes[0].average_temperature = self._sim.states[u].average_temperature
        # It is assumed that the cotyledons are dropped at time of first square. compute changes in AbscisedLeafWeight, state.leaf_weight, state.leaf_nitrogen and CumPlantNLoss caused by the abscission of the cotyledons.
        cdef double cotylwt = 0.20  # cotylwt is the leaf weight of the cotyledons.
        self._sim.states[u].abscised_leaf_weight += cotylwt
        self._sim.states[u].leaf_weight -= cotylwt
        self._sim.states[u].cumulative_nitrogen_loss += cotylwt * self._sim.states[u].leaf_nitrogen / self._sim.states[u].leaf_weight
        self._sim.states[u].leaf_nitrogen -= cotylwt * self._sim.states[u].leaf_nitrogen / self._sim.states[u].leaf_weight

    def _phenology(self, u):
        """
        This is is the main function for simulating events of phenology and abscission in the cotton plant. It is called each day from DailySimulation().

        CottonPhenology() calls PreFruitingNode(), DaysToFirstSquare(), _create_first_square(), AddVegetativeBranch(), AddFruitingBranch(), AddFruitingNode(), SimulateFruitingSite(), LeafAbscission(), FruitingSitesAbscission().
        """
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
        delayVegByCStress = self._sim.cultivar_parameters[27] + self._sim.states[u].carbon_stress * (vpheno[3] + vpheno[4] * self._sim.states[u].carbon_stress)
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
        # The following section is executed if the first square has not yet been formed. Function DaysToFirstSquare() is called to compute the  number of days to 1st square, and function PreFruitingNode() is called to simulate the formation of prefruiting nodes.
        if self._sim.first_square <= 0:
            self.DaysTo1stSqare = self._days_to_first_square(u)
            PreFruitingNode(self._sim.states[u], stemNRatio, self._sim.cultivar_parameters)
            # When first square is formed, FirstSquare is assigned the day of year.
            # Function _create_first_square() is called for formation of first square.
            if self._sim.states[u].kday >= <int>self.DaysTo1stSqare:
                self._sim.first_square = self._sim.states[u].daynum
                self._create_first_square(u, stemNRatio)
            # if a first square has not been formed, call LeafAbscission() and exit.
            else:
                LeafAbscission(self._sim, u)
                return
        # The following is executed after the appearance of the first square.
        # If there are only one or two vegetative branches, and if plant population allows it, call AddVegetativeBranch() to decide if a new vegetative branch is to be added. Note that dense plant populations (large per_plant_area) prevent new vegetative branch formation.
        if self._sim.states[u].number_of_vegetative_branches == 1 and self._sim.per_plant_area >= vpheno[5]:
            AddVegetativeBranch(self._sim.states[u], delayVegByCStress, stemNRatio, self.DaysTo1stSqare, self._sim.cultivar_parameters, PhenDelayByNStress)
        if self._sim.states[u].number_of_vegetative_branches == 2 and self._sim.per_plant_area >= vpheno[6]:
            AddVegetativeBranch(self._sim.states[u], delayVegByCStress, stemNRatio, self.DaysTo1stSqare, self._sim.cultivar_parameters, PhenDelayByNStress)
        # The maximum number of nodes per fruiting branch (nidmax) is affected by plant density. It is computed as a function of density_factor.
        cdef int nidmax  # maximum number of nodes per fruiting branch.
        nidmax = <int>(vpheno[7] * self._sim.density_factor + 0.5)
        if nidmax > 5:
            nidmax = 5
        # Start loop over all existing vegetative branches.
        # Call AddFruitingBranch() to decide if a new node (and a new fruiting branch) is to be added on this stem.
        for k in range(self._sim.states[u].number_of_vegetative_branches):
            if self._sim.states[u].vegetative_branches[k].number_of_fruiting_branches < 30:
                AddFruitingBranch(self._sim.states[u], k, delayVegByCStress, stemNRatio, self._sim.density_factor, self._sim.cultivar_parameters, PhenDelayByNStress)
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
        LeafAbscission(self._sim, u)

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
        if 0 < self._sim.day_emerge <= self._sim.day_start + u:
            self._sim.states[u].kday = (self._sim.day_start + u) - self._sim.day_emerge + 1
            self._calc_light_interception(u)
            self._column_shading(u)
        else:
            self._sim.states[u].kday = 0
            self._sim.states[u].light_interception = 0
            self.relative_radiation_received_by_a_soil_column = [1] * 20
        # The following functions are executed each day (also before emergence).
        self._daily_climate(u)  # computes climate variables for today.
        SoilTemperature(self._sim, u,
                        self.relative_radiation_received_by_a_soil_column)  # executes all modules of soil and canopy temperature.
        SoilProcedures(self._sim, u)  # executes all other soil processes.
        SoilNitrogen(self._sim, u)  # computes nitrogen transformations in the soil.
        SoilSum(self._sim.states[u], self._sim.row_space)  # computes totals of water and N in the soil.
        # The following is executed each day after plant emergence:
        if self._sim.states[u].daynum >= self._sim.day_emerge and isw > 0:
            # If this day is after emergence, assign to isw the value of 2.
            isw = 2
            self._sim.states[u].day_inc = PhysiologicalAge(self._sim.states[
                                                               u].hours)  # physiological days increment for this day. computes physiological age
            self._defoliate(u)  # effects of defoliants applied.
            self._stress(u)  # computes water stress factors.
            self._get_net_photosynthesis(u)  # computes net photosynthesis.
            self._growth(u)  # executes all modules of plant growth.
            self._phenology(u)  # executes all modules of plant phenology.
            PlantNitrogen(self._sim, u)  # computes plant nitrogen allocation.
            CheckDryMatterBal(self._sim.states[u])  # checks plant dry matter balance.
        # Check if the date to stop simulation has been reached, or if this is the last day with available weather data. Simulation will also stop when no leaves remain on the plant.
        if self._sim.states[u].daynum >= LastDayWeatherData or (
                self._sim.states[u].kday > 10 and self._sim.states[u].leaf_area_index < 0.0002):
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
        self._sim.density_factor = exp(
            self._sim.cultivar_parameters[1] * (5 - self._sim.plant_population / 10000))
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
