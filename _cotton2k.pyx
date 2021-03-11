from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdlib cimport malloc
from libc.math cimport exp

from datetime import datetime, date

from _global cimport *
from _structs cimport *
from _io cimport *

cdef extern from "Cottonmodel.h":

    cdef cppclass C2KApp:
        C2KApp() except +
        void DailySimulation(Simulation &)


def _date2doy(d):
    if isinstance(d, str):
        d = datetime.strptime(d, "%Y-%m-%d")
    if isinstance(d, date):
        return d.timetuple().tm_yday
    elif isinstance(d, int) and d > 0:
        return d
    else:
        return 0

cdef void read_filenames(
    vector[string] &filenames,
    string actual_weather_filename,
    string predicted_weather_filename,
    string soil_hydraulic_filename,
    string soil_init_filename,
    string agricultural_input_filename,
):
    filenames[0] = actual_weather_filename
    filenames[1] = predicted_weather_filename
    filenames[2] = soil_hydraulic_filename
    filenames[3] = soil_init_filename
    filenames[4] = agricultural_input_filename

cdef void read_soil_map_config(Simulation &sim, int day_start, int day_stop, int frequency):
    global SoilMapFreq
    SoilMapFreq = frequency
    sim.day_start_soil_maps = day_start
    sim.day_stop_soil_maps = day_stop

cdef void initialize_switch(Simulation &sim):
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

cdef void InitializeGrid(Simulation &sim):
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
    PlantPopulation = PlantsPerM / sim.row_space * 1000000;
    PerPlantArea = 1000000 / PlantPopulation;
    DensityFactor = exp(VarPar[1] * (5 - PlantPopulation / 10000));
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
        sumwk = sumwk + wk(k, sim.row_space);
        if sim.plant_row_column == 0 and sumwk > PlantRowLocation:
            if (sumwk - PlantRowLocation) > (0.5 * wk(k, sim.row_space)):
                sim.plant_row_column = k - 1
            else:
                sim.plant_row_column = k


cdef class _Simulation:
    cdef Simulation _sim

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
        return (self._sim.states[i] for i in range(self._sim.day_finish - self._sim.day_start + 1))

    def run(self):
        app = new C2KApp()
        app.DailySimulation(self._sim)

    def read_input(self, profile, description, **kwargs):
        """This is the main function for reading input."""
        cdef vector[string] filenames = [b'', b'', b'', b'', b'']
        InitializeGlobal()
        profile_name = profile.encode("utf-8")
        self._sim.profile_name = profile_name
        read_filenames(
            filenames,
            kwargs.get("actual_weather_filename", "").encode("UTF-8"),
            kwargs.get("predicted_weathre_filename", "").encode("UTF-8"),
            kwargs.get("soil_hydraulic_filename", "").encode("UTF-8"),
            kwargs.get("soil_init_filename", "").encode("UTF-8"),
            kwargs.get("agricultural_input_filename", "").encode("UTF-8"),
        )
        read_soil_map_config(
            self._sim,
            _date2doy(kwargs.get("soil_map_start_date", 0)),
            _date2doy(kwargs.get("soil_map_stop_date", 0)),
            kwargs.get("soil_map_frequency", 999),
        )
        initialize_switch(self._sim)
        _description = description.encode("utf-8")
        self._sim.states = <State *> malloc(sizeof(State) * self._sim.day_finish - self._sim.day_start + 1)
        LastDayOfActualWeather = OpenClimateFile(filenames[0], filenames[1], self._sim.day_start, self._sim.climate)
        InitializeGrid(self._sim)
        ReadSoilImpedance(self._sim)
        InitSoil(filenames[3])
        ReadAgriculturalInput(self._sim, filenames[4])
        InitializeSoilData(self._sim, filenames[2])
        InitializeSoilTemperature()
        InitializeRootData(self._sim)
        # initialize some variables at the start of simulation.
        SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN
        PlantWeightAtStart = TotalRootWeight + TotalStemWeight + TotalLeafWeight + ReserveC
