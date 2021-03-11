from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdlib cimport malloc

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

cdef void read_mulch_config(
    Simulation &sim,
    int indicator,
    int day_start,
    int day_stop,
    double short_wave_transmissivity,
    double long_wave_transmissivity,
):
    sim.mulch_indicator = indicator
    sim.day_start_mulch = day_start
    sim.day_end_mulch = day_stop if day_stop >= 0 else sim.day_finish
    sim.mulch_transmissivity_short_wave = short_wave_transmissivity
    sim.mulch_transmissivity_long_wave = long_wave_transmissivity

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
    def output_flags(self):
        return OutIndex

    @output_flags.setter
    def output_flags(self, flags):
        for n in range(24):
            OutIndex[n] = 0
        for i, flag in enumerate(flags):
            OutIndex[i + 1] = flag

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
        # call DataOutput here because global variables varnish after run
        DataOutput(self._sim)

    def read_input(self, profile, description, **kwargs):
        """This is the main function for reading input."""
        cdef vector[string] filenames = [b'', b'', b'', b'', b'']
        InitializeGlobal()
        profile_name = profile.encode("utf-8")
        self._sim.profile_name = profile_name
        read_mulch_config(
            self._sim,
            kwargs.get("mulch_indicator", 0),
            _date2doy(kwargs.get("mulch_start_date", 0)),
            _date2doy(kwargs.get("mulch_stop_date", 0)),
            kwargs.get("mulch_short_wave_transmissivity", 0),
            kwargs.get("mulch_long_wave_transmissivity", 0),
        )
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
        OpenOutputFiles(_description, <char *>self._sim.profile_name, self._sim.day_emerge, self._sim.year)
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
