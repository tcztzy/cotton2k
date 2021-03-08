from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdlib cimport malloc

from datetime import datetime, date

cdef extern void WriteInitialInputData(Simulation &, bool, double, double, double, const char *, int, const char *, const char *, const char *, const char *, const char *, const char *)

cdef extern from "global.h":
    void InitializeGlobal()
    int isw
    int Kday
    int SoilMapFreq
    int OutIndex[24]
    double PlantPopulation
    double TotalSoilNo3N
    double TotalSoilNh4N
    double TotalSoilUreaN
    double SoilNitrogenAtStart
    double PlantWeightAtStart
    double TotalRootWeight
    double TotalStemWeight
    double TotalLeafWeight
    double ReserveC

cdef extern from "State.h":
    ctypedef struct State:
        double plant_height
        double plant_weight
        double lint_yield

cdef extern from "Climate.h":
    ctypedef struct ClimateStruct:
        pass

cdef extern from "Simulation.h":
    ctypedef struct Simulation:
        const char *profile_name
        int year
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_emerge
        unsigned int day_plant
        unsigned int day_start_soil_maps
        unsigned int day_stop_soil_maps
        unsigned int day_start_co2
        unsigned int day_end_co2
        double co2_enrichment_factor
        unsigned int day_start_mulch
        unsigned int day_end_mulch
        unsigned int mulch_indicator
        double mulch_transmissivity_short_wave
        double mulch_transmissivity_long_wave
        double latitude
        double longitude
        double elevation
        double row_space
        State *states
        ClimateStruct climate[400]

cdef extern from "GettingInput_1.cpp":
    Simulation ReadProfileFile(const char *, vector[string] &)
    void ReadCalibrationData()
    void InitializeGrid(Simulation &)
    double PlantsPerM
    double SkipRowWidth
    int nSiteNum
    string SiteName
    int nVarNum
    string VarName

cdef extern from "GettingInput_2.cpp":
    void ReadSoilImpedance(Simulation &)
    void InitSoil(const string &)
    void InitializeSoilTemperature()
    void InitializeSoilData(Simulation &, const string &)
    void InitializeRootData(Simulation &)

cdef extern from "gettingInput_3.cpp":
    int OpenClimateFile(const string &, const string &, const int &, ClimateStruct[400])
    void ReadAgriculturalInput(Simulation &, const string &)

cdef extern from "Output.h":
    void DataOutput(Simulation &)
    void OpenOutputFiles(const char *, const char *, int, int);

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

cdef void read_site_config(
    Simulation &sim,
    double latitude,
    double longitude,
    double elevation,
    int site_id,
):
    sim.latitude = latitude
    sim.longitude = longitude
    sim.elevation = elevation
    global nSiteNum
    nSiteNum = site_id

cdef void read_plant_config(
    Simulation &sim,
    double row_space,
    double skip_row,
    double plants_per_meter,
    int var_id):
    global SkipRowWidth, PlantsPerM, nVarNum
    sim.row_space = row_space
    SkipRowWidth = skip_row
    PlantsPerM = plants_per_meter
    nVarNum = var_id

cdef void read_soil_map_config(Simulation &sim, int day_start, int day_stop, int frequency):
    global SoilMapFreq
    SoilMapFreq = frequency
    sim.day_start_soil_maps = day_start
    sim.day_stop_soil_maps = day_stop

def read_output_flags(output_flags):
    for n in range(24):
        OutIndex[n] = 0
    for i, flag in enumerate(output_flags):
        OutIndex[i+1] = flag

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
        return datetime.strptime(f"{self.year} {j}", "%Y %j").date()

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
        self._sim.day_start = d.timetuple().tm_yday

    @property
    def end_date(self):
        return self._doy2date(self._sim.day_finish)

    @end_date.setter
    def end_date(self, d):
        self._sim.day_finish = d.timetuple().tm_yday

    @property
    def emerge_date(self):
        return self._doy2date(self._sim.day_emerge)

    @emerge_date.setter
    def emerge_date(self, d):
        self._sim.day_emerge = d.timetuple().tm_yday

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
        self._sim = ReadProfileFile(profile_name, filenames)
        read_site_config(
            self._sim,
            kwargs.get("latitude", 0),
            kwargs.get("longitude", 0),
            kwargs.get("elevation", 0),
            kwargs["site_id"],
        )
        read_plant_config(
            self._sim,
            kwargs.get("row_space", 0),
            kwargs.get("skip_row", 0),
            kwargs.get("plants_per_meter", 0),
            kwargs["var_id"]
        )
        read_soil_map_config(
            self._sim,
            _date2doy(kwargs.get("soil_map_start_date", 0)),
            _date2doy(kwargs.get("soil_map_stop_date", 0)),
            kwargs.get("soil_map_frequency", 999),
        )
        read_output_flags(kwargs.get("output_flags", ()))
        initialize_switch(self._sim)
        _description = description.encode("utf-8")
        OpenOutputFiles(_description, <char *>self._sim.profile_name, self._sim.day_emerge, self._sim.year)
        self._sim.states = <State *> malloc(sizeof(State) * self._sim.day_finish - self._sim.day_start + 1)
        ReadCalibrationData()
        LastDayOfActualWeather = OpenClimateFile(filenames[0], filenames[1], self._sim.day_start, self._sim.climate)
        InitializeGrid(self._sim)
        ReadSoilImpedance(self._sim)
        WriteInitialInputData(self._sim, OutIndex[1], PlantsPerM, SkipRowWidth, PlantPopulation, filenames[0].c_str(), LastDayOfActualWeather, filenames[1].c_str(), filenames[4].c_str(), filenames[3].c_str(), filenames[2].c_str(), SiteName.c_str(), VarName.c_str())
        InitSoil(filenames[3])
        ReadAgriculturalInput(self._sim, filenames[4])
        InitializeSoilData(self._sim, filenames[2])
        InitializeSoilTemperature()
        InitializeRootData(self._sim)
        # initialize some variables at the start of simulation.
        SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN
        PlantWeightAtStart = TotalRootWeight + TotalStemWeight + TotalLeafWeight + ReserveC
