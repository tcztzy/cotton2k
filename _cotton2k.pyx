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
