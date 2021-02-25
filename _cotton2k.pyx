from libcpp cimport bool
from libcpp.string cimport string
from libc.stdlib cimport malloc

from datetime import datetime

cdef extern void WriteInitialInputData(Simulation &, bool, double, double, double, const char *, int, const char *, const char *, const char *, const char *, const char *, const char *)

cdef extern from "global.h":
    void InitializeGlobal()
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
        int year
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_emerge
        State *states
        ClimateStruct climate[400]

cdef extern from "GettingInput_1.cpp":
    Simulation ReadProfileFile(const char *, string &, string &, string &, string &, string &)
    void ReadCalibrationData()
    void InitializeGrid(Simulation &)
    double PlantsPerM
    double SkipRowWidth
    string SiteName
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

cdef extern from "Cottonmodel.h":

    cdef cppclass C2KApp:
        C2KApp() except +
        void DailySimulation(Simulation &)

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


def read_input(str profile):
    """This is the main function for reading input."""
    cdef string ActWthFileName
    cdef string PrdWthFileName
    cdef string SoilHydFileName
    cdef string SoilInitFileName
    cdef string AgrInputFileName
    InitializeGlobal()
    _sim = ReadProfileFile(profile.encode("utf-8"), ActWthFileName, PrdWthFileName, SoilHydFileName, SoilInitFileName, AgrInputFileName)
    _sim.states = <State *> malloc(sizeof(State) * _sim.day_finish - _sim.day_start + 1)
    ReadCalibrationData()
    LastDayOfActualWeather = OpenClimateFile(ActWthFileName, PrdWthFileName, _sim.day_start, _sim.climate)
    InitializeGrid(_sim)
    ReadSoilImpedance(_sim)
    WriteInitialInputData(_sim, OutIndex[1], PlantsPerM, SkipRowWidth, PlantPopulation, ActWthFileName.c_str(), LastDayOfActualWeather, PrdWthFileName.c_str(), AgrInputFileName.c_str(), SoilInitFileName.c_str(), SoilHydFileName.c_str(), SiteName.c_str(), VarName.c_str())
    InitSoil(SoilInitFileName)
    ReadAgriculturalInput(_sim, AgrInputFileName)
    InitializeSoilData(_sim, SoilHydFileName)
    InitializeSoilTemperature()
    InitializeRootData(_sim)
    # initialize some variables at the start of simulation.
    SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN
    PlantWeightAtStart = TotalRootWeight + TotalStemWeight + TotalLeafWeight + ReserveC
    sim = _Simulation()
    sim._sim = _sim
    return sim
