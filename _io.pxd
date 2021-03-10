from libcpp.string cimport string

from _structs cimport Simulation, ClimateStruct

cdef extern from "GettingInput_1.cpp":
    void InitializeGrid(Simulation &)
    double PlantsPerM
    double SkipRowWidth

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
