from libcpp.string cimport string

from _structs cimport Simulation, ClimateStruct

cdef extern from "GettingInput_2.cpp":
    void InitializeSoilTemperature()
    void InitializeSoilData(Simulation &, const string &)
    void InitializeRootData(Simulation &)
    double rnnh4[14]
    double rnno3[14]
    double oma[14]
    double h2oint[14]

cdef extern from "gettingInput_3.cpp":
    int OpenClimateFile(const string &, const string &, const int &, ClimateStruct[400])
    void ReadAgriculturalInput(Simulation &, const string &)
