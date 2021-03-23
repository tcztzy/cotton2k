from libcpp.string cimport string

from _structs cimport Simulation, ClimateStruct

cdef extern from "GettingInput_2.cpp":
    void InitializeSoilTemperature()
    void InitializeSoilData(Simulation &, const string &, unsigned int)
    void InitializeRootData(Simulation &)
    double rnnh4[14]
    double rnno3[14]
    double oma[14]
    double h2oint[14]
    double psisfc
    double psidra
    double ldepth[9]
    double condfc[9]
    double pclay[9]
    double psand[9]

cdef extern from "gettingInput_3.cpp":
    void ReadAgriculturalInput(Simulation &, const string &)
