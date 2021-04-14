from libcpp.string cimport string

from _cotton2k._structs cimport cSimulation, ClimateStruct

cdef extern from "GettingInput_2.cpp":
    void InitializeSoilTemperature()
    void InitializeSoilData(cSimulation &, unsigned int)
    void InitializeRootData(cSimulation &)
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
    double LayerDepth
