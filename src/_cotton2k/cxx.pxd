from .climate cimport ClimateStruct
from .irrigation cimport Irrigation
from .state cimport cState

cdef extern from "Simulation.hpp":
    ctypedef struct cSimulation "Simulation":
        int year
        unsigned int day_emerge
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_plant
        unsigned int day_topping
        unsigned int day_defoliate
        double latitude
        double longitude
        double elevation
        double row_space
        double plant_population
        double per_plant_area
        double density_factor
        unsigned int first_bloom
        unsigned int first_square
        unsigned int plant_row_column
        double cultivar_parameters[61]
        ClimateStruct climate[400]
        Irrigation irrigation[150]
        cState *states

cdef extern from "global.h":
    ctypedef struct NitrogenFertilizer:
        int day
        int mthfrt
        int ksdr
        int lsdr
        double amtamm
        double amtnit
        double amtura
    void InitializeGlobal()
    unsigned int ncurve
    int inrim
    int isw
    int LastDayWeatherData
    const int maxl
    const int maxk
    int nl
    int nk
    double gh2oc[10]
    double tstbd[10][10]
    double impede[10][10]
    double SitePar[21]
    double TotalSoilNo3N
    double TotalSoilNh4N
    double TotalSoilUreaN
    double TotalRootWeight
    double ReserveC
    double PlantRowLocation
    double RatioImplicit
    double conmax
    double airdr[9]
    double thetas[9]
    double alpha[9]
    double vanGenuchtenBeta[9]
    double SaturatedHydCond[9]
    double BulkDensity[9]
    double DefoliantAppRate[5]
    int DefoliationDate[5]
    int DefoliationMethod[5]
    NitrogenFertilizer NFertilizer[150]
    int NumNitApps
    int NumIrrigations