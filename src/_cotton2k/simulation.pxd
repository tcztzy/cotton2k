from libc.stdlib cimport malloc
from libc.math cimport exp
from libc.stdint cimport uint32_t

from _cotton2k._global cimport (
    BulkDensity,
    DefoliantAppRate,
    DefoliationDate,
    DefoliationMethod,
    InitializeGlobal,
    LastDayWeatherData,
    NFertilizer,
    NitrogenFertilizer,
    NumIrrigations,
    NumNitApps,
    PlantRowLocation,
    RatioImplicit,
    SaturatedHydCond,
    SitePar,
    airdr,
    alpha,
    conmax,
    isw,
    maxk,
    maxl,
    nk,
    nl,
    thetas,
    vanGenuchtenBeta,
)
from _cotton2k.climate cimport ClimateStruct
from _cotton2k.fruiting_site cimport Stage
from _cotton2k.irrigation cimport Irrigation
from _cotton2k.state cimport cState, Hour

cdef extern from "Simulation.hpp":
    ctypedef struct cSimulation "Simulation":
        int year
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_emerge
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
        unsigned int plant_row_column
        cState *states
        double cultivar_parameters[61]
        ClimateStruct climate[400]
        Irrigation irrigation[150]


cdef extern from "CottonPhenology.h":
    void CottonPhenology(cSimulation &, uint32_t)

cdef extern from "DailyClimate.h":
    void DayClim(cSimulation &, uint32_t u)

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

cdef extern from "PlantGrowth.h":
    void CheckDryMatterBal(cState &)
    void Defoliate(cSimulation &, uint32_t)
    void GetNetPhotosynthesis(cSimulation &, uint32_t, const double &)
    double PhysiologicalAge(Hour[24])
    void PlantGrowth(cSimulation &, const uint32_t &, const int &, const double &)
    void Stress(cSimulation &, unsigned int)

cdef extern from "PlantNitrogen.h":
    void PlantNitrogen(cSimulation &, uint32_t)

cdef extern from "SoilNitrogen.h":
    void SoilNitrogen(cSimulation &, unsigned int)

cdef extern from "SoilProcedures.h":
    void SoilProcedures(cSimulation &, uint32_t)
    void SoilSum(cState &, double)

cdef extern from "SoilTemperature.h":
    void ColumnShading(cState &, double[20], double, double, unsigned int)
    void SoilTemperature(cSimulation &, uint32_t, double[20])
