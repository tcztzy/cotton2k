from libc.stdint cimport uint32_t

from .cxx cimport (
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
    cSimulation,
)
from .fruiting_site cimport Stage
from .state cimport cState, Hour

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
    void PlantGrowth(cState &, double, double, double, double[61], int, unsigned int, unsigned int, unsigned int, unsigned int)
    void Stress(cState &, double)

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
