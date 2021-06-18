from libc.stdint cimport uint32_t

from .cxx cimport (
BulkDensity,
DeepSoilTemperature,
DefoliantAppRate,
DefoliationDate,
DefoliationMethod,
FieldCapacity,
PercentDefoliation,
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
SoilTemp,
airdr,
alpha,
conmax,
isw,
maxk,
maxl,
nk,
nl,
thad,
thetas,
vanGenuchtenBeta,
cSimulation,
)
from .fruiting_site cimport Stage
from .state cimport cState, Hour
from .soil cimport cSoilCell

cdef extern:
    double daytmp(cSimulation &, uint32_t, double, double, uint32_t, double, double)
    void AverageAirTemperatures(Hour[24], double &, double &, double &)
    double tdewhour(cSimulation &, uint32_t, uint32_t, double, double, double, double, double, double, double, double)
    double SimulateRunoff(cSimulation &, uint32_t, double, double, uint32_t)
    void EvapoTranspiration(cState &, double, double, double, double, double)

cdef extern from "CottonPhenology.h":
    void SimulateFruitingSite(cSimulation &, uint32_t, int, int, int, int &, const double &)
    void AddFruitingNode(cState &, int, int, double, double, double, double[61], double)
    void PreFruitingNode(cState &, double, double[61])

cdef extern from "FruitAbscission.h":
    void FruitingSitesAbscission(cSimulation &, uint32_t)

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

cdef extern from "LeafAbscission.h":
    void LeafAbscission(cSimulation &, uint32_t)

cdef extern from "PlantGrowth.h":
    void CheckDryMatterBal(cState &)
    double PhysiologicalAge(Hour[24])
    void DryMatterBalance(cState &, double &, double &, double &, double &, double)
    void ActualFruitGrowth(cState &)
    void ActualLeafGrowth(cState &)

cdef extern from "PlantNitrogen.h":
    void PlantNitrogen(cSimulation &, uint32_t)

cdef extern from "RootGrowth.h":
    double PotentialRootGrowth(cSoilCell[40][20], int, int, double)
    void ComputeActualRootGrowth(cState &, double, double, double, int, unsigned int, unsigned int)

cdef extern from "SoilNitrogen.h":
    void SoilNitrogen(cSimulation &, unsigned int)

cdef extern from "SoilProcedures.h":
    void SoilProcedures(cSimulation &, uint32_t)
    void SoilSum(cState &, double)

cdef extern from "SoilTemperature.h":
    void SoilTemperatureInit(cSimulation &)
    void PredictEmergence(cSimulation &, unsigned int, int)
    void SoilHeatFlux(cState &, double, int, int, int, int, double)
    void CanopyBalance(int, int, double, double, double, double, double, double, double, double &, const int &)
    void SoilSurfaceBalance(cState &, int, int, double, double, double, double, double, double &, double &, double &, double, double)
    double SensibleHeatTransfer(double, double, double, double)