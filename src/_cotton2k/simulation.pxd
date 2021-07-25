from libc.stdint cimport uint32_t

from .cxx cimport (
BulkDensity,
DefoliantAppRate,
DefoliationDate,
DefoliationMethod,
FieldCapacity,
PercentDefoliation,
InitializeGlobal,
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
from .fruiting_site cimport Stage, FruitingSite
from .state cimport cState, cHour
from .soil cimport cSoilCell, cSoil

cdef extern:
    double daytmp(cSimulation &, uint32_t, double, double, double, double)
    void AverageAirTemperatures(cHour[24], double &, double &, double &)
    double tdewhour(cSimulation &, uint32_t, double, double, double, double, double, double, double, double)
    double SimulateRunoff(cSimulation &, uint32_t, double, double, uint32_t)
    void EvapoTranspiration(cState &, double, double, double, double, double)

cdef extern from "FruitAbscission.h":
    void FruitingSitesAbscission(cSimulation &, uint32_t)

cdef extern from "GettingInput_2.cpp":
    void InitializeSoilTemperature()
    void InitializeSoilData(cSimulation &, unsigned int)
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
    void PreFruitLeafAbscission(cState &, double, unsigned int, unsigned int, unsigned int, double)
    void MainStemLeafAbscission(cState &, int, int, double, unsigned int, unsigned int)
    void DefoliationLeafAbscission(cState &, unsigned int)

cdef extern from "PlantGrowth.h":
    void DryMatterBalance(cState &, double &, double &, double &, double &, double)
    void ActualFruitGrowth(cState &)
    void ActualLeafGrowth(cState &)

cdef extern from "SoilNitrogen.h":
    void UreaHydrolysis(cSoilCell &, int, int)
    void MineralizeNitrogen(cSoilCell &, int, int, const int &, const int &, double)
    void Nitrification(cSoilCell &, int, int, double)
    void Denitrification(cSoilCell &, int, int, double)

cdef extern from "SoilProcedures.h":
    void ApplyFertilizer(cSimulation &, unsigned int)
    void RootsCapableOfUptake(cSoil &)
    double AveragePsi(const cState &, double)
    void WaterTable(cSimulation &, unsigned int)
    void WaterUptake(cSimulation &, unsigned int)
    void GravityFlow(cSoilCell[40][20], double, double)
    void CapillaryFlow(cSimulation &, unsigned int)
    void DripFlow(cSoilCell[40][20], double, double)

cdef extern from "SoilTemperature.h":
    void SoilTemperatureInit(cSimulation &)
    void PredictEmergence(cSimulation &, unsigned int, int)
    void SoilHeatFlux(cState &, double, int, int, int, int, double)
    void CanopyBalance(int, int, double, double, double, double, double, double, double, double &, const int &)
    double SensibleHeatTransfer(double, double, double, double)
    double ThermalCondSoil(double, double, int)
