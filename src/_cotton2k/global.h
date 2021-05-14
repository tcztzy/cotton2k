//   global.h
#pragma once

#include "Irrigation.h"
//
//  definition of global variables
//  ==============================
//    For dictionary of global variables see file "global.cpp"
////    Constants    ////
const int maxl = 40;
const int maxk = 20;
const double pi = 3.14159;
////    Structures    ////
typedef struct scratch
{
    double amitri, ep, es;
} scratch;
extern scratch Scratch21[400];
typedef struct NitrogenFertilizer
{
        int day, mthfrt, ksdr, lsdr;
        double amtamm, amtnit, amtura;
} NitrogenFertilizer;
extern NitrogenFertilizer NFertilizer[150];
////    Integers    ////
extern int DayStartPredIrrig, DayStopPredIrrig,
    inrim, IrrigMethod, isw,
    LastDayWeatherData, LastIrrigation, LastTaprootLayer,
    LocationColumnDrip, LocationLayerDrip,
    MainStemNodes, MinDaysBetweenIrrig,
    nk, nl, noitr, NumAbscisedLeaves,
    NumIrrigations, NumNitApps,
    NumSheddingTags, NumWaterTableData, WaterTableLayer;
extern int CultivationDate[5], DayWaterTableInput[20], DefoliationDate[5], DefoliationMethod[5],
    LateralRootFlag[maxl], SoilHorizonNum[maxl];
extern unsigned int ncurve;// number of input soil-moisture curves in the impedance table.
////    Double    ////
extern double ActualBollGrowth, ActualBurrGrowth,
    ActualSquareGrowth, ActualStemGrowth,
    AverageLwp, AverageLwpMin, AverageSoilPsi,
    BurrNConc, BurrNitrogen, BurrWeightGreenBolls, BurrWeightOpenBolls,
    CarbonAllocatedForRootGrowth, conmax,
    CottonWeightOpenBolls, CumFertilizerN,
    CumNetPhotosynth, CumNitrogenUptake,
    CumWaterAdded, CumWaterDrained, DayTimeTemp, dclay,
    DeepSoilTemperature, DepthLastRootLayer, dsand,
    ElCondSatSoilToday, Gintot, GreenBollsLost,
    InitialTotalSoilWater, IrrigationDepth,
    LightIntercept, LwpMax, LwpMin, MaxIrrigation, MineralizedOrganicN,
    NetPhotosynthesis, NightTimeTemp,
    PercentDefoliation, PetioleNitrogen, PetioleNO3NConc,
    PlantRowLocation, PotGroAllBolls, PotGroAllBurrs,
    PotGroAllLeaves, PotGroAllPetioles, PotGroAllRoots, PotGroAllSquares, PotGroStem,
    RatioImplicit, ReserveC, RootNitrogen, RootWeightLoss,
    SeedNitrogen, SoilNitrogenLoss,
    StemNConc, SumNO3N90, SupplyNH4N, SupplyNO3N,
    TapRootLength, TotalActualLeafGrowth, TotalActualPetioleGrowth,
    TotalPetioleWeight, TotalRootWeight, TotalSoilNh4N,
    TotalSoilNitrogen, TotalSoilNo3N, TotalSoilUreaN, TotalSoilWater, TotalSquareWeight;

extern double AbscissionLag[20], airdr[9],
    alpha[9], vanGenuchtenBeta[9], BulkDensity[9],
    ClayVolumeFraction[maxl], CultivationDepth[5],
    DefoliantAppRate[5], ElCondSatSoil[20],
    FieldCapacity[maxl], FoliageTemp[maxk],
    FreshOrganicNitrogen[maxl][maxk], gh2oc[10],
    HeatCapacitySoilSolid[maxl], HeatCondDrySoil[maxl], HumusNitrogen[maxl][maxk],
    HumusOrganicMatter[maxl][maxk],impede[10][10],
    LevelsOfWaterTable[20], LwpMinX[3], LwpX[3],
    MarginalWaterContent[maxl], MaxWaterCapacity[maxl],
    NO3FlowFraction[maxl], PetioleWeightPreFru[9], PoreSpace[maxl],
    PotGroLeafAreaPreFru[9], PotGroLeafWeightPreFru[9], PotGroPetioleWeightPreFru[9], rlat1[maxl], rlat2[maxl],
    RootImpede[maxl][maxk], SandVolumeFraction[maxl], SaturatedHydCond[9], ShedByCarbonStress[20],
    ShedByNitrogenStress[20], ShedByWaterStress[20], SitePar[21], SoilPsi[maxl][maxk],
    SoilTemp[maxl][maxk], SoilTempDailyAvrg[maxl][maxk], StemWeight[365],
    thad[maxl], thetar[maxl], thetas[9], thts[maxl], tstbd[10][10],
    VolNh4NContent[maxl][maxk], VolUreaNContent[maxl][maxk],
    VolWaterContent[maxl][maxk];

void InitializeGlobal();

extern "C"
{
    double dl(unsigned int);
    double wk(unsigned int, double);
}
