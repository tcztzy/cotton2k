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
    double amitri, ep;
} scratch;
extern scratch Scratch21[400];
typedef struct NitrogenFertilizer
{
        int day, mthfrt, ksdr, lsdr;
        double amtamm, amtnit, amtura;
} NitrogenFertilizer;
extern NitrogenFertilizer NFertilizer[150];
////    Integers    ////
extern int DayStartPredIrrig, DayStopPredIrrig, isw,
    LastDayWeatherData, LastIrrigation, LastTaprootLayer,
    LocationColumnDrip, LocationLayerDrip,
    MainStemNodes, MinDaysBetweenIrrig,
    nk, nl, noitr, NumAbscisedLeaves,
    NumIrrigations, NumNitApps,
    NumSheddingTags, NumWaterTableData, WaterTableLayer;
extern int CultivationDate[5], DayWaterTableInput[20], DefoliationDate[5], DefoliationMethod[5],
    LateralRootFlag[maxl], SoilHorizonNum[maxl];
////    Double    ////
extern double AverageLwp, AverageLwpMin, AverageSoilPsi,
    conmax, CumFertilizerN, CumNetPhotosynth, CumNitrogenUptake,
    CumWaterAdded, CumWaterDrained, DayTimeTemp, dclay,
    DepthLastRootLayer, dsand,
    ElCondSatSoilToday, GreenBollsLost,
    InitialTotalSoilWater, IrrigationDepth,
    LwpMax, LwpMin, MineralizedOrganicN,
    NetPhotosynthesis, NightTimeTemp,
    PercentDefoliation,
    PlantRowLocation, PotGroAllBolls, PotGroAllBurrs,
    PotGroAllLeaves, PotGroAllPetioles, PotGroAllRoots, PotGroAllSquares, PotGroStem,
    RatioImplicit, ReserveC, RootWeightLoss, SoilNitrogenLoss,
    SumNO3N90, TapRootLength, TotalSoilNh4N,
    TotalSoilNitrogen, TotalSoilNo3N, TotalSoilUreaN, TotalSoilWater;

extern double AbscissionLag[20], airdr[9],
    alpha[9], vanGenuchtenBeta[9], BulkDensity[9],
    ClayVolumeFraction[maxl], CultivationDepth[5],
    DefoliantAppRate[5], ElCondSatSoil[20],
    FieldCapacity[maxl], FoliageTemp[maxk],
    FreshOrganicNitrogen[maxl][maxk],
    HeatCapacitySoilSolid[maxl], HeatCondDrySoil[maxl], HumusNitrogen[maxl][maxk],
    HumusOrganicMatter[maxl][maxk],
    LevelsOfWaterTable[20], LwpMinX[3], LwpX[3],
    MarginalWaterContent[maxl], MaxWaterCapacity[maxl],
    NO3FlowFraction[maxl], PetioleWeightPreFru[9], PoreSpace[maxl],
    PotGroLeafAreaPreFru[9], PotGroLeafWeightPreFru[9], PotGroPetioleWeightPreFru[9],
    RootImpede[maxl][maxk], SandVolumeFraction[maxl], SaturatedHydCond[9], ShedByCarbonStress[20],
    ShedByNitrogenStress[20], ShedByWaterStress[20], SitePar[21], SoilPsi[maxl][maxk],
    SoilTemp[maxl][maxk], SoilTempDailyAvrg[maxl][maxk], StemWeight[365],
    thad[maxl], thetar[maxl], thetas[9], thts[maxl],
    VolNh4NContent[maxl][maxk], VolUreaNContent[maxl][maxk];

void InitializeGlobal();

extern "C"
{
    double dl(unsigned int);
    double wk(unsigned int, double);
}
