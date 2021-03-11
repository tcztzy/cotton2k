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
        int kday, mainStemNodes;
        double amitri, averageSoilPsi, avrgDailyTemp, burrNConc, burrWeightOpenBolls,
            cottonWeightGreenBolls, cottonWeightOpenBolls,
            cumFertilizerN, cumNetPhotosynth, cumNitrogenUptake,
            cumPlantNLoss, cumWaterAdded, cumWaterDrained,
            dayTimeTemp, deadwt, deepSoilTemperature,
            ep, es, extraCarbon, fruitGrowthRatio,
            greenBollsLost, gbw, gintot, h2obal,
            leafAreaIndex, leafNitrogen, lightIntercept, lwpMin,
            mineralizedOrganicN,
            netPhotosynthesis, nightTimeTemp, nStressFruiting,
            nStressRoots, nStressVeg, petioleNO3NConc, petioleNitrogen,
            reserveC, rootNitrogen, rootWeightLoss,
            seedNitrogen, soilNitrogenLoss, squareNConc, squareNitrogen,
            stemNConc, stemNitrogen, sumNO3N90, supplyNH4N, supplyNO3N,
            tapRootLength, totalLeafWeight, totalPetioleWeight, totalRequiredN,
            totalRootWeight, totalSoilNh4N, totalSoilNo3N, totalSoilUreaN,
            totalSoilWater, totalSquareWeight, totalStemWeight;
        double soilPsi[maxl][maxk], soilTempDailyAvrg[maxl][maxk],
            volWaterContent[maxl][maxk], volNh4NContent[maxl][maxk], volNo3NContent[maxl][maxk];
} scratch;
extern scratch Scratch21[400];
typedef struct NitrogenFertilizer
{
        int day, mthfrt, ksdr, lsdr;
        double amtamm, amtnit, amtura;
} NitrogenFertilizer;
extern NitrogenFertilizer NFertilizer[150];
////    Integers    ////
extern int DayFirstDef,
    DayStartPredIrrig, DayStopPredIrrig,
    inrim, IrrigMethod, isw, Kday,
    LastDayWeatherData, LastIrrigation, LastTaprootLayer,
    LocationColumnDrip, LocationLayerDrip,
    MainStemNodes, MinDaysBetweenIrrig,
    nk, nl, noitr, NumAbscisedLeaves,
    NumIrrigations, NumNitApps, NumPreFruNodes,
    NumSheddingTags, NumWaterTableData, WaterTableLayer;
extern int CultivationDate[5], DayWaterTableInput[20], DefoliationDate[5], DefoliationMethod[5],
    LateralRootFlag[maxl], SoilHorizonNum[maxl];
////    Double    ////
extern double ActualBollGrowth, ActualBurrGrowth,
    ActualSquareGrowth, ActualStemGrowth, addwtbl,
    AverageLwp, AverageLwpMin, AverageSoilPsi, AvrgDailyTemp,
    BurrNConc, BurrNitrogen, BurrWeightGreenBolls, BurrWeightOpenBolls,
    CarbonAllocatedForRootGrowth, conmax,
    CottonWeightGreenBolls, CottonWeightOpenBolls, CumFertilizerN,
    CumNetPhotosynth, CumNitrogenUptake,
    CumWaterAdded, CumWaterDrained, DayTimeTemp, dclay,
    DeepSoilTemperature, DensityFactor, DepthLastRootLayer, dsand,
    ElCondSatSoilToday, FruitGrowthRatio, ginp, Gintot, GreenBollsLost,
    InitialTotalSoilWater, IrrigationDepth,
    LeafAreaIndex, LeafNitrogen, LeafWeightAreaRatio,
    LightIntercept, LwpMax, LwpMin, MaxIrrigation, MineralizedOrganicN,
    NetPhotosynthesis, NightTimeTemp, NStressFruiting, NStressRoots, NStressVeg,
    PercentDefoliation, PerPlantArea, PetioleNitrogen, PetioleNO3NConc,
    PlantPopulation, PlantRowLocation, PlantWeightAtStart, PotGroAllBolls, PotGroAllBurrs,
    PotGroAllLeaves, PotGroAllPetioles, PotGroAllRoots, PotGroAllSquares, PotGroStem,
    RatioImplicit, ReserveC, RootNitrogen, RootWeightLoss,
    SeedNitrogen, SoilNitrogenAtStart, SoilNitrogenLoss, SquareNConc,
    SquareNitrogen, StemNConc, StemNitrogen, SumNO3N90, SupplyNH4N, SupplyNO3N,
    TapRootLength, TotalActualLeafGrowth, TotalActualPetioleGrowth, TotalLeafArea,
    TotalLeafWeight, TotalPetioleWeight, TotalRootWeight, TotalSoilNh4N,
    TotalSoilNitrogen, TotalSoilNo3N, TotalSoilUreaN, TotalSoilWater, TotalSquareWeight, TotalStemWeight;

extern double AbscissionLag[20], AgeOfPreFruNode[9], airdr[9],
    alpha[9], vanGenuchtenBeta[9], BulkDensity[9],
    ClayVolumeFraction[maxl], CultivationDepth[5],
    DefoliantAppRate[5], DelayNewFruBranch[3], ElCondSatSoil[20],
    FieldCapacity[maxl], FoliageTemp[maxk],
    FreshOrganicNitrogen[maxl][maxk], gh2oc[10],
    HeatCapacitySoilSolid[maxl], HeatCondDrySoil[maxl], HumusNitrogen[maxl][maxk],
    HumusOrganicMatter[maxl][maxk],impede[10][10], LeafAreaPreFru[9],
     LeafWeightPreFru[9], LevelsOfWaterTable[20], LwpMinX[3], LwpX[3],
    MarginalWaterContent[maxl], MaxWaterCapacity[maxl],
    NO3FlowFraction[maxl], PetioleWeightPreFru[9], PoreSpace[maxl],
    PotGroLeafAreaPreFru[9], PotGroLeafWeightPreFru[9], PotGroPetioleWeightPreFru[9], rlat1[maxl], rlat2[maxl],
    RootImpede[maxl][maxk], SandVolumeFraction[maxl], SaturatedHydCond[9], ShedByCarbonStress[20],
    ShedByNitrogenStress[20], ShedByWaterStress[20], SitePar[21], SoilPsi[maxl][maxk],
    SoilTemp[maxl][maxk], SoilTempDailyAvrg[maxl][maxk], StemWeight[365],
    thad[maxl], thetar[maxl], thetas[9], thts[maxl], tstbd[10][10],
    VarPar[61], VolNh4NContent[maxl][maxk], VolUreaNContent[maxl][maxk],
    VolWaterContent[maxl][maxk];

void InitializeGlobal();

extern "C"
{
    double dl(unsigned int);
    double wk(unsigned int, double);
}
