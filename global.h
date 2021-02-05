//   global.h
#pragma once

#include <exception>
#include <filesystem>
#include "stdafx.h"
#include "Irrigation.h"

namespace fs = std::filesystem;
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
            cumEvaporation, cumFertilizerN, cumNetPhotosynth, cumNitrogenUptake,
            cumPlantNLoss, cumTranspiration, cumWaterAdded, cumWaterDrained,
            dayTimeTemp, deadwt, deepSoilTemperature,
            ep, es, extraCarbon, fruitGrowthRatio,
            greenBollsLost, gbw, gintot, h2obal,
            leafAreaIndex, leafNConc, leafNitrogen, lightIntercept, lintYield, lwpMin,
            mineralizedOrganicN,
            netPhotosynthesis, nightTimeTemp, nitrogenStress, nStressFruiting,
            nStressRoots, nStressVeg, numGreenBolls, numSquares, numOpenBolls,
            petioleNConc, petioleNO3NConc, petioleNitrogen, PixInPlants, plantWeight,
            reserveC, rn, rootNConc, rootNitrogen, rootWeightLoss,
            seedNConc, seedNitrogen, soilNitrogenLoss, squareNConc, squareNitrogen,
            stemNConc, stemNitrogen, sumNO3N90, supplyNH4N, supplyNO3N,
            tapRootLength, totalLeafWeight, totalPetioleWeight, totalRequiredN,
            totalRootWeight, totalSoilNh4N, totalSoilNo3N, totalSoilUreaN,
            totalSoilWater, totalSquareWeight, totalStemWeight, waterStressStem;
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
    nk, nl, noitr, NumAbscisedLeaves, NumFruitSites,
    NumIrrigations, NumNitApps, NumPreFruNodes,
    NumSheddingTags, NumVegBranches, NumWaterTableData, SoilMapFreq, WaterTableLayer;
extern int CultivationDate[5], DayWaterTableInput[20], DefoliationDate[5], DefoliationMethod[5],
    LateralRootFlag[maxl], NumFruitBranches[3], NumNodes[3][30], OutIndex[24], pixday[10], pixmth[10],
    RootColNumLeft[maxl], RootColNumRight[maxl], SoilHorizonNum[maxl];
////    Boolean    ////
extern bool bPollinSwitch;
////    Double    ////
extern double ActualBollGrowth, ActualBurrGrowth,
    ActualSquareGrowth, ActualSoilEvaporation, ActualStemGrowth, ActualTranspiration, addwtbl,
    AppliedWater, AverageLwp, AverageLwpMin, AverageSoilPsi, AvrgDailyTemp,
    BloomWeightLoss, BurrNConc, BurrNitrogen, BurrWeightGreenBolls, BurrWeightOpenBolls,
    CarbonAllocatedForRootGrowth, conmax,
    CottonWeightGreenBolls, CottonWeightOpenBolls, CumEvaporation, CumFertilizerN,
    CumNetPhotosynth, CumNitrogenUptake, CumPlantNLoss, CumTranspiration,
    CumWaterAdded, CumWaterDrained,
    DayTimeTemp, dclay,
    DeepSoilTemperature, DensityFactor, DepthLastRootLayer, dsand,
    ElCondSatSoilToday, ExtraCarbon, FruitGrowthRatio,
    ginp, Gintot, GreenBollsLost,
    InitialTotalSoilWater, IrrigationDepth,
    LeafAreaIndex, LeafNConc, LeafNitrogen, LeafWeightAreaRatio,
    LightIntercept, LintYield, LwpMax, LwpMin,
    MaxIrrigation, MineralizedOrganicN,
    NetPhotosynthesis, NightTimeTemp, NitrogenStress, NStressFruiting,
    NStressRoots, NStressVeg, NumGreenBolls, NumOpenBolls, NumSquares,
    PercentDefoliation, PerPlantArea, PetioleNConc, PetioleNitrogen, PetioleNO3NConc,
    pixcon, pixda, pixdn, pixdz, PixInPlants, PlantPopulation, PlantRowLocation,
    PlantWeight, PlantWeightAtStart, PotGroAllBolls, PotGroAllBurrs,
    PotGroAllLeaves, PotGroAllPetioles, PotGroAllRoots, PotGroAllSquares, PotGroStem,
    RatioImplicit, ReferenceTransp, ReserveC, Rn, RootNConc, RootNitrogen,
    RootWeightLoss,
    SeedNConc, SeedNitrogen, SoilNitrogenAtStart, SoilNitrogenLoss, SquareNConc,
    SquareNitrogen, StemNConc, StemNitrogen, SumNO3N90, SupplyNH4N, SupplyNO3N,
    TapRootLength, TotalActualLeafGrowth, TotalActualPetioleGrowth, TotalLeafArea,
    TotalLeafWeight, TotalPetioleWeight, TotalRequiredN, TotalRootWeight, TotalSoilNh4N,
    TotalSoilNitrogen, TotalSoilNo3N, TotalSoilUreaN, TotalSoilWater, TotalSquareWeight,
    TotalStemWeight,
    WaterStressStem;

extern double AbscissionLag[20], AgeOfPreFruNode[9], airdr[9], AirTemp[24], albedo[24],
    alpha[9], vanGenuchtenBeta[9], BulkDensity[9],
    ClayVolumeFraction[maxl], CloudTypeCorr[24],
    CultivationDepth[5],
    DefoliantAppRate[5], DelayNewFruBranch[3], DelayNewNode[3][30],
    DewPointTemp[24], dl[maxl],
    ElCondSatSoil[20], es1hour[24], es2hour[24],
    FieldCapacity[maxl], FoliageTemp[maxk], FreshOrganicMatter[maxl][maxk],
    FreshOrganicNitrogen[maxl][maxk], gh2oc[10],
    HeatCapacitySoilSolid[maxl], HeatCondDrySoil[maxl], HumusNitrogen[maxl][maxk],
    HumusOrganicMatter[maxl][maxk],
    impede[10][10], LeafAreaMainStem[3][30], LeafAreaPreFru[9],
    LeafWeightMainStem[3][30], LeafWeightPreFru[9],
    LevelsOfWaterTable[20], LwpMinX[3], LwpX[3],
    MarginalWaterContent[maxl], MaxWaterCapacity[maxl], MulchTemp[maxk],
    NO3FlowFraction[maxl],
    PetioleWeightMainStem[3][30], PetioleWeightPreFru[9],
    pixppa[10], PoreSpace[maxl],
    PotGroLeafAreaMainStem[3][30],
    PotGroLeafAreaPreFru[9], PotGroLeafWeightMainStem[3][30],
    PotGroLeafWeightPreFru[9],
    PotGroPetioleWeightMainStem[3][30],
    PotGroPetioleWeightPreFru[9], PotGroRoots[maxl][maxk],
    ReferenceETP[24], RelativeHumidity[24], rlat1[maxl], rlat2[maxl],
    RootGroFactor[maxl][maxk], RootImpede[maxl][maxk],
    SandVolumeFraction[maxl], SaturatedHydCond[9], ShedByCarbonStress[20],
    ShedByNitrogenStress[20], ShedByWaterStress[20], SitePar[21], SoilPsi[maxl][maxk],
    SoilTemp[maxl][maxk], SoilTempDailyAvrg[maxl][maxk],
    StemWeight[365],
    thad[maxl], thetar[maxl], thetas[9], thts[maxl], tstbd[10][10],
    VarPar[61], VolNh4NContent[maxl][maxk], VolNo3NContent[maxl][maxk], VolUreaNContent[maxl][maxk],
    VolWaterContent[maxl][maxk],
    WindSpeed[24], wk[maxk];

void InitializeGlobal();
