//   global.h
#pragma once
#include <string>
//
//  definition of global variables
//  ==============================
//    For dictionary of global variables see file "global.cpp"
////    Constants    ////
const int maxl = 40;
const int maxk = 20;
const double pi = 3.14159;
extern unsigned int version;
////    Structures    ////
extern struct scratch {
    int daynum, kday, firstBloom, firstSquare, lastTaprootLayer, mainStemNodes,
        numAbscisedLeaves, numFruitSites, numLayersWithRoots, numPreFruNodes,
        numSheddingTags, numVegBranches;
    int fruitingCode[3][30][5], lateralRootFlag[maxl], numFruitBranches[3],
        numNodes[3][30], rootColNumLeft[maxl], rootColNumRight[maxl];
    std::string date;
    double abscisedFruitSites, abscisedLeafWeight, amitri, averageLwp,
        averageLwpMin, averageSoilPsi, avrgDailyTemp, bloomWeightLoss,
        burrNConc, burrNitrogen, burrWeightGreenBolls, burrWeightOpenBolls,
        carbonStress, cottonWeightGreenBolls, cottonWeightOpenBolls,
        cumEvaporation, cumFertilizerN, cumNetPhotosynth, cumNitrogenUptake,
        cumPlantNLoss, cumTranspiration, cumWaterAdded, cumWaterDrained,
        dayTimeTemp, deadwt, deepSoilTemperature, ep, es, extraCarbon,
        fruitGrowthRatio, greenBollsLost, gbw, gintot, h2obal, leafAreaIndex,
        leafNConc, leafNitrogen, lightIntercept, lintYield, lwpMin,
        mineralizedOrganicN, netPhotosynthesis, nightTimeTemp, nitrogenStress,
        nStressFruiting, nStressRoots, nStressVeg, numGreenBolls, numSquares,
        numOpenBolls, petioleNConc, petioleNO3NConc, petioleNitrogen,
        PixInPlants, plantHeight, plantWeight, rad, rain, reserveC, rn,
        rootNConc, rootNitrogen, rootWeightLoss, runoff, seedNConc,
        seedNitrogen, soilNitrogenLoss, squareNConc, squareNitrogen, stemNConc,
        stemNitrogen, sumNO3N90, supplyNH4N, supplyNO3N, tapRootLength, tmax,
        tmin, totalPetioleWeight, totalRequiredN,
        totalRootWeight, totalSoilNh4N, totalSoilNo3N, totalSoilUreaN,
        totalSoilWater, totalSquareWeight, totalStemWeight, waterStress,
        waterStressStem, wind;
    double abscissionLag[20], ageOfBoll[3][30][5], ageOfPreFruNode[9],
        ageOfSite[3][30][5], avrgNodeTemper[3][30][5], bollWeight[3][30][5],
        burrWeight[3][30][5], delayNewFruBranch[3], delayNewNode[3][30],
        foliageTemp[maxk], freshOrganicMatter[maxl][maxk],
        freshOrganicNitrogen[maxl][maxk], fruitFraction[3][30][5],
        humusOrganicMatter[maxl][maxk], leafAge[3][30][5],
        leafAreaMainStem[3][30], leafAreaNodes[3][30][5], leafAreaPreFru[9],
        leafWeightMainStem[3][30], leafWeightNodes[3][30][5],
        leafWeightPreFru[9], lwpMinX[3], lwpX[3], mulchTemp[maxk],
        nhum[maxl][maxk], petioleWeightMainStem[3][30],
        petioleWeightNodes[3][30][5], petioleWeightPreFru[9],
        potGroLeafAreaMainStem[3][30], potGroLeafWeightMainStem[3][30],
        potGroPetioleWeightMainStem[3][30], rootAge[maxl][maxk],
        rootsv[maxl][maxk], rootWeight[maxl][maxk][3],
        rootWtCapblUptake[maxl][maxk], shedByCarbonStress[20],
        shedByNitrogenStress[20], shedByWaterStress[20], soilPsi[maxl][maxk],
        soilTemp[maxl][maxk], soilTempDailyAvrg[maxl][maxk],
        squareWeight[3][30][5], volWaterContent[maxl][maxk],
        volNh4NContent[maxl][maxk], volNo3NContent[maxl][maxk],
        volUreaNContent[maxl][maxk];
} Scratch21[400];
//
extern struct Climstruct {
    int nDay;
    double Rad, Tmax, Tmin, Rain, Wind, Tdew;
} Clim[400];
extern struct NitrogenFertilizer {
    int day, mthfrt, ksdr, lsdr;
    double amtamm, amtnit, amtura;
} NFertilizer[150];
extern struct Irrigation {
    int day, method, LocationColumnDrip, LocationLayerDrip;
    double amount;
} Irrig[150];
////    Strings    ////
extern std::string AgrInputFileName, Date,
    PlantmapFileName, ProfileName, SoilHydFileName,
    SoilInitFileName;
////    Integers    ////
extern int DayEmerge, DayEndCO2, DayEndMulch, DayFinish, DayFirstDef,
    DayOfSimulation, Daynum, DayPlant, DayStart, DayStartCO2, DayStartMulch,
    DayStartPlantMaps, DayStartPredIrrig, DayStartSoilMaps, DayStopPlantMaps,
    DayStopPredIrrig, DayStopSoilMaps, FirstBloom, FirstSquare, inrim,
    IrrigMethod, isw, iyear, Kday, KdayAdjust, LastDayWeatherData,
    LastIrrigation, LastTaprootLayer, LocationColumnDrip, LocationLayerDrip,
    MainStemNodes, MinDaysBetweenIrrig, MulchIndicator, ncurve, nk, nl, noitr,
    NumAbscisedLeaves, NumAdjustDays, NumFruitSites, NumIrrigations,
    NumLayersWithRoots, NumNitApps, NumPreFruNodes, NumRootAgeGroups,
    NumSheddingTags, NumVegBranches, NumWaterTableData, PlantMapFreq,
    PlantRowColumn, SoilMapFreq, WaterTableLayer;
extern int CultivationDate[5], DayWaterTableInput[20], DefoliationDate[5],
    DefoliationMethod[5], FruitingCode[3][30][5], LateralRootFlag[maxl],
    MapDataDate[30], NumFruitBranches[3], NumNodes[3][30], NodeLayer[3][30],
    NodeLayerPreFru[9], OutIndex[24],
    pixday[10], pixmth[10], RootColNumLeft[maxl], RootColNumRight[maxl],
    SoilHorizonNum[maxl];
////    Boolean    ////
extern bool bEnd, bPollinSwitch, nadj[5];
////    Double    ////
extern double AbscisedFruitSites, AbscisedLeafWeight, ActualBollGrowth,
    ActualBurrGrowth, ActualSquareGrowth, ActualSoilEvaporation,
    ActualStemGrowth, ActualTranspiration, addwtbl, AdjAddHeightRate,
    AdjAddMSNodesRate, AdjAddSitesRate, AdjGreenBollAbsc, AdjSquareAbsc,
    AppliedWater, AverageLwp, AverageLwpMin, AverageSoilPsi, AvrgDailyTemp,
    BloomWeightLoss, BurrNConc, BurrNitrogen, BurrWeightGreenBolls,
    BurrWeightOpenBolls, CarbonAllocatedForRootGrowth, CarbonStress,
    CO2EnrichmentFactor, conmax, CottonWeightGreenBolls, CottonWeightOpenBolls,
    CumEvaporation, CumFertilizerN, CumNetPhotosynth, CumNitrogenUptake,
    CumPlantNLoss, CumTranspiration, CumWaterAdded, CumWaterDrained,
    DailyRootLoss, DayInc, DayLength, DayTimeTemp, dclay, DeepSoilTemperature,
    DensityFactor, DepthLastRootLayer, dsand, ElCondSatSoilToday, Elevation,
    ExtraCarbon, FruitGrowthRatio, ginp, Gintot, GreenBollsLost,
    InitialTotalSoilWater, IrrigationDepth, Latitude, LeafAreaIndex, LeafNConc,
    LeafNitrogen, LeafWeightAreaRatio, LightIntercept, LintYield, Longitude,
    LwpMax, LwpMin, MaxIrrigation, MineralizedOrganicN, MulchTranLW,
    MulchTranSW, NetPhotosynthesis, NightTimeTemp, NitrogenStress,
    NStressFruiting, NStressRoots, NStressVeg, NumGreenBolls, NumOpenBolls,
    NumSquares, PercentDefoliation, PerPlantArea, PetioleNConc, PetioleNitrogen,
    PetioleNO3NConc, pixcon, pixda, pixdn, pixdz, PixInPlants, PlantHeight,
    PlantPopulation, PlantRowLocation, PlantWeight, PlantWeightAtStart,
    PotGroAllBolls, PotGroAllBurrs, PotGroAllLeaves, PotGroAllPetioles,
    PotGroAllRoots, PotGroAllSquares, PotGroStem, RatioImplicit,
    ReferenceTransp, ReserveC, Rn, RootNConc, RootNitrogen, RootWeightLoss,
    RowSpace, SeedNConc, SeedNitrogen, SoilNitrogenAtStart, SoilNitrogenLoss,
    SquareNConc, SquareNitrogen, StemNConc, StemNitrogen, SumNO3N90, SupplyNH4N,
    SupplyNO3N, TapRootLength, TotalActualLeafGrowth, TotalActualPetioleGrowth,
    TotalPetioleWeight, TotalRequiredN,
    TotalRootWeight, TotalSoilNh4N, TotalSoilNitrogen, TotalSoilNo3N,
    TotalSoilUreaN, TotalSoilWater, TotalSquareWeight, TotalStemWeight,
    WaterStress, WaterStressStem;

extern double AbscissionLag[20], ActualRootGrowth[maxl][maxk],
    AgeOfBoll[3][30][5], AgeOfPreFruNode[9], AgeOfSite[3][30][5], airdr[9],
    AirTemp[24], albedo[24], alpha[9], AvrgNodeTemper[3][30][5], beta[9],
    BollWeight[3][30][5], BulkDensity[9], BurrWeight[3][30][5], cgind[3],
    ClayVolumeFraction[maxl], CloudCoverRatio[24], CloudTypeCorr[24],
    CultivationDepth[5], DefoliantAppRate[5], DelayNewFruBranch[3],
    DelayNewNode[3][30], DewPointTemp[24], dl[maxl], ElCondSatSoil[20],
    es1hour[24], es2hour[24], FieldCapacity[maxl], FoliageTemp[maxk],
    FreshOrganicMatter[maxl][maxk], FreshOrganicNitrogen[maxl][maxk],
    FruitFraction[3][30][5], gh2oc[10], HeatCapacitySoilSolid[maxl],
    HeatCondDrySoil[maxl], HumusNitrogen[maxl][maxk],
    HumusOrganicMatter[maxl][maxk], impede[10][10], LeafAge[3][30][5],
    LeafAreaMainStem[3][30], LeafArea[20], LeafAreaIndexes[20], LeafAreaNodes[3][30][5], LeafAreaPreFru[9], LeafNitrogenLayer[20],
    LeafWeightLayer[20], LeafWeightMainStem[3][30], LeafWeightNodes[3][30][5], LeafWeightPreFru[9], LightInterceptLayer[20],
    LevelsOfWaterTable[20], LwpMinX[3], LwpX[3], MarginalWaterContent[maxl],
    MapDataAllSiteNum[30], MapDataGreenBollNum[30], MapDataMainStemNodes[30],
    MapDataPlantHeight[30], MapDataSquareNum[30], MaxWaterCapacity[maxl],
    MulchTemp[maxk], NO3FlowFraction[maxl], PetioleWeightMainStem[3][30],
    PetioleWeightNodes[3][30][5], PetioleWeightPreFru[9], pixppa[10],
    PoreSpace[maxl], PotGroBolls[3][30][5], PotGroBurrs[3][30][5],
    PotGroLeafAreaMainStem[3][30], PotGroLeafAreaNodes[3][30][5],
    PotGroLeafAreaPreFru[9], PotGroLeafWeightMainStem[3][30],
    PotGroLeafWeightNodes[3][30][5], PotGroLeafWeightPreFru[9],
    PotGroPetioleWeightMainStem[3][30], PotGroPetioleWeightNodes[3][30][5],
    PotGroPetioleWeightPreFru[9], PotGroRoots[maxl][maxk],
    PotGroSquares[3][30][5], Radiation[24], ReferenceETP[24],
    RelativeHumidity[24], rlat1[maxl], rlat2[maxl], RootAge[maxl][maxk],
    RootGroFactor[maxl][maxk], RootImpede[maxl][maxk],
    RootWeight[maxl][maxk][3], RootWtCapblUptake[maxl][maxk], rracol[maxk],
    SandVolumeFraction[maxl], SaturatedHydCond[9], ShedByCarbonStress[20],
    ShedByNitrogenStress[20], ShedByWaterStress[20], SitePar[21],
    SoilPsi[maxl][maxk], SoilTemp[maxl][maxk], SoilTempDailyAvrg[maxl][maxk],
    SquareWeight[3][30][5], StemWeight[365], thad[maxl], thetar[maxl],
    thetas[9], thts[maxl], tstbd[10][10], VarPar[61],
    VolNh4NContent[maxl][maxk], VolNo3NContent[maxl][maxk],
    VolUreaNContent[maxl][maxk], VolWaterContent[maxl][maxk], WindSpeed[24],
    wk[maxk];
extern double light_intercept_parameter, light_intercept_parameters[20];
void WriteStateVariables(bool bAdjusting);
double TotalLeafWeight();           // total leaf weight, g per plant.
double TotalLeafArea();