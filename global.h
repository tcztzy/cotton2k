//   global.h
#pragma once
//
//  definition of global variables
//  ==============================
//    For dictionary of global variables see file "global.cpp"
////    Constants    ////
const int maxl = 40;
const int maxk = 20;
const double pi = 3.14159;
enum CLIMATE_METRIC {
    TMAX = 0,
    TMIN = 1,
    IRRD = 2,
    RAIN = 3,
    WIND = 4,
    TDEW = 5,
};
////    Structures    ////
//
extern struct Climstruct {
    int nDay;
    double Rad, Tmax, Tmin, Rain, Wind, Tdew;
} Clim[400];
extern struct Irrigation {
    int day, method, LocationColumnDrip, LocationLayerDrip;
    double amount;
} Irrig[150];
////    Integers    ////
extern int DayEmerge, DayEndMulch, DayFinish, DayFirstDef,
    DayOfSimulation, Daynum, DayPlant, DayStart, DayStartMulch,
    DayStartPredIrrig, DayStopPredIrrig, FirstBloom, FirstSquare, inrim,
    IrrigMethod, isw, iyear, Kday, KdayAdjust,
    LastIrrigation, LastTaprootLayer, LocationColumnDrip, LocationLayerDrip,
    MainStemNodes, MinDaysBetweenIrrig, MulchIndicator, ncurve, nk, nl, noitr,
    NumAbscisedLeaves, NumAdjustDays, NumFruitSites, NumIrrigations,
    NumLayersWithRoots, NumPreFruNodes, NumRootAgeGroups,
    NumSheddingTags, NumVegBranches, PlantRowColumn,
    WaterTableLayer;
extern int CultivationDate[5], DefoliationDate[5],
    DefoliationMethod[5], FruitingCode[3][30][5], LateralRootFlag[maxl],
    NumFruitBranches[3], NumNodes[3][30], NodeLayer[3][30],
    NodeLayerPreFru[9], RootColNumLeft[maxl],
    RootColNumRight[maxl], SoilHorizonNum[maxl];
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
    conmax, CottonWeightGreenBolls, CottonWeightOpenBolls,
    CumEvaporation, CumFertilizerN, CumNetPhotosynth, CumNitrogenUptake,
    CumPlantNLoss, CumTranspiration, CumWaterAdded, CumWaterDrained,
    DailyRootLoss, DayInc, DayTimeTemp, dclay, DeepSoilTemperature,
    DensityFactor, DepthLastRootLayer, dsand, ElCondSatSoilToday, ExtraCarbon,
    FruitGrowthRatio, ginp, Gintot, GreenBollsLost,
    InitialTotalSoilWater, IrrigationDepth, LeafAreaIndex, LeafNConc,
    LeafNitrogen, LeafWeightAreaRatio, LightIntercept, LintYield,
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
    TotalPetioleWeight, TotalRequiredN, TotalRootWeight, TotalSoilNh4N,
    TotalSoilNitrogen, TotalSoilNo3N, TotalSoilUreaN, TotalSoilWater,
    TotalSquareWeight, TotalStemWeight, WaterStress, WaterStressStem;

extern double AbscissionLag[20], ActualRootGrowth[maxl][maxk],
    AgeOfBoll[3][30][5], AgeOfPreFruNode[9], AgeOfSite[3][30][5], airdr[9],
    AirTemp[24], albedo[24], alpha[9], AvrgNodeTemper[3][30][5],
    AverageLeafAge[20], beta[9], BollWeight[3][30][5], BulkDensity[9],
    BurrWeight[3][30][5], cgind[3], ClayVolumeFraction[maxl],
    CloudCoverRatio[24], CloudTypeCorr[24], CultivationDepth[5],
    DefoliantAppRate[5], DelayNewFruBranch[3], DelayNewNode[3][30],
    DewPointTemp[24], dl[maxl], es1hour[24], es2hour[24],
    FieldCapacity[maxl], FoliageTemp[maxk], FreshOrganicMatter[maxl][maxk],
    FreshOrganicNitrogen[maxl][maxk], FruitFraction[3][30][5], gh2oc[10],
    HeatCapacitySoilSolid[maxl], HeatCondDrySoil[maxl],
    HumusNitrogen[maxl][maxk], HumusOrganicMatter[maxl][maxk], impede[10][10],
    LeafAge[3][30][5], LeafAreaMainStem[3][30], LeafArea[20],
    LeafAreaIndexes[20], LeafAreaNodes[3][30][5], LeafAreaPreFru[9],
    LeafWeightLayer[20], LeafWeightMainStem[3][30], LeafWeightNodes[3][30][5],
    LeafWeightPreFru[9], LightInterceptLayer[20],
    LwpMinX[3], LwpX[3], MarginalWaterContent[maxl],
    MaxWaterCapacity[maxl], MulchTemp[maxk],
    NO3FlowFraction[maxl], PetioleWeightMainStem[3][30],
    PetioleWeightNodes[3][30][5], PetioleWeightPreFru[9],
    PoreSpace[maxl], PotGroBolls[3][30][5], PotGroBurrs[3][30][5],
    PotGroLeafAreaMainStem[3][30], PotGroLeafAreaNodes[3][30][5],
    PotGroLeafAreaPreFru[9], PotGroLeafWeightMainStem[3][30],
    PotGroLeafWeightNodes[3][30][5], PotGroLeafWeightPreFru[9],
    PotGroPetioleWeightMainStem[3][30], PotGroPetioleWeightNodes[3][30][5],
    PotGroPetioleWeightPreFru[9], PotGroRoots[maxl][maxk],
    PotGroSquares[3][30][5], Radiation[24], ReferenceETP[24],
    RelativeHumidity[24], rlat1[maxl], rlat2[maxl], RootAge[maxl][maxk],
    RootGroFactor[maxl][maxk], RootImpede[maxl][maxk],
    RootWeight[maxl][maxk][3], RootWtCapblUptake[maxl][maxk],
    SandVolumeFraction[maxl], SaturatedHydCond[9], ShedByCarbonStress[20],
    ShedByNitrogenStress[20], ShedByWaterStress[20],
    SoilPsi[maxl][maxk], SoilTemp[maxl][maxk], SoilTempDailyAvrg[maxl][maxk],
    SquareWeight[3][30][5], StemWeight[365], thad[maxl], thetar[maxl],
    thetas[9], thts[maxl], tstbd[10][10], VarPar[61],
    VolNh4NContent[maxl][maxk], VolNo3NContent[maxl][maxk],
    VolUreaNContent[maxl][maxk], VolWaterContent[maxl][maxk], WindSpeed[24],
    wk[maxk];
void WriteStateVariables(bool bAdjusting);
double TotalLeafWeight();  // total leaf weight, g per plant.
double TotalLeafArea();
