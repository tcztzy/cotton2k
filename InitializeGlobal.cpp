// File InitializeGlobal.cpp
//
#include "global.h"

///////////////////////////////////////////////////////////////////////////
void InitializeGlobal()
//     This function initializes many "global" variables at the start of a
//  simulation. It is called from ReadInput(). Note that initialization
//  is needed at the start of each simulation (NOT at start of the run).
{
    addwtbl = 0;
    AppliedWater = 0;
    AverageLwp = 0;
    AverageLwpMin = 0;

    BloomWeightLoss = 0;
    BurrNConc = 0;
    BurrNitrogen = 0;
    BurrWeightGreenBolls = 0;
    BurrWeightOpenBolls = 0;

    CarbonAllocatedForRootGrowth = 0;
    CottonWeightGreenBolls = 0;
    CottonWeightOpenBolls = 0;
    CarbonStress = 1;
    CumEvaporation = 0;
    CumFertilizerN = 0;
    CumNetPhotosynth = 0;
    CumNitrogenUptake = 0;
    CumPlantNLoss = 0;
    CumTranspiration = 0;
    CumWaterAdded = 0;
    CumWaterDrained = 0;

    ExtraCarbon = 0;
    FruitGrowthRatio = 1;

    ginp = 0.35;
    Gintot = 0.35;
    GreenBollsLost = 0;

    LastIrrigation = 0;
    LeafAreaIndex = 0.001;
    LeafNConc = .056;
    LeafNitrogen = 0.0112;
    LintYield = 0;

    MaxIrrigation = 0;
    MineralizedOrganicN = 0;

    NitrogenStress = 1;
    NumAbscisedLeaves = 0;
    NumOpenBolls = 0;
    NumPreFruNodes = 1;
    NumSheddingTags = 0;
    NumVegBranches = 1;
    NumWaterTableData = 0;
    NStressFruiting = 1;
    NStressRoots = 1;
    NStressVeg = 1;

    PercentDefoliation = 0;
    PetioleNConc = 0;
    PetioleNitrogen = 0;
    PetioleNO3NConc = 0;
    pixcon = 0;
    pixda = 1;
    pixdn = 1;
    pixdz = 1;
    PixInPlants = 0;
    PlantWeight = 0;
    PotGroStem = 0;

    ReserveC = 0.06;
    RootNConc = 0.026;
    RootNitrogen = 0.0052;
    RootWeightLoss = 0;

    SeedNConc = 0;
    SeedNitrogen = 0;
    SoilNitrogenLoss = 0;
    SquareNConc = 0;
    SquareNitrogen = 0;
    StemNConc = 0.036;
    StemNitrogen = 0.0072;
    SumNO3N90 = 0;
    SupplyNH4N = 0;
    SupplyNO3N = 0;

    TotalActualLeafGrowth = 0;
    TotalActualPetioleGrowth = 0;
    TotalLeafArea = 0;
    TotalLeafWeight = 0.20;
    TotalPetioleWeight = 0;
    TotalRequiredN = 0;
    TotalSquareWeight = 0;
    TotalStemWeight = 0.2;

    WaterStressStem = 1;
    WaterTableLayer = 1000;
//
    for (int i = 0; i < 3; i++) {
        DelayNewFruBranch[i] = 0;
        LwpMinX[i] = 0;
        LwpX[i] = 0;
        NumFruitBranches[i] = 0;
        for (int j = 0; j < 30; j++) {
            DelayNewNode[i][j] = 0;
            LeafAreaMainStem[i][j] = 0;
            LeafWeightMainStem[i][j] = 0;
            NumNodes[i][j] = 0;
            PetioleWeightMainStem[i][j] = 0;
            PotGroLeafAreaMainStem[i][j] = 0;
            PotGroLeafWeightMainStem[i][j] = 0;
            PotGroPetioleWeightMainStem[i][j] = 0;
            for (int k = 0; k < 5; k++) {
                AgeOfBoll[i][j][k] = 0;
                AgeOfSite[i][j][k] = 0;
                AvrgNodeTemper[i][j][k] = 0;
                BollWeight[i][j][k] = 0;
                BurrWeight[i][j][k] = 0;
                FruitingCode[i][j][k] = 0;
                FruitFraction[i][j][k] = 0;
                LeafAge[i][j][k] = 0;
                LeafAreaNodes[i][j][k] = 0;
                LeafWeightNodes[i][j][k] = 0;
                PetioleWeightNodes[i][j][k] = 0;
                PotGroBolls[i][j][k] = 0;
                PotGroBurrs[i][j][k] = 0;
                PotGroLeafAreaNodes[i][j][k] = 0;
                PotGroLeafWeightNodes[i][j][k] = 0;
                PotGroPetioleWeightNodes[i][j][k] = 0;
                PotGroSquares[i][j][k] = 0;
                SquareWeight[i][j][k] = 0;
            }
        }
    }
//
    for (int i = 0; i < 5; i++) {
        DefoliationDate[i] = 0;
        DefoliationMethod[i] = 0;
        DefoliantAppRate[i] = 0;
    }
//
    for (int i = 0; i < 9; i++) {
        AgeOfPreFruNode[i] = 0;
        LeafAreaPreFru[i] = 0;
        PotGroLeafAreaPreFru[i] = 0;
        PotGroLeafWeightPreFru[i] = 0;
        PotGroPetioleWeightPreFru[i] = 0;
        LeafWeightPreFru[i] = 0;
        PetioleWeightPreFru[i] = 0;
    }
//
    for (int i = 0; i < 20; i++) {
        AbscissionLag[i] = 0;
        DayWaterTableInput[i] = 0;
        ElCondSatSoil[i] = 0;
        LevelsOfWaterTable[i] = 0;
        ShedByCarbonStress[i] = 0;
        ShedByNitrogenStress[i] = 0;
        ShedByWaterStress[i] = 0;
    }
//
    for (int k = 0; k < maxk; k++) {
        FoliageTemp[k] = 295;
        MulchTemp[k] = 295;
    }
//
    for (int l = 0; l < maxl; l++) {
        rlat1[l] = 0;
        rlat2[l] = 0;
        for (int k = 0; k < maxk; k++) {
            RootImpede[l][k] = 0;
        }
    }
//
    for (int i = 0; i < 365; i++) {
        StemWeight[i] = 0;
    }
}
//////////////////////////////////////////////////////////
