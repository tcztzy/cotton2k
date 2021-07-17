// File InitializeGlobal.cpp
//
#include "global.h"

///////////////////////////////////////////////////////////////////////////
void InitializeGlobal()
//     This function initializes many "global" variables at the start of a
//  simulation. It is called from ReadInput(). Note that initialization
//  is needed at the start of each simulation (NOT at start of the run).
{
    AverageLwp = 0;
    AverageLwpMin = 0;

    BurrNConc = 0;

    CumFertilizerN = 0;
    CumNetPhotosynth = 0;
    CumNitrogenUptake = 0;
    CumWaterAdded = 0;
    CumWaterDrained = 0;

    GreenBollsLost = 0;

    LastIrrigation = 0;

    MaxIrrigation = 0;
    MineralizedOrganicN = 0;

    NumAbscisedLeaves = 0;
    NumSheddingTags = 0;
    NumWaterTableData = 0;

    PercentDefoliation = 0;
    PetioleNitrogen = 0;
    PotGroStem = 0;

    ReserveC = 0.06;
    RootNitrogen = 0.0052;
    RootWeightLoss = 0;

    SoilNitrogenLoss = 0;
    StemNConc = 0.036;
    SumNO3N90 = 0;

    WaterTableLayer = 1000;
//
    for (int i = 0; i < 3; i++) {
        LwpMinX[i] = 0;
        LwpX[i] = 0;
    }
//
    for (int i = 0; i < 5; i++) {
        DefoliationDate[i] = 0;
        DefoliationMethod[i] = 0;
        DefoliantAppRate[i] = 0;
    }
//
    for (int i = 0; i < 9; i++) {
        PotGroLeafAreaPreFru[i] = 0;
        PotGroLeafWeightPreFru[i] = 0;
        PotGroPetioleWeightPreFru[i] = 0;
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
    }
//
    for (int l = 0; l < maxl; l++) {
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
