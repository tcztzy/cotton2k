// File   PlantAdjustment.cpp
//
//   List of functions in this file:
//       WriteStateVariables()
//       GoBack()
//
#include "global.h"
#include "GeneralFunctions.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
void WriteStateVariables(bool bAdjusting, const string &Date, const int &Daynum, const int &u,
                         const int &FirstBloom, const int &FirstSquare, const int &NumLayersWithRoots,
                         const double &PlantHeight, const double &AbscisedFruitSites, const double &AbscisedLeafWeight,
                         const double &WaterStress, const ClimateStruct Clim[400])
//     This function stores all state or rate variables, needed for output, in the structure Scratch21. It is called from DailySimulation() and DailyOutput().
//     Each record is an array cell, starting from day of start of simulation. 
//
{
//     Variables needed for output:
    Scratch21[u].kday = Kday;
//
    Scratch21[u].abscisedFruitSites = AbscisedFruitSites;
    Scratch21[u].averageSoilPsi = AverageSoilPsi;
    Scratch21[u].avrgDailyTemp = AvrgDailyTemp;
    Scratch21[u].burrNConc = BurrNConc;
    Scratch21[u].burrWeightOpenBolls = BurrWeightOpenBolls;
    Scratch21[u].carbonStress = CarbonStress;
    Scratch21[u].cottonWeightGreenBolls = CottonWeightGreenBolls;
    Scratch21[u].cottonWeightOpenBolls = CottonWeightOpenBolls;
    Scratch21[u].cumEvaporation = CumEvaporation;
    Scratch21[u].cumFertilizerN = CumFertilizerN;
    Scratch21[u].cumNetPhotosynth = CumNetPhotosynth;
    Scratch21[u].cumNitrogenUptake = CumNitrogenUptake;
    Scratch21[u].cumTranspiration = CumTranspiration;
    Scratch21[u].cumWaterAdded = CumWaterAdded;
    Scratch21[u].cumWaterDrained = CumWaterDrained;
    Scratch21[u].deadwt = AbscisedLeafWeight + BloomWeightLoss + GreenBollsLost + RootWeightLoss;
    Scratch21[u].date = Date;
    Scratch21[u].daynum = Daynum;
    Scratch21[u].dayTimeTemp = DayTimeTemp;
    Scratch21[u].gbw = CottonWeightGreenBolls + BurrWeightGreenBolls;
//     h2obal is computed as the water balance in mm. It should always be zero.
//  The "positive" amount is the initial water in the soil slab, plus water
//  added by rain and irrigation, and also water added from the water-table.
//  The "negative" is the present total soil water in the soil slab, and cumulative
//  amounts lost by transpiration, evaporation and drainage.
    Scratch21[u].h2obal = InitialTotalSoilWater + CumWaterAdded
                                            + addwtbl - TotalSoilWater - CumTranspiration
                                            - CumEvaporation - CumWaterDrained;
    Scratch21[u].leafAreaIndex = LeafAreaIndex;
    Scratch21[u].leafNConc = LeafNConc;
    Scratch21[u].lightIntercept = LightIntercept;
    Scratch21[u].lintYield = LintYield;
    Scratch21[u].lwpMin = LwpMin;
    Scratch21[u].mainStemNodes = MainStemNodes;
    Scratch21[u].mineralizedOrganicN = MineralizedOrganicN;
    Scratch21[u].netPhotosynthesis = NetPhotosynthesis;
    Scratch21[u].nightTimeTemp = NightTimeTemp;
    Scratch21[u].nitrogenStress = NitrogenStress;
    Scratch21[u].nStressFruiting = NStressFruiting;
    Scratch21[u].nStressVeg = NStressVeg;
    Scratch21[u].numGreenBolls = NumGreenBolls;
    Scratch21[u].numOpenBolls = NumOpenBolls;
    Scratch21[u].numSquares = NumSquares;
    Scratch21[u].petioleNConc = PetioleNConc;
    Scratch21[u].petioleNO3NConc = PetioleNO3NConc;
    Scratch21[u].plantHeight = PlantHeight;
    Scratch21[u].plantWeight = PlantWeight;
    Scratch21[u].rad = Clim[u].Rad;
    Scratch21[u].rain = Clim[u].Rain;
    Scratch21[u].reserveC = ReserveC;
    Scratch21[u].rn = Rn;
    Scratch21[u].rootNConc = RootNConc;
    Scratch21[u].seedNConc = SeedNConc;
    Scratch21[u].soilNitrogenLoss = SoilNitrogenLoss;
    Scratch21[u].stemNConc = StemNConc;
    Scratch21[u].sumNO3N90 = SumNO3N90;
    Scratch21[u].tmax = Clim[u].Tmax;
    Scratch21[u].tmin = Clim[u].Tmin;
    Scratch21[u].totalLeafWeight = TotalLeafWeight;
    Scratch21[u].totalPetioleWeight = TotalPetioleWeight;
    Scratch21[u].totalRootWeight = TotalRootWeight;
    Scratch21[u].totalSquareWeight = TotalSquareWeight;
    Scratch21[u].totalSoilWater = TotalSoilWater;
    Scratch21[u].totalStemWeight = TotalStemWeight;
    Scratch21[u].waterStress = WaterStress;
    Scratch21[u].waterStressStem = WaterStressStem;
    Scratch21[u].wind = Clim[u].Wind;
//
    for (int l = 0; l < maxl; l++)
        for (int k = 0; k < maxk; k++) {
            Scratch21[u].rootWtCapblUptake[l][k] = RootWtCapblUptake[l][k];
            Scratch21[u].soilPsi[l][k] = SoilPsi[l][k];
            Scratch21[u].soilTempDailyAvrg[l][k] = SoilTempDailyAvrg[l][k];
            Scratch21[u].volNh4NContent[l][k] = VolNh4NContent[l][k];
            Scratch21[u].volNo3NContent[l][k] = VolNo3NContent[l][k];
            Scratch21[u].volWaterContent[l][k] = VolWaterContent[l][k];
        }
}
