// File   PlantAdjustment.cpp
//
//   List of functions in this file:
//       WriteStateVariables()
//       GoBack()
//
#include "global.h"
#include "GeneralFunctions.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
void WriteStateVariables(bool bAdjusting, const string &Date, const int &Daynum, const int &DayOfSimulation,
                         const int &FirstBloom, const int &FirstSquare, const int &NumLayersWithRoots,
                         const double &PlantHeight, const double &AbscisedFruitSites, const double &AbscisedLeafWeight,
                         const double &WaterStress, const ClimateStruct Clim[400])
//     This function stores all state or rate variables, needed for output, in the structure Scratch21. It is called from DailySimulation() and DailyOutput().
//     Each record is an array cell, starting from day of start of simulation. 
//
{
//     Variables needed for output:
    Scratch21[DayOfSimulation - 1].kday = Kday;
//
    Scratch21[DayOfSimulation - 1].abscisedFruitSites = AbscisedFruitSites;
    Scratch21[DayOfSimulation - 1].averageSoilPsi = AverageSoilPsi;
    Scratch21[DayOfSimulation - 1].avrgDailyTemp = AvrgDailyTemp;
    Scratch21[DayOfSimulation - 1].burrNConc = BurrNConc;
    Scratch21[DayOfSimulation - 1].burrWeightOpenBolls = BurrWeightOpenBolls;
    Scratch21[DayOfSimulation - 1].carbonStress = CarbonStress;
    Scratch21[DayOfSimulation - 1].cottonWeightGreenBolls = CottonWeightGreenBolls;
    Scratch21[DayOfSimulation - 1].cottonWeightOpenBolls = CottonWeightOpenBolls;
    Scratch21[DayOfSimulation - 1].cumEvaporation = CumEvaporation;
    Scratch21[DayOfSimulation - 1].cumFertilizerN = CumFertilizerN;
    Scratch21[DayOfSimulation - 1].cumNetPhotosynth = CumNetPhotosynth;
    Scratch21[DayOfSimulation - 1].cumNitrogenUptake = CumNitrogenUptake;
    Scratch21[DayOfSimulation - 1].cumTranspiration = CumTranspiration;
    Scratch21[DayOfSimulation - 1].cumWaterAdded = CumWaterAdded;
    Scratch21[DayOfSimulation - 1].cumWaterDrained = CumWaterDrained;
    Scratch21[DayOfSimulation - 1].deadwt = AbscisedLeafWeight + BloomWeightLoss + GreenBollsLost + RootWeightLoss;
    Scratch21[DayOfSimulation - 1].date = Date;
    Scratch21[DayOfSimulation - 1].daynum = Daynum;
    Scratch21[DayOfSimulation - 1].dayTimeTemp = DayTimeTemp;
    Scratch21[DayOfSimulation - 1].gbw = CottonWeightGreenBolls + BurrWeightGreenBolls;
//     h2obal is computed as the water balance in mm. It should always be zero.
//  The "positive" amount is the initial water in the soil slab, plus water
//  added by rain and irrigation, and also water added from the water-table.
//  The "negative" is the present total soil water in the soil slab, and cumulative
//  amounts lost by transpiration, evaporation and drainage.
    Scratch21[DayOfSimulation - 1].h2obal = InitialTotalSoilWater + CumWaterAdded
                                            + addwtbl - TotalSoilWater - CumTranspiration
                                            - CumEvaporation - CumWaterDrained;
    Scratch21[DayOfSimulation - 1].leafAreaIndex = LeafAreaIndex;
    Scratch21[DayOfSimulation - 1].leafNConc = LeafNConc;
    Scratch21[DayOfSimulation - 1].lightIntercept = LightIntercept;
    Scratch21[DayOfSimulation - 1].lintYield = LintYield;
    Scratch21[DayOfSimulation - 1].lwpMin = LwpMin;
    Scratch21[DayOfSimulation - 1].mainStemNodes = MainStemNodes;
    Scratch21[DayOfSimulation - 1].mineralizedOrganicN = MineralizedOrganicN;
    Scratch21[DayOfSimulation - 1].netPhotosynthesis = NetPhotosynthesis;
    Scratch21[DayOfSimulation - 1].nightTimeTemp = NightTimeTemp;
    Scratch21[DayOfSimulation - 1].nitrogenStress = NitrogenStress;
    Scratch21[DayOfSimulation - 1].nStressFruiting = NStressFruiting;
    Scratch21[DayOfSimulation - 1].nStressVeg = NStressVeg;
    Scratch21[DayOfSimulation - 1].numGreenBolls = NumGreenBolls;
    Scratch21[DayOfSimulation - 1].numOpenBolls = NumOpenBolls;
    Scratch21[DayOfSimulation - 1].numSquares = NumSquares;
    Scratch21[DayOfSimulation - 1].petioleNConc = PetioleNConc;
    Scratch21[DayOfSimulation - 1].petioleNO3NConc = PetioleNO3NConc;
    Scratch21[DayOfSimulation - 1].plantHeight = PlantHeight;
    Scratch21[DayOfSimulation - 1].plantWeight = PlantWeight;
    Scratch21[DayOfSimulation - 1].rad = GetFromClim(Clim, "rad", Daynum);
    Scratch21[DayOfSimulation - 1].rain = GetFromClim(Clim, "rain", Daynum);
    Scratch21[DayOfSimulation - 1].reserveC = ReserveC;
    Scratch21[DayOfSimulation - 1].rn = Rn;
    Scratch21[DayOfSimulation - 1].rootNConc = RootNConc;
    Scratch21[DayOfSimulation - 1].seedNConc = SeedNConc;
    Scratch21[DayOfSimulation - 1].soilNitrogenLoss = SoilNitrogenLoss;
    Scratch21[DayOfSimulation - 1].stemNConc = StemNConc;
    Scratch21[DayOfSimulation - 1].sumNO3N90 = SumNO3N90;
    Scratch21[DayOfSimulation - 1].tmax = GetFromClim(Clim, "tmax", Daynum);
    Scratch21[DayOfSimulation - 1].tmin = GetFromClim(Clim, "tmin", Daynum);
    Scratch21[DayOfSimulation - 1].totalLeafWeight = TotalLeafWeight;
    Scratch21[DayOfSimulation - 1].totalPetioleWeight = TotalPetioleWeight;
    Scratch21[DayOfSimulation - 1].totalRootWeight = TotalRootWeight;
    Scratch21[DayOfSimulation - 1].totalSquareWeight = TotalSquareWeight;
    Scratch21[DayOfSimulation - 1].totalSoilWater = TotalSoilWater;
    Scratch21[DayOfSimulation - 1].totalStemWeight = TotalStemWeight;
    Scratch21[DayOfSimulation - 1].waterStress = WaterStress;
    Scratch21[DayOfSimulation - 1].waterStressStem = WaterStressStem;
    Scratch21[DayOfSimulation - 1].wind = GetFromClim(Clim, "wind", Daynum);
//
    for (int l = 0; l < maxl; l++)
        for (int k = 0; k < maxk; k++) {
            Scratch21[DayOfSimulation - 1].rootWtCapblUptake[l][k] = RootWtCapblUptake[l][k];
            Scratch21[DayOfSimulation - 1].soilPsi[l][k] = SoilPsi[l][k];
            Scratch21[DayOfSimulation - 1].soilTempDailyAvrg[l][k] = SoilTempDailyAvrg[l][k];
            Scratch21[DayOfSimulation - 1].volNh4NContent[l][k] = VolNh4NContent[l][k];
            Scratch21[DayOfSimulation - 1].volNo3NContent[l][k] = VolNo3NContent[l][k];
            Scratch21[DayOfSimulation - 1].volWaterContent[l][k] = VolWaterContent[l][k];
        }
}
