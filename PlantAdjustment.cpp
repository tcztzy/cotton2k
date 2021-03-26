// File   PlantAdjustment.cpp
//
//   List of functions in this file:
//       WriteStateVariables()
//       GoBack()
//
#include "global.h"
#include "Simulation.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
void WriteStateVariables(Simulation &sim, unsigned int u)
//     This function stores all state or rate variables, needed for output, in the structure Scratch21. It is called from DailySimulation().
//     Each record is an array cell, starting from day of start of simulation.
//
{
    State &state = sim.states[u];
//     Variables needed for output:
    Scratch21[u].kday = Kday;
//
    Scratch21[u].averageSoilPsi = AverageSoilPsi;
    Scratch21[u].avrgDailyTemp = AvrgDailyTemp;
    Scratch21[u].burrNConc = BurrNConc;
    Scratch21[u].burrWeightOpenBolls = BurrWeightOpenBolls;
    Scratch21[u].cottonWeightGreenBolls = CottonWeightGreenBolls;
    Scratch21[u].cottonWeightOpenBolls = CottonWeightOpenBolls;
    Scratch21[u].cumFertilizerN = CumFertilizerN;
    Scratch21[u].cumNetPhotosynth = CumNetPhotosynth;
    Scratch21[u].cumNitrogenUptake = CumNitrogenUptake;
    Scratch21[u].cumWaterAdded = CumWaterAdded;
    Scratch21[u].cumWaterDrained = CumWaterDrained;
    Scratch21[u].deadwt = sim.states[u].abscised_leaf_weight + state.bloom_weight_loss + GreenBollsLost + RootWeightLoss;
    Scratch21[u].dayTimeTemp = DayTimeTemp;
    Scratch21[u].gbw = CottonWeightGreenBolls + BurrWeightGreenBolls;
    Scratch21[u].lightIntercept = LightIntercept;
    Scratch21[u].lwpMin = LwpMin;
    Scratch21[u].mainStemNodes = MainStemNodes;
    Scratch21[u].mineralizedOrganicN = MineralizedOrganicN;
    Scratch21[u].netPhotosynthesis = NetPhotosynthesis;
    Scratch21[u].nightTimeTemp = NightTimeTemp;
    Scratch21[u].nStressFruiting = NStressFruiting;
    Scratch21[u].nStressVeg = NStressVeg;
    Scratch21[u].petioleNO3NConc = PetioleNO3NConc;
    Scratch21[u].reserveC = ReserveC;
    Scratch21[u].soilNitrogenLoss = SoilNitrogenLoss;
    Scratch21[u].stemNConc = StemNConc;
    Scratch21[u].sumNO3N90 = SumNO3N90;
    Scratch21[u].totalLeafWeight = TotalLeafWeight;
    Scratch21[u].totalPetioleWeight = TotalPetioleWeight;
    Scratch21[u].totalRootWeight = TotalRootWeight;
    Scratch21[u].totalSquareWeight = TotalSquareWeight;
    Scratch21[u].totalSoilWater = TotalSoilWater;
    Scratch21[u].totalStemWeight = TotalStemWeight;
//
    for (int l = 0; l < maxl; l++)
        for (int k = 0; k < maxk; k++) {
            Scratch21[u].soilPsi[l][k] = SoilPsi[l][k];
            Scratch21[u].soilTempDailyAvrg[l][k] = SoilTempDailyAvrg[l][k];
            Scratch21[u].volNh4NContent[l][k] = VolNh4NContent[l][k];
            Scratch21[u].volWaterContent[l][k] = VolWaterContent[l][k];
        }
}
