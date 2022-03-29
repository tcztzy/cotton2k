//  PlantGrowth_1.cpp
//
//   functions in this file:
// PhysiologicalAge()
// LeafResistance()
// PlantGrowth()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//////////////////////////////////////////////////
double PhysiologicalAge()  // computes physiological age
//     This function returns the daily 'physiological age' increment,
//  based on hourly temperatures. It is called each day by SimulateThisDay().
//     The following global variable is used here:
//        AirTemp[] = array of hourly temperatures.
{
    //     The following constant Parameters are used in this function:
    const double p1 = 12;  // threshold temperature, C
    const double p2 =
        14;  // temperature, C, above p1, for one physiological day.
    const double p3 = 1.5;  // maximum value of a physiological day.
    //     The threshold value is assumed to be 12 C (p1). One physiological day
    //     is
    //  equivalent to a day with an average temperature of 26 C, and therefore
    //  the heat units are divided by 14 (p2).
    //     A linear relationship is assumed between temperature and heat unit
    //  accumulation in the range of 12 C (p1) to 33 C (p2*p3+p1). the effect of
    //  temperatures higher than 33 C is assumed to be equivalent to that of 33
    //  C.
    double dayfd =
        0;  // the daily contribution to physiological age (return value).
    for (int ihr = 0; ihr < 24; ihr++) {
        double tfd = (AirTemp[ihr] - p1) /
                     p2;  // the hourly contribution to physiological age.
        if (tfd < 0) tfd = 0;
        if (tfd > p3) tfd = p3;
        dayfd += tfd;
    }
    return dayfd / 24;
}
//////////////////////////////////
double LeafResistance(double agel)
//     This function computes and returns the resistance of leaves of cotton
// plants to transpiration. It is assumed to be a function of leaf age.
// It is called from LeafWaterPotential().
//     The input argument (agel) is leaf age in physiological days.
{
    //     The following constant parameters are used:
    const double afac = 160;   // factor used for computing leaf resistance.
    const double agehi = 94;   // higher limit for leaf age.
    const double agelo = 48;   // lower limit for leaf age.
    const double rlmin = 0.5;  // minimum leaf resistance.
                               //
    double leafResistance;
    if (agel <= agelo)
        leafResistance = rlmin;
    else if (agel >= agehi)
        leafResistance = rlmin + (agehi - agelo) * (agehi - agelo) / afac;
    else {
        double ax = 2 * agehi - agelo;  // intermediate variable
        leafResistance = rlmin + (agel - agelo) * (ax - agel) / afac;
    }
    return leafResistance;
}
////////////////////////////////////////////////////////////////////////////
void PlantGrowth()
//     This function simulates the potential and actual growth of cotton plants.
//  It is called from SimulateThisDay(), and it calls the following functions:
//    ActualFruitGrowth(), ActualLeafGrowth(), ActualRootGrowth(),
//    AddPlantHeight(), DryMatterBalance(), PotentialFruitGrowth(),
//    PotentialLeafGrowth(), PotentialRootGrowth(), PotentialStemGrowth().
//
//     The following global variables are referenced here:
//        ActualStemGrowth, DayInc, FirstSquare, FruitingCode, Kday, pixdz,
//        PerPlantArea, RowSpace, WaterStressStem.
//
//     The following global variables are set here:
//        LeafAreaIndex, PlantHeight, PotGroAllRoots, PotGroStem, StemWeight,
//        TotalLeafArea, TotalLeafWeight, TotalPetioleWeight, TotalStemWeight.
{
    //     Call PotentialLeafGrowth() to compute potential growth rate of
    //     leaves.
    PotentialLeafGrowth();
    //     If it is after first square, call PotentialFruitGrowth() to compute
    //     potential
    //  growth rate of squares and bolls.
    if (FruitingCode[0][0][0] > 0) PotentialFruitGrowth();
    //     Active stem tissue (stemnew) is the difference between
    //     TotalStemWeight
    //  and the value of StemWeight(kkday).
    int voldstm =
        32;  // constant parameter (days for stem tissue to become "old")
    int kkday = Kday - voldstm;  // age of young stem tissue
    if (kkday < 1) kkday = 1;
    double stemnew = TotalStemWeight -
                     StemWeight[kkday];  // dry weight of active stem tissue.
    //     Call PotentialStemGrowth() to compute PotGroStem, potential growth
    //     rate of stems.
    //  The effect of temperature is introduced, by multiplying potential growth
    //  rate by DayInc. Stem growth is also affected by water stress
    //  (WaterStressStem) and possible PIX application (pixdz).   PotGroStem is
    //  limited by (maxstmgr * PerPlantArea) g per plant per day.
    PotGroStem =
        PotentialStemGrowth(stemnew) * DayInc * WaterStressStem * pixdz;
    double maxstmgr =
        0.067;  // maximum posible potential stem growth, g dm-2 day-1.
    if (PotGroStem > maxstmgr * PerPlantArea)
        PotGroStem = maxstmgr * PerPlantArea;
    //	   Call PotentialRootGrowth() to compute potential growth rate of roots.
    double sumpdr;  // total potential growth rate of roots in g per slab. this
                    // is computed in PotentialRootGrowth() and used in
                    // ActualRootGrowth().
    sumpdr = PotentialRootGrowth();
    //     Total potential growth rate of roots is converted from g per
    //  slab (sumpdr) to g per plant (PotGroAllRoots).
    PotGroAllRoots = sumpdr * 100 * PerPlantArea / RowSpace;
    //     Limit PotGroAllRoots to (maxrtgr*PerPlantArea) g per plant per day.
    double maxrtgr =
        0.045;  // maximum possible potential root growth, g dm-2 day-1.
    if (PotGroAllRoots > maxrtgr * PerPlantArea)
        PotGroAllRoots = maxrtgr * PerPlantArea;
    //     Call DryMatterBalance() to compute carbon balance, allocation of
    //     carbon to
    //  plant parts, and carbon stress. DryMatterBalance() also computes and
    //  returns the values of the following arguments:
    //     cdleaf is carbohydrate requirement for leaf growth, g per plant per
    //     day. cdpet is carbohydrate requirement for petiole growth, g per
    //     plant per day. cdroot is carbohydrate requirement for root growth, g
    //     per plant per day. cdstem is carbohydrate requirement for stem
    //     growth, g per plant per day.
    double cdstem, cdleaf, cdpet, cdroot;
    DryMatterBalance(cdstem, cdleaf, cdpet, cdroot);
    //     If it is after first square, call ActualFruitGrowth() to compute
    //     actual
    //  growth rate of squares and bolls.
    if (FruitingCode[0][0][0] > 0) ActualFruitGrowth();
    //     Initialize TotalLeafWeight. It is assumed that cotyledons fall off
    //  at time of first square. Also initialize TotalLeafArea and
    //  TotalPetioleWeight.
    TotalPetioleWeight = 0;
    //     Call ActualLeafGrowth to compute actual growth rate of leaves and
    //     compute leaf area index.
    ActualLeafGrowth();
    LeafAreaIndex = TotalLeafArea() / PerPlantArea;
    //     Add ActualStemGrowth to TotalStemWeight, and define StemWeight(Kday)
    //     for this day.
    TotalStemWeight += ActualStemGrowth;
    StemWeight[Kday] = TotalStemWeight;
    //     Plant density affects growth in height of tall plants.
    double htdenf = 55;  // minimum plant height for plant density affecting
                         // growth in height.
    double z1;           // intermediate variable to compute denf2.
    z1 = (PlantHeight - htdenf) / htdenf;
    if (z1 < 0) z1 = 0;
    if (z1 > 1) z1 = 1;
    double denf2;  // effect of plant density on plant growth in height.
    denf2 = 1 + z1 * (DensityFactor - 1);
    //     Call AddPlantHeight to compute PlantHeight.
    PlantHeight += AddPlantHeight(denf2);
    //     Call ActualRootGrowth() to compute actual root growth.
    ComputeActualRootGrowth(sumpdr);
}
