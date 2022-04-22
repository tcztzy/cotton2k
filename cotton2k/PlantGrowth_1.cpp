//  PlantGrowth_1.cpp
//
//   functions in this file:
// PhysiologicalAge()
// LeafResistance()
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
