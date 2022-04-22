//  PlantGrowth_1.cpp
//
//   functions in this file:
// LeafResistance()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
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
