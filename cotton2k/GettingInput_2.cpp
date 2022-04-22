//  File GettingInput_2.cpp
//     Consists of the following functions:
// InitializeSoilData()
// ReadSoilHydraulicData()
// InitializeSoilTemperature()
// form()
//

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#include <math.h>
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//////////////////////////////////////////////////////////////////////
double form(double c0, double d0, double g0)
//     This function computes the aggregation factor for 2 mixed soil materials.
//  Arguments referenced:
//     c0	- heat conductivity of first material
//     d0	- heat conductivity of second material
//     g0	- shape factor for these materials
//
{
    double form0 = (2 / (1 + (c0 / d0 - 1) * g0) +
                    1 / (1 + (c0 / d0 - 1) * (1 - 2 * g0))) /
                   3;
    return form0;
}
