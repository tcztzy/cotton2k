//  RootGrowth_2.cpp
//
//   functions in this file:
// InitiateLateralRoots()
// LateralRootGrowthLeft()
// LateralRootGrowthRight()
//
#include <cstdint>
#include <cmath>
#include <numeric>
#include "global.h"
#include "Simulation.hpp"

using namespace std;

extern "C"
{
    double SoilTemOnRootGrowth(double);
}

//////////////////////////////
void InitiateLateralRoots()
//     This function initiates lateral root growth. It is called from ComputeActualRootGrowth().
//
//     The following global variables are referenced here:
//       DepthLastRootLayer, dl, LastTaprootLayer, TapRootLength
//     The following global variable is set here:      LateralRootFlag
{
    const double distlr = 12; // the minimum distance, in cm, from the tip of the taproot, for a lateral root to be able to grow.
    double sdl;               // distance of a layer from tip of taproot, cm.
    sdl = TapRootLength - DepthLastRootLayer;
    //     Loop on soil layers, from the lowest layer with roots upward:
    for (int l = LastTaprootLayer; l >= 0; l--)
    {
        //     Compute distance from tip of taproot.
        sdl += dl(l);
        //     If a layer is marked for a lateral (LateralRootFlag[l] = 1) and its
        //  distance from the tip is larger than distlr - initiate a lateral
        //  (LateralRootFlag[l] = 2).
        if (sdl > distlr && LateralRootFlag[l] == 1)
            LateralRootFlag[l] = 2;
    }
}
