//  RootGrowth_2.cpp
//
//   functions in this file:
// InitiateLateralRoots()
// LateralRootGrowthLeft()
// LateralRootGrowthRight()
// RootSummation()
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

//////////////////////////////
double RootCultivation(SoilCell soil_cells[40][20], int NumRootAgeGroups, double cultivation_depth, double DailyRootLoss, double row_space)
//     This function is executed on the day of soil cultivation. It is called from
//  ActualRootGrowth(). It has been adapted from GOSSYM. It is assumed that the roots in the
//  upper soil layers, as defined by the depth of cultivation, are destroyed, with the
//  exception of the soil cells that are within 15 cm of the plant row.
//
//     The following global variables are referenced here:
//       dl, CultivationDepth, NumRootAgeGroups, nk, nl, PlantRowLocation, wk
//     The following global variables are set here:
//       RootWeight, DailyRootLoss
//     The argument j - is the serial number of this cultivation.
{
    //     The depth of cultivation (CultivationDepth) is in cm. The number of layers affected by
    //  it (lcult) is determined. Loop for all columns, and check if this column is more than 15 cm
    //  from the plant row: if this is true, destroy all roots and add their weight to DailyRootLoss.
    int lcult = 0;    // number of soil layers affected by cultivation.
    double sdpth = 0; // sum depth to the end of the layer.
    for (int l = 0; l < nl; l++)
    {
        sdpth += dl(l);
        if (sdpth >= cultivation_depth)
        {
            lcult = l;
            break;
        }
    }
    //
    double sumwk = 0; // sum of column widths from edge of slab to this column.
    for (int k = 0; k < nk; k++)
    {
        sumwk += wk(k, row_space);
        if (fabs(sumwk - PlantRowLocation) >= 15)
        {
            for (int l = 0; l < lcult; l++)
                for (int i = 0; i < NumRootAgeGroups; i++)
                {
                    DailyRootLoss += soil_cells[l][k].root.weight[i];
                    soil_cells[l][k].root.weight[i] = 0;
                }
        }
    }
    return DailyRootLoss;
}

//////////////////////////////
void RootSummation(State &state, int NumRootAgeGroups, double row_space, double per_plant_area)
//     This function has been added for compatibility with GOSSYM root routines.
//  It is called from ActualRootGrowth(). It summarizes root data, in a form ready
//  for output or plotting. Sums of root weights for cells, for age groups and for
//  the total slab are calculated. state.root_weight is calculated in g per plant.
//
//     The following global variables are referenced here:
//  dl, Kday, nk, nl, NumRootAgeGroups,
//  RootWeight, RootWtCapblUptake, RowSpace, wk
{
    //     Compute the total root weight (of all age classes) for all soil cells as
    double roots = 0; // total weight of roots of all classes, g per slab.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++)
            roots += accumulate(state.soil.cells[l][k].root.weight, state.soil.cells[l][k].root.weight + NumRootAgeGroups, double(0));
    //     Convert total root weight from g per slab to g per plant.
    state.root_weight = roots * 100 * per_plant_area / row_space;
}
