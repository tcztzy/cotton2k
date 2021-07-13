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
void LateralRootGrowthLeft(State &state, int l, int NumRootAgeGroups, unsigned int plant_row_column, double row_space)
//     This function computes the elongation of the lateral roots
//  in a soil layer(l) to the left. It is called from ActualRootGrowth().
//     It calls function SoilTemOnRootGrowth().
//
//     The following global variables are referenced here:
//       NumRootAgeGroups, PlantRowColumn, PoreSpace, RootGroFactor,
//       SoilTempDailyAvrg.
//     The following global variables are set here:
//       RootAge, RootColNumLeft, RootWeight.
//     The argument used:     l - layer number in the soil slab.
{
    //     The following constant parameters are used:
    const double p1 = 0.10;   // constant parameter.
    const double rlatr = 3.6; // potential growth rate of lateral roots, cm/day.
    const double rtran = 0.2; // the ratio of root mass transferred to a new soil
    // soil cell, when a lateral root grows into it.
    //     On its initiation, lateral root length is assumed to be equal to the
    //  width of a soil column soil cell at the location of the taproot.
    if (rlat1[l] <= 0)
        rlat1[l] = wk(plant_row_column, row_space);
    double stday; // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][plant_row_column] - 273.161;
    double temprg; // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
    //     Define the column with the tip of this lateral root (ktip)
    int ktip = 0;     // column with the tips of the laterals to the left
    double sumwk = 0; // summation of columns width
    for (int k = plant_row_column; k >= 0; k--)
    {
        sumwk += wk(k, row_space);
        if (sumwk >= rlat1[l])
        {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the left. Potential
    //  growth rate (u) is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if (state.soil.cells[l][ktip].water_content < PoreSpace[l])
    {
        rlat1[l] += rlatr * temprg * (1 - p1 + state.soil.cells[l][ktip].root.growth_factor * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran) of
        //	mass of roots is transferred to the new soil cell.
        if (rlat1[l] > sumwk && ktip > 0)
        {
            int newktip; // column into which the tip of the lateral grows to left.
            newktip = ktip - 1;
            for (int i = 0; i < NumRootAgeGroups; i++)
            {
                double tran = state.soil.cells[l][ktip].root.weight[i] * rtran;
                state.soil.cells[l][ktip].root.weight[i] -= tran;
                state.soil.cells[l][newktip].root.weight[i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer idefi
            if (state.soil.cells[l][newktip].root.age == 0)
                state.soil.cells[l][newktip].root.age = 0.01;
            if (newktip < state.soil.layers[l].number_of_left_columns_with_root)
                state.soil.layers[l].number_of_left_columns_with_root = newktip;
        }
    }
}

//////////////////////////////
void LateralRootGrowthRight(State &state, int l, int NumRootAgeGroups, unsigned int plant_row_column, double row_space)
//     This function computes the elongation of the lateral roots
//  in a soil layer(l) to the right. It is called from ActualRootGrowth().
//     It calls function SoilTemOnRootGrowth().
//
//     The following global variables are referenced here:
//  nk, NumRootAgeGroups, PlantRowColumn, PoreSpace, RootGroFactor,
//  SoilTempDailyAvrg.
//     The following global variables are set here:
//  RootAge, RootColNumRight, RootWeight.
//     The argument used:      l - layer number in the soil slab.
{
    //     The following constant parameters are used:
    const double p1 = 0.10;   // constant parameter.
    const double rlatr = 3.6; // potential growth rate of lateral roots, cm/day.
    const double rtran = 0.2; // the ratio of root mass transferred to a new soil
    // soil cell, when a lateral root grows into it.
    //     On its initiation, lateral root length is assumed to be equal to the width
    //  of a soil column soil cell at the location of the taproot.
    int klocp1 = plant_row_column + 1;
    if (rlat2[l] <= 0)
        rlat2[l] = wk(klocp1, row_space);
    double stday; // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][klocp1] - 273.161;
    double temprg; // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
    // define the column with the tip of this lateral root (ktip)
    int ktip = 0; // column with the tips of the laterals to the right
    double sumwk = 0;
    for (int k = klocp1; k < nk; k++)
    {
        sumwk += wk(k, row_space);
        if (sumwk >= rlat2[l])
        {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the right. Potential
    //  growth rate is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if (state.soil.cells[l][ktip].water_content < PoreSpace[l])
    {
        rlat2[l] += rlatr * temprg * (1 - p1 + state.soil.cells[l][ktip].root.growth_factor * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran) of
        //	mass of roots is transferred to the new soil cell.
        if (rlat2[l] > sumwk && ktip < nk - 1)
        {
            int newktip;        // column into which the tip of the lateral grows to left.
            newktip = ktip + 1; // column into which the tip of the lateral grows to left.
            for (int i = 0; i < NumRootAgeGroups; i++)
            {
                double tran = state.soil.cells[l][ktip].root.weight[i] * rtran;
                state.soil.cells[l][ktip].root.weight[i] -= tran;
                state.soil.cells[l][newktip].root.weight[i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if (state.soil.cells[l][newktip].root.age == 0)
                state.soil.cells[l][newktip].root.age = 0.01;
            if (newktip > state.soil.layers[l].number_of_right_columns_with_root)
                state.soil.layers[l].number_of_right_columns_with_root = newktip;
        }
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
