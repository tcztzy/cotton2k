//  RootGrowth_2.cpp
//
//   functions in this file:
// TapRootGrowth()
// InitiateLateralRoots()
// LateralRootGrowthLeft()
// LateralRootGrowthRight()
// RootDeath()
// RootCultivation()
// RootSummation()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//
//////////////////////////////
void TapRootGrowth()
//     This function computes the elongation of the taproot. It is
//  called from ActualRootGrowth(). It calls SoilTemOnRootGrowth().
//
//     The following global variables are referenced here:
//       dl, nl, NumRootAgeGroups, PlantRowColumn, PoreSpace, RootGroFactor,
//       SoilTempDailyAvrg, VolWaterContent.
//     The following global variables are set here:
//       DepthLastRootLayer, LastTaprootLayer, NumLayersWithRoots, RootAge,
//       RootColNumLeft, RootColNumRight, RootWeight, TapRootLength.
{
    //     The following constant parameters are used:
    const double p1 = 0.10;  // constant parameter.
    const double rtapr = 4;  // potential growth rate of the taproot, cm/day.
    //     It is assumed that taproot elongation takes place irrespective
    //  of the supply of carbon to the roots. This elongation occurs in the
    //  two columns of the slab where the plant is located.
    //     Tap root elongation does not occur in water logged soil (water
    //     table).
    int klocp1 = PlantRowColumn +
                 1;  // the second column in which taproot growth occurs.
    if (VolWaterContent[LastTaprootLayer][PlantRowColumn] >=
            PoreSpace[LastTaprootLayer] ||
        VolWaterContent[LastTaprootLayer][klocp1] >=
            PoreSpace[LastTaprootLayer])
        return;
    //     Average soil resistance (avres) is computed at the root tip.
    // avres = average value of RootGroFactor for the two soil cells at the tip
    // of the taproot.
    double avres = 0.5 * (RootGroFactor[LastTaprootLayer][PlantRowColumn] +
                          RootGroFactor[LastTaprootLayer][klocp1]);
    //     It is assumed that a linear empirical function of avres controls the
    //     rate of
    //  taproot elongation. The potential elongation rate of the taproot is also
    //  modified by soil temperature (SoilTemOnRootGrowth function), soil
    //  resistance, and soil moisture near the root tip.
    //     Actual growth is added to the taproot length TapRootLength.
    double stday;  // daily average soil temperature (C) at root tip.
    stday = 0.5 * (SoilTempDailyAvrg[LastTaprootLayer][PlantRowColumn] +
                   SoilTempDailyAvrg[LastTaprootLayer][klocp1]) -
            273.161;
    double addtaprt;  // added taproot length, cm
    addtaprt = rtapr * (1 - p1 + avres * p1) * SoilTemOnRootGrowth(stday);
    TapRootLength += addtaprt;
    //     DepthLastRootLayer, the depth (in cm) to the end of the last layer
    //     with
    //  roots, is used to check if the taproot reaches a new soil layer.
    //  When the new value of TapRootLength is greater than DepthLastRootLayer -
    //  it means that the roots penetrate to a new soil layer. In this case, and
    //  if this is not the last layer in the slab, the following is executed:
    //     LastTaprootLayer and DepthLastRootLayer are incremented. If this is a
    //     new layer with
    //  roots, NumLayersWithRoots is also redefined and two soil cells of the
    //  new layer are defined as containing roots (by initializing
    //  RootColNumLeft and RootColNumRight).
    if (LastTaprootLayer > nl - 2 || TapRootLength <= DepthLastRootLayer)
        return;
    //     The following is executed when the taproot reaches a new soil layer.
    LastTaprootLayer++;
    DepthLastRootLayer += dl[LastTaprootLayer];
    if (LastTaprootLayer > NumLayersWithRoots - 1) {
        NumLayersWithRoots = LastTaprootLayer + 1;
        if (NumLayersWithRoots > nl) NumLayersWithRoots = nl;
    }
    if (RootColNumLeft[LastTaprootLayer] == 0 ||
        RootColNumLeft[LastTaprootLayer] > PlantRowColumn)
        RootColNumLeft[LastTaprootLayer] = PlantRowColumn;
    if (RootColNumRight[LastTaprootLayer] == 0 ||
        RootColNumRight[LastTaprootLayer] < klocp1)
        RootColNumRight[LastTaprootLayer] = klocp1;
    //     RootAge is initialized for these soil cells.
    RootAge[LastTaprootLayer][PlantRowColumn] = 0.01;
    RootAge[LastTaprootLayer][klocp1] = 0.01;
    //     Some of the mass of class 1 roots is transferred downwards to
    //  the new cells. The transferred mass is proportional to 2 cm of
    //  layer width, but it is not more than half the existing mass in the
    //  last layer.
    for (int i = 0; i < NumRootAgeGroups; i++) {
        double tran;  // root mass transferred to the cell below when the
                      // elongating taproot reaches a new soil layer.
        // first column
        tran = RootWeight[LastTaprootLayer - 1][PlantRowColumn][i] * 2 /
               dl[LastTaprootLayer - 1];
        if (tran > 0.5 * RootWeight[LastTaprootLayer - 1][PlantRowColumn][i])
            tran = 0.5 * RootWeight[LastTaprootLayer - 1][PlantRowColumn][i];
        RootWeight[LastTaprootLayer][PlantRowColumn][i] += tran;
        RootWeight[LastTaprootLayer - 1][PlantRowColumn][i] -= tran;
        // second column
        tran = RootWeight[LastTaprootLayer - 1][klocp1][i] * 2 /
               dl[LastTaprootLayer - 1];
        if (tran > 0.5 * RootWeight[LastTaprootLayer - 1][klocp1][i])
            tran = 0.5 * RootWeight[LastTaprootLayer - 1][klocp1][i];
        RootWeight[LastTaprootLayer][klocp1][i] += tran;
        RootWeight[LastTaprootLayer - 1][klocp1][i] -= tran;
    }
}
//////////////////////////////
void InitiateLateralRoots()
//     This function initiates lateral root growth. It is called from
//     ComputeActualRootGrowth().
//
//     The following global variables are referenced here:
//       DepthLastRootLayer, dl, LastTaprootLayer, TapRootLength
//     The following global variable is set here:      LateralRootFlag
{
    const double distlr =
        12;  // the minimum distance, in cm, from the tip of the taproot, for a
             // lateral root to be able to grow.
    double sdl;  // distance of a layer from tip of taproot, cm.
    sdl = TapRootLength - DepthLastRootLayer;
    //     Loop on soil layers, from the lowest layer with roots upward:
    for (int l = LastTaprootLayer; l >= 0; l--) {
        //     Compute distance from tip of taproot.
        sdl += dl[l];
        //     If a layer is marked for a lateral (LateralRootFlag[l] = 1) and
        //     its
        //  distance from the tip is larger than distlr - initiate a lateral
        //  (LateralRootFlag[l] = 2).
        if (sdl > distlr && LateralRootFlag[l] == 1) LateralRootFlag[l] = 2;
    }
}
//////////////////////////////
void LateralRootGrowthLeft(int l)
//     This function computes the elongation of the lateral roots
//  in a soil layer(l) to the left. It is called from ActualRootGrowth().
//     It calls function SoilTemOnRootGrowth().
//
//     The following global variables are referenced here:
//       NumRootAgeGroups, PlantRowColumn, PoreSpace, RootGroFactor,
//       SoilTempDailyAvrg, VolWaterContent, wk.
//     The following global variables are set here:
//       RootAge, RootColNumLeft, RootWeight.
//     The argument used:     l - layer number in the soil slab.
{
    //     The following constant parameters are used:
    const double p1 = 0.10;  // constant parameter.
    const double rlatr =
        3.6;  // potential growth rate of lateral roots, cm/day.
    const double rtran =
        0.2;  // the ratio of root mass transferred to a new soil
              // soil cell, when a lateral root grows into it.
    //     On its initiation, lateral root length is assumed to be equal to the
    //  width of a soil column soil cell at the location of the taproot.
    if (rlat1[l] <= 0) rlat1[l] = wk[PlantRowColumn];
    double stday;  // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][PlantRowColumn] - 273.161;
    double temprg;  // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
    //     Define the column with the tip of this lateral root (ktip)
    int ktip = 0;      // column with the tips of the laterals to the left
    double sumwk = 0;  // summation of columns width
    for (int k = PlantRowColumn; k >= 0; k--) {
        sumwk += wk[k];
        if (sumwk >= rlat1[l]) {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the left. Potential
    //  growth rate (u) is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if (VolWaterContent[l][ktip] < PoreSpace[l]) {
        rlat1[l] += rlatr * temprg * (1 - p1 + RootGroFactor[l][ktip] * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran)
        //     of
        //	mass of roots is transferred to the new soil cell.
        if (rlat1[l] > sumwk && ktip > 0) {
            int newktip;  // column into which the tip of the lateral grows to
                          // left.
            newktip = ktip - 1;
            for (int i = 0; i < NumRootAgeGroups; i++) {
                double tran = RootWeight[l][ktip][i] * rtran;
                RootWeight[l][ktip][i] -= tran;
                RootWeight[l][newktip][i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if (RootAge[l][newktip] == 0) RootAge[l][newktip] = 0.01;
            if (newktip < RootColNumLeft[l]) RootColNumLeft[l] = newktip;
        }
    }
}
//////////////////////////////
void LateralRootGrowthRight(int l)
//     This function computes the elongation of the lateral roots
//  in a soil layer(l) to the right. It is called from ActualRootGrowth().
//     It calls function SoilTemOnRootGrowth().
//
//     The following global variables are referenced here:
//  nk, NumRootAgeGroups, PlantRowColumn, PoreSpace, RootGroFactor,
//  SoilTempDailyAvrg, VolWaterContent, wk.
//     The following global variables are set here:
//  RootAge, RootColNumRight, RootWeight.
//     The argument used:      l - layer number in the soil slab.
{
    //     The following constant parameters are used:
    const double p1 = 0.10;  // constant parameter.
    const double rlatr =
        3.6;  // potential growth rate of lateral roots, cm/day.
    const double rtran =
        0.2;  // the ratio of root mass transferred to a new soil
              // soil cell, when a lateral root grows into it.
    //     On its initiation, lateral root length is assumed to be equal to the
    //     width
    //  of a soil column soil cell at the location of the taproot.
    int klocp1 = PlantRowColumn + 1;
    if (rlat2[l] <= 0) rlat2[l] = wk[klocp1];
    double stday;  // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][klocp1] - 273.161;
    double temprg;  // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
    // define the column with the tip of this lateral root (ktip)
    int ktip = 0;  // column with the tips of the laterals to the right
    double sumwk = 0;
    for (int k = klocp1; k < nk; k++) {
        sumwk += wk[k];
        if (sumwk >= rlat2[l]) {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the right. Potential
    //  growth rate is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if (VolWaterContent[l][ktip] < PoreSpace[l]) {
        rlat2[l] += rlatr * temprg * (1 - p1 + RootGroFactor[l][ktip] * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran)
        //     of
        //	mass of roots is transferred to the new soil cell.
        if (rlat2[l] > sumwk && ktip < nk - 1) {
            int newktip;  // column into which the tip of the lateral grows to
                          // left.
            newktip =
                ktip +
                1;  // column into which the tip of the lateral grows to left.
            for (int i = 0; i < NumRootAgeGroups; i++) {
                double tran = RootWeight[l][ktip][i] * rtran;
                RootWeight[l][ktip][i] -= tran;
                RootWeight[l][newktip][i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if (RootAge[l][newktip] == 0) RootAge[l][newktip] = 0.01;
            if (newktip > RootColNumRight[l]) RootColNumRight[l] = newktip;
        }
    }
}
//////////////////////////////
void RootDeath(int l, int k) {
    //     This function computes the death of root tissue in each soil cell
    //     containing roots.
    //  When root age reaches a threshold thdth(i), a proportion dth(i) of the
    //  roots in class i dies. The mass of dead roots is added to DailyRootLoss.
    //     It has been adapted from GOSSYM, but the threshold age for this
    //     process is based on
    //  the time from when the roots first grew into each soil cell.
    //     It is assumed that root death rate is greater in dry soil, for all
    //     root classes except
    //  class 1. Root death rate is increased to the maximum value in soil
    //  saturated with water.
    //
    //     The following global variables are referenced here:
    //       RootAge, PoreSpace, SoilPsi, VolWaterContent
    //     The following global variables are set here:
    //       RootWeight, DailyRootLoss
    //     The arguments k, l - are column and layer numbers.
    //
    //     The constant parameters are used:
    const double aa =
        0.008;  // a parameter in the equation for computing dthfac.
    const double dth[3] = {
        0.0001, 0.0002,
        0.0001};  // the daily proportion of death of root tissue.
    const double dthmax =
        0.10;  // a parameter in the equation for computing dthfac.
    const double psi0 =
        -14.5;  // a parameter in the equation for computing dthfac.
    const double thdth[3] = {
        30.0, 50.0, 100.0};  // the time threshold, from the initial
                             // penetration of roots to a soil cell, after
                             // which death of root tissue of class i may occur.
                             //
    for (int i = 0; i < 3; i++) {
        if (RootAge[l][k] > thdth[i]) {
            double dthfac;  // the computed proportion of roots dying in each
                            // class.
            dthfac = dth[i];
            if (VolWaterContent[l][k] >= PoreSpace[l])
                dthfac = dthmax;
            else {
                if (i <= 1 && SoilPsi[l][k] <= psi0)
                    dthfac += aa * (psi0 - SoilPsi[l][k]);
                if (dthfac > dthmax) dthfac = dthmax;
            }
            DailyRootLoss += RootWeight[l][k][i] * dthfac;
            RootWeight[l][k][i] -= RootWeight[l][k][i] * dthfac;
        }
    }
}
//////////////////////////////
void RootCultivation(int j)
//     This function is executed on the day of soil cultivation. It is called
//     from
//  ActualRootGrowth(). It has been adapted from GOSSYM. It is assumed that the
//  roots in the upper soil layers, as defined by the depth of cultivation, are
//  destroyed, with the exception of the soil cells that are within 15 cm of the
//  plant row.
//
//     The following global variables are referenced here:
//       dl, CultivationDepth, NumRootAgeGroups, nk, nl, PlantRowLocation, wk
//     The following global variables are set here:
//       RootWeight, DailyRootLoss
//     The argument j - is the serial number of this cultivation.
{
    //     The depth of cultivation (CultivationDepth) is in cm. The number of
    //     layers affected by
    //  it (lcult) is determined. Loop for all columns, and check if this column
    //  is more than 15 cm from the plant row: if this is true, destroy all
    //  roots and add their weight to DailyRootLoss.
    int lcult = 0;     // number of soil layers affected by cultivation.
    double sdpth = 0;  // sum depth to the end of the layer.
    for (int l = 0; l < nl; l++) {
        sdpth += dl[l];
        if (sdpth >= CultivationDepth[j]) {
            lcult = l;
            break;
        }
    }
    //
    double sumwk = 0;  // sum of column widths from edge of slab to this column.
    for (int k = 0; k < nk; k++) {
        sumwk += wk[k];
        if (fabs(sumwk - PlantRowLocation) >= 15) {
            for (int l = 0; l < lcult; l++)
                for (int i = 0; i < NumRootAgeGroups; i++) {
                    DailyRootLoss += RootWeight[l][k][i];
                    RootWeight[l][k][i] = 0;
                }
        }
    }
}
//////////////////////////////
void RootSummation()
//     This function has been added for compatibility with GOSSYM root routines.
//  It is called from ActualRootGrowth(). It summarizes root data, in a form
//  ready for output or plotting. Sums of root weights for cells, for age groups
//  and for the total slab are calculated. TotalRootWeight is calculated in g
//  per plant.
//
//     The following global variables are referenced here:
//  dl, Kday, nk, nl, NumLayersWithRoots, NumRootAgeGroups, PerPlantArea,
//  RootWeight, RootWtCapblUptake, RowSpace, wk
//     The following global variable is set here:     TotalRootWeight
{
    //     Compute the root weight (of all age classes) for all soil cells as
    //  rootsv, and the total as roots.
    double roots = 0;  // total weight of roots of all classes, g per slab.
    double rootsv[maxl][maxk];  // total dry weight of roots in a cell, g per
                                // soil cell.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++) {
            rootsv[l][k] = 0;
            for (int i = 0; i < 3; i++) rootsv[l][k] += RootWeight[l][k][i];
            roots += rootsv[l][k];
        }
    //
    for (int l = 0; l < maxl; l++)
        for (int k = 0; k < maxk; k++)
            Scratch21[DayOfSimulation - 1].rootsv[l][k] = rootsv[l][k];
    //     Convert total root weight from g per slab to g per plant.
    TotalRootWeight = roots * 100 * PerPlantArea / RowSpace;
}
