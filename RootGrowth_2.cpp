//  RootGrowth_2.cpp
//
//   functions in this file:
// RedistRootNewGrowth()
// TapRootGrowth()
// InitiateLateralRoots()
// LateralRootGrowthLeft()
// LateralRootGrowthRight()
// RootAging()
// RootDeath()
// RootCultivation()
// RootSummation()
//
#include <numeric>
#include <fstream>
#include "global.h"
#include "Simulation.h"

using namespace std;

extern "C"
{
    double SoilTemOnRootGrowth(double);
}

//////////////////////////////////////////////////
void RedistRootNewGrowth(Simulation &sim, uint32_t u, int l, int k, double addwt)
//     This function computes the redistribution of new growth of
//  roots into adjacent soil cells. It is called from ActualRootGrowth().
//     Redistribution is affected by the factors rgfdn, rgfsd, rgfup.
//  and the values of RootGroFactor(l,k) in this soil cell and in the adjacent
//  cells. The values of ActualRootGrowth(l,k) for this and for the adjacent
//  soil cells are computed. The code of this module is based, with
//  major changes, on the code of GOSSYM.
//
//     The following arguments are referenced:
//       addwt - actual growth rate of roots in this soil cell.
//       k, l - column and layer numbers.
//     The following global variables are referenced here:
//       dl, nk, nl, PlantRowColumn, RootGroFactor, wk.
//     The following global variables are set here:
//       ActualRootGrowth, DepthLastRootLayer, LastTaprootLayer,
//       RootAge, RootColNumLeft, RootColNumRight, TapRootLength.
{
    State &state = sim.states[u];
    //     The following constant parameters are used. These are relative factors for root growth
    // to adjoining cells, downwards, sideways, and upwards, respectively. These factors are
    // relative to the volume of the soil cell from which growth originates.
    const double rgfdn = 900;
    const double rgfsd = 600;
    const double rgfup = 10;
    //     Set the number of layer above and below this layer, and
    //  the number of columns to the right and to the left of this column.
    int lm1, lp1; // layer above and below layer l.
    if (l == nl - 1)
        lp1 = l;
    else
        lp1 = l + 1;
    if (l == 0)
        lm1 = l;
    else
        lm1 = l - 1;
    //
    int km1, kp1; // column to the left and to the right of column k.
    if (k == nk - 1)
        kp1 = k;
    else
        kp1 = k + 1;
    if (k == 0)
        km1 = k;
    else
        km1 = k - 1;
    //     Compute proportionality factors (efac1, efacl, efacr, efacu, efacd) as the product
    //  of RootGroFactor and the geotropic factors in the respective soil cells. Note that
    //  the geotropic factors are relative to the volume of the soil cell.
    //     Compute the sum srwp of the proportionality factors.
    double efac1; // product of RootGroFactor and geotropic factor for this cell.
    double efacd; // as efac1 for the cell below this cell.
    double efacl; // as efac1 for the cell to the left of this cell.
    double efacr; // as efac1 for the cell to the right of this cell.
    double efacu; // as efac1 for the cell above this cell.
    double srwp;  // sum of all efac values.
    efac1 = dl[l] * wk[k] * state.root[l][k].growth_factor;
    efacl = rgfsd * state.root[l][km1].growth_factor;
    efacr = rgfsd * state.root[l][kp1].growth_factor;
    efacu = rgfup * state.root[lm1][k].growth_factor;
    efacd = rgfdn * state.root[lp1][k].growth_factor;
    srwp = efac1 + efacl + efacr + efacu + efacd;
    //     If srwp is very small, all the added weight will be in the
    //  same soil soil cell, and execution of this function is ended.
    if (srwp < 1e-10)
    {
        sim.states[u].root[l][k].actual_growth = addwt;
        return;
    }
    //     Allocate the added dry matter to this and the adjoining
    //  soil cells in proportion to the EFAC factors.
    sim.states[u].root[l][k].actual_growth += addwt * efac1 / srwp;
    sim.states[u].root[l][km1].actual_growth += addwt * efacl / srwp;
    sim.states[u].root[l][kp1].actual_growth += addwt * efacr / srwp;
    sim.states[u].root[lm1][k].actual_growth += addwt * efacu / srwp;
    sim.states[u].root[lp1][k].actual_growth += addwt * efacd / srwp;
    //     If roots are growing into new soil soil cells, initialize
    //  their RootAge to 0.01.
    if (sim.states[u].root[l][km1].age == 0)
        sim.states[u].root[l][km1].age = 0.01;
    if (sim.states[u].root[l][kp1].age == 0)
        sim.states[u].root[l][kp1].age = 0.01;
    if (sim.states[u].root[lm1][k].age == 0)
        sim.states[u].root[lm1][k].age = 0.01;
    //     If this new compartmment is in a new layer with roots, also
    //  initialize its RootColNumLeft and RootColNumRight values.
    if (sim.states[u].root[lp1][k].age == 0 && efacd > 0)
    {
        sim.states[u].root[lp1][k].age = 0.01;
        if (RootColNumLeft[lp1] == 0 || k < RootColNumLeft[lp1])
            RootColNumLeft[lp1] = k;
        if (RootColNumRight[lp1] == 0 || k > RootColNumRight[lp1])
            RootColNumRight[lp1] = k;
    }
    //     If this is in the location of the taproot, and the roots reach a new soil layer,
    //  update the taproot parameters TapRootLength, DepthLastRootLayer, and LastTaprootLayer.
    if (k == sim.plant_row_column || k == sim.plant_row_column + 1)
        if (lp1 > LastTaprootLayer && efacd > 0)
        {
            TapRootLength = DepthLastRootLayer + 0.01;
            DepthLastRootLayer += dl[lp1];
            LastTaprootLayer = lp1;
        }
    //     Update NumLayersWithRoots, if necessary, and the values of RootColNumLeft and RootColNumRight for
    //  this layer.
    if (sim.states[u].number_of_layers_with_root <= l && efacd > 0)
        sim.states[u].number_of_layers_with_root = l + 1;
    if (km1 < RootColNumLeft[l])
        RootColNumLeft[l] = km1;
    if (kp1 > RootColNumRight[l])
        RootColNumRight[l] = kp1;
}

//////////////////////////////
void TapRootGrowth(Simulation &sim, uint32_t u, const int &NumRootAgeGroups)
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
    State &state = sim.states[u];
    //     The following constant parameters are used:
    const double p1 = 0.10;                // constant parameter.
    const double rtapr = 4;                // potential growth rate of the taproot, cm/day.
                                           //     It is assumed that taproot elongation takes place irrespective
                                           //  of the supply of carbon to the roots. This elongation occurs in the
                                           //  two columns of the slab where the plant is located.
                                           //     Tap root elongation does not occur in water logged soil (water table).
    int klocp1 = sim.plant_row_column + 1; // the second column in which taproot growth occurs.
    if (VolWaterContent[LastTaprootLayer][sim.plant_row_column] >= PoreSpace[LastTaprootLayer] || VolWaterContent[LastTaprootLayer][klocp1] >= PoreSpace[LastTaprootLayer])
        return;
    //     Average soil resistance (avres) is computed at the root tip.
    // avres = average value of RootGroFactor for the two soil cells at the tip of the taproot.
    double avres = 0.5 * (state.root[LastTaprootLayer][sim.plant_row_column].growth_factor + state.root[LastTaprootLayer][klocp1].growth_factor);
    //     It is assumed that a linear empirical function of avres controls the rate of
    //  taproot elongation. The potential elongation rate of the taproot is also modified by
    //  soil temperature (SoilTemOnRootGrowth function), soil resistance, and soil
    //  moisture near the root tip.
    //     Actual growth is added to the taproot length TapRootLength.
    double stday; // daily average soil temperature (C) at root tip.
    stday = 0.5 * (SoilTempDailyAvrg[LastTaprootLayer][sim.plant_row_column] + SoilTempDailyAvrg[LastTaprootLayer][klocp1]) - 273.161;
    double addtaprt; // added taproot length, cm
    addtaprt = rtapr * (1 - p1 + avres * p1) * SoilTemOnRootGrowth(stday);
    TapRootLength += addtaprt;
    //     DepthLastRootLayer, the depth (in cm) to the end of the last layer with
    //  roots, is used to check if the taproot reaches a new soil layer.
    //  When the new value of TapRootLength is greater than DepthLastRootLayer - it means that the
    //  roots penetrate to a new soil layer. In this case, and if this
    //  is not the last layer in the slab, the following is executed:
    //     LastTaprootLayer and DepthLastRootLayer are incremented. If this is a new layer with
    //  roots, NumLayersWithRoots is also redefined and two soil cells of the new layer
    //  are defined as containing roots (by initializing RootColNumLeft and RootColNumRight).
    if (LastTaprootLayer > nl - 2 || TapRootLength <= DepthLastRootLayer)
        return;
    //     The following is executed when the taproot reaches a new soil layer.
    LastTaprootLayer++;
    DepthLastRootLayer += dl[LastTaprootLayer];
    if (LastTaprootLayer > sim.states[u].number_of_layers_with_root - 1)
    {
        sim.states[u].number_of_layers_with_root = LastTaprootLayer + 1;
        if (sim.states[u].number_of_layers_with_root > nl)
            sim.states[u].number_of_layers_with_root = nl;
    }
    if (RootColNumLeft[LastTaprootLayer] == 0 ||
        RootColNumLeft[LastTaprootLayer] > sim.plant_row_column)
        RootColNumLeft[LastTaprootLayer] = sim.plant_row_column;
    if (RootColNumRight[LastTaprootLayer] == 0 ||
        RootColNumRight[LastTaprootLayer] < klocp1)
        RootColNumRight[LastTaprootLayer] = klocp1;
    //     RootAge is initialized for these soil cells.
    sim.states[u].root[LastTaprootLayer][sim.plant_row_column].age = 0.01;
    sim.states[u].root[LastTaprootLayer][klocp1].age = 0.01;
    //     Some of the mass of class 1 roots is transferred downwards to
    //  the new cells. The transferred mass is proportional to 2 cm of
    //  layer width, but it is not more than half the existing mass in the
    //  last layer.
    for (int i = 0; i < NumRootAgeGroups; i++)
    {
        double tran; // root mass transferred to the cell below when the elongating taproot
        // reaches a new soil layer.
        // first column
        tran = sim.states[u].root[LastTaprootLayer - 1][sim.plant_row_column].weight[i] * 2 / dl[LastTaprootLayer - 1];
        if (tran > 0.5 * sim.states[u].root[LastTaprootLayer - 1][sim.plant_row_column].weight[i])
            tran = 0.5 * sim.states[u].root[LastTaprootLayer - 1][sim.plant_row_column].weight[i];
        sim.states[u].root[LastTaprootLayer][sim.plant_row_column].weight[i] += tran;
        sim.states[u].root[LastTaprootLayer - 1][sim.plant_row_column].weight[i] -= tran;
        // second column
        tran = sim.states[u].root[LastTaprootLayer - 1][klocp1].weight[i] * 2 / dl[LastTaprootLayer - 1];
        if (tran > 0.5 * sim.states[u].root[LastTaprootLayer - 1][klocp1].weight[i])
            tran = 0.5 * sim.states[u].root[LastTaprootLayer - 1][klocp1].weight[i];
        sim.states[u].root[LastTaprootLayer][klocp1].weight[i] += tran;
        sim.states[u].root[LastTaprootLayer - 1][klocp1].weight[i] -= tran;
    }
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
        sdl += dl[l];
        //     If a layer is marked for a lateral (LateralRootFlag[l] = 1) and its
        //  distance from the tip is larger than distlr - initiate a lateral
        //  (LateralRootFlag[l] = 2).
        if (sdl > distlr && LateralRootFlag[l] == 1)
            LateralRootFlag[l] = 2;
    }
}

//////////////////////////////
void LateralRootGrowthLeft(Simulation &sim, uint32_t u, int l, const int &NumRootAgeGroups)
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
    State &state = sim.states[u];
    //     The following constant parameters are used:
    const double p1 = 0.10;   // constant parameter.
    const double rlatr = 3.6; // potential growth rate of lateral roots, cm/day.
    const double rtran = 0.2; // the ratio of root mass transferred to a new soil
    // soil cell, when a lateral root grows into it.
    //     On its initiation, lateral root length is assumed to be equal to the
    //  width of a soil column soil cell at the location of the taproot.
    if (rlat1[l] <= 0)
        rlat1[l] = wk[sim.plant_row_column];
    double stday; // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][sim.plant_row_column] - 273.161;
    double temprg; // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
    //     Define the column with the tip of this lateral root (ktip)
    int ktip = 0;     // column with the tips of the laterals to the left
    double sumwk = 0; // summation of columns width
    for (int k = sim.plant_row_column; k >= 0; k--)
    {
        sumwk += wk[k];
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
    if (VolWaterContent[l][ktip] < PoreSpace[l])
    {
        rlat1[l] += rlatr * temprg * (1 - p1 + state.root[l][ktip].growth_factor * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran) of
        //	mass of roots is transferred to the new soil cell.
        if (rlat1[l] > sumwk && ktip > 0)
        {
            int newktip; // column into which the tip of the lateral grows to left.
            newktip = ktip - 1;
            for (int i = 0; i < NumRootAgeGroups; i++)
            {
                double tran = sim.states[u].root[l][ktip].weight[i] * rtran;
                sim.states[u].root[l][ktip].weight[i] -= tran;
                sim.states[u].root[l][newktip].weight[i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if (sim.states[u].root[l][newktip].age == 0)
                sim.states[u].root[l][newktip].age = 0.01;
            if (newktip < RootColNumLeft[l])
                RootColNumLeft[l] = newktip;
        }
    }
}

//////////////////////////////
void LateralRootGrowthRight(Simulation &sim, uint32_t u, int l, const int &NumRootAgeGroups)
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
    State &state = sim.states[u];
    //     The following constant parameters are used:
    const double p1 = 0.10;   // constant parameter.
    const double rlatr = 3.6; // potential growth rate of lateral roots, cm/day.
    const double rtran = 0.2; // the ratio of root mass transferred to a new soil
    // soil cell, when a lateral root grows into it.
    //     On its initiation, lateral root length is assumed to be equal to the width
    //  of a soil column soil cell at the location of the taproot.
    int klocp1 = sim.plant_row_column + 1;
    if (rlat2[l] <= 0)
        rlat2[l] = wk[klocp1];
    double stday; // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][klocp1] - 273.161;
    double temprg; // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
    // define the column with the tip of this lateral root (ktip)
    int ktip = 0; // column with the tips of the laterals to the right
    double sumwk = 0;
    for (int k = klocp1; k < nk; k++)
    {
        sumwk += wk[k];
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
    if (VolWaterContent[l][ktip] < PoreSpace[l])
    {
        rlat2[l] += rlatr * temprg * (1 - p1 + state.root[l][ktip].growth_factor * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran) of
        //	mass of roots is transferred to the new soil cell.
        if (rlat2[l] > sumwk && ktip < nk - 1)
        {
            int newktip;        // column into which the tip of the lateral grows to left.
            newktip = ktip + 1; // column into which the tip of the lateral grows to left.
            for (int i = 0; i < NumRootAgeGroups; i++)
            {
                double tran = sim.states[u].root[l][ktip].weight[i] * rtran;
                sim.states[u].root[l][ktip].weight[i] -= tran;
                sim.states[u].root[l][newktip].weight[i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if (sim.states[u].root[l][newktip].age == 0)
                sim.states[u].root[l][newktip].age = 0.01;
            if (newktip > RootColNumRight[l])
                RootColNumRight[l] = newktip;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void RootAging(Root &root, int l, int k)
//     This function is called from ActualRootGrowth(). It updates the variable celage(l,k)
//  for the age of roots in each soil cell containing roots. When root age reaches a threshold
//  thtrn(i), a transformation of root tissue from class i to class i+1 occurs. The proportion
//  transformed is trn(i).
//     It has been adapted from the code of GOSSYM, but the threshold
//  age for this process is based on the time from when the roots first
//  grew into each soil cell (whereas the time from emergence was used
//  in GOSSYM). Note: only 3 root age groups are assumed here.
//
//     The following global variable is referenced here:       SoilTempDailyAvrg.
//     The following global variables are set here:        RootAge, RootWeight.
//     The arguments k, l - are column and layer numbers.
{
    //     The following constant parameters are used:
    const double thtrn[2] = {20.0, 40.0};   // the time threshold, from the initial
                                            // penetration of roots to a soil cell, after which some of the root
                                            // mass of class i may be transferred into the next class (i+1).
    const double trn[2] = {0.0060, 0.0050}; //  the daily proportion of this transfer.
                                            //
    double stday;                           // daily average soil temperature (c) of soil cell.
    stday = SoilTempDailyAvrg[l][k] - 273.161;
    root.age += SoilTemOnRootGrowth(stday);
    //
    for (int i = 0; i < 2; i++)
    {
        if (root.age > thtrn[i])
        {
            double xtr; // root mass transferred from one class to the next.
            xtr = trn[i] * root.weight[i];
            root.weight[i + 1] += xtr;
            root.weight[i] -= xtr;
        }
    }
}

//////////////////////////////
double RootDeath(Root &root, int l, int k, double DailyRootLoss)
{
    //     This function computes the death of root tissue in each soil cell containing roots.
    //  When root age reaches a threshold thdth(i), a proportion dth(i) of the roots in class i
    //  dies. The mass of dead roots is added to DailyRootLoss.
    //     It has been adapted from GOSSYM, but the threshold age for this process is based on
    //  the time from when the roots first grew into each soil cell.
    //     It is assumed that root death rate is greater in dry soil, for all root classes except
    //  class 1. Root death rate is increased to the maximum value in soil saturated with water.
    //
    //     The following global variables are referenced here:
    //       RootAge, PoreSpace, SoilPsi, VolWaterContent
    //     The following global variables are set here:
    //       RootWeight, DailyRootLoss
    //     The arguments k, l - are column and layer numbers.
    //
    //     The constant parameters are used:
    const double aa = 0.008;                        // a parameter in the equation for computing dthfac.
    const double dth[3] = {0.0001, 0.0002, 0.0001}; // the daily proportion of death of root tissue.
    const double dthmax = 0.10;                     // a parameter in the equation for computing dthfac.
    const double psi0 = -14.5;                      // a parameter in the equation for computing dthfac.
    const double thdth[3] = {30.0, 50.0, 100.0};    // the time threshold, from the initial
    // penetration of roots to a soil cell, after
    // which death of root tissue of class i may occur.
    //
    for (int i = 0; i < 3; i++)
    {
        if (root.age > thdth[i])
        {
            double dthfac; // the computed proportion of roots dying in each class.
            dthfac = dth[i];
            if (VolWaterContent[l][k] >= PoreSpace[l])
                dthfac = dthmax;
            else
            {
                if (i <= 1 && SoilPsi[l][k] <= psi0)
                    dthfac += aa * (psi0 - SoilPsi[l][k]);
                if (dthfac > dthmax)
                    dthfac = dthmax;
            }
            DailyRootLoss += root.weight[i] * dthfac;
            root.weight[i] -= root.weight[i] * dthfac;
        }
    }
    return DailyRootLoss;
}

//////////////////////////////
double RootCultivation(Root root[40][20], int j, const int &NumRootAgeGroups, double DailyRootLoss)
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
        sdpth += dl[l];
        if (sdpth >= CultivationDepth[j])
        {
            lcult = l;
            break;
        }
    }
    //
    double sumwk = 0; // sum of column widths from edge of slab to this column.
    for (int k = 0; k < nk; k++)
    {
        sumwk += wk[k];
        if (fabs(sumwk - PlantRowLocation) >= 15)
        {
            for (int l = 0; l < lcult; l++)
                for (int i = 0; i < NumRootAgeGroups; i++)
                {
                    DailyRootLoss += root[l][k].weight[i];
                    root[l][k].weight[i] = 0;
                }
        }
    }
    return DailyRootLoss;
}

//////////////////////////////
void RootSummation(Simulation &sim, uint32_t u, const int &NumRootAgeGroups)
//     This function has been added for compatibility with GOSSYM root routines.
//  It is called from ActualRootGrowth(). It summarizes root data, in a form ready
//  for output or plotting. Sums of root weights for cells, for age groups and for
//  the total slab are calculated. TotalRootWeight is calculated in g per plant.
//
//     The following global variables are referenced here:
//  dl, Kday, nk, nl, NumRootAgeGroups, PerPlantArea,
//  OutIndex, RootWeight, RootWtCapblUptake, RowSpace, wk
//     The following global variable is set here:     TotalRootWeight
{
    State &state = sim.states[u];
    //     Compute the total root weight (of all age classes) for all soil cells as
    double roots = 0;          // total weight of roots of all classes, g per slab.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++)
            roots += accumulate(sim.states[u].root[l][k].weight, sim.states[u].root[l][k].weight+3, double(0));
    //     Convert total root weight from g per slab to g per plant.
    TotalRootWeight = roots * 100 * PerPlantArea / sim.row_space;
}
