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
#include <fstream>
#include "global.h"

using namespace std;

extern "C" {
    double SoilTemOnRootGrowth(double);
}

//////////////////////////////////////////////////
tuple<int> RedistRootNewGrowth(int l, int k, double addwt, int NumLayersWithRoots, double RootAge[40][20])
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
    efac1 = dl[l] * wk[k] * RootGroFactor[l][k];
    efacl = rgfsd * RootGroFactor[l][km1];
    efacr = rgfsd * RootGroFactor[l][kp1];
    efacu = rgfup * RootGroFactor[lm1][k];
    efacd = rgfdn * RootGroFactor[lp1][k];
    srwp = efac1 + efacl + efacr + efacu + efacd;
//     If srwp is very small, all the added weight will be in the
//  same soil soil cell, and execution of this function is ended.
    if (srwp < 1e-10) {
        ActualRootGrowth[l][k] = addwt;
        return make_tuple(NumLayersWithRoots);
    }
//     Allocate the added dry matter to this and the adjoining
//  soil cells in proportion to the EFAC factors.
    ActualRootGrowth[l][k] += addwt * efac1 / srwp;
    ActualRootGrowth[l][km1] += addwt * efacl / srwp;
    ActualRootGrowth[l][kp1] += addwt * efacr / srwp;
    ActualRootGrowth[lm1][k] += addwt * efacu / srwp;
    ActualRootGrowth[lp1][k] += addwt * efacd / srwp;
//     If roots are growing into new soil soil cells, initialize
//  their RootAge to 0.01.
    if (RootAge[l][km1] == 0)
        RootAge[l][km1] = 0.01;
    if (RootAge[l][kp1] == 0)
        RootAge[l][kp1] = 0.01;
    if (RootAge[lm1][k] == 0)
        RootAge[lm1][k] = 0.01;
//     If this new compartmment is in a new layer with roots, also
//  initialize its RootColNumLeft and RootColNumRight values.
    if (RootAge[lp1][k] == 0 && efacd > 0) {
        RootAge[lp1][k] = 0.01;
        if (RootColNumLeft[lp1] == 0 || k < RootColNumLeft[lp1])
            RootColNumLeft[lp1] = k;
        if (RootColNumRight[lp1] == 0 || k > RootColNumRight[lp1])
            RootColNumRight[lp1] = k;
    }
//     If this is in the location of the taproot, and the roots reach a new soil layer, 
//  update the taproot parameters TapRootLength, DepthLastRootLayer, and LastTaprootLayer.
    if (k == PlantRowColumn || k == PlantRowColumn + 1)
        if (lp1 > LastTaprootLayer && efacd > 0) {
            TapRootLength = DepthLastRootLayer + 0.01;
            DepthLastRootLayer += dl[lp1];
            LastTaprootLayer = lp1;
        }
//     Update NumLayersWithRoots, if necessary, and the values of RootColNumLeft and RootColNumRight for
//  this layer.
    if (NumLayersWithRoots <= l && efacd > 0)
        NumLayersWithRoots = l + 1;
    if (km1 < RootColNumLeft[l])
        RootColNumLeft[l] = km1;
    if (kp1 > RootColNumRight[l])
        RootColNumRight[l] = kp1;
    return make_tuple(NumLayersWithRoots);
}

//////////////////////////////
tuple<int>
TapRootGrowth(const int &NumRootAgeGroups, int NumLayersWithRoots, double RootWeight[40][20][3], double RootAge[40][20])
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
    const double p1 = 0.10; // constant parameter.
    const double rtapr = 4; // potential growth rate of the taproot, cm/day.
//     It is assumed that taproot elongation takes place irrespective
//  of the supply of carbon to the roots. This elongation occurs in the
//  two columns of the slab where the plant is located.
//     Tap root elongation does not occur in water logged soil (water table).
    int klocp1 = PlantRowColumn + 1; // the second column in which taproot growth occurs.
    if (VolWaterContent[LastTaprootLayer][PlantRowColumn] >= PoreSpace[LastTaprootLayer]
        || VolWaterContent[LastTaprootLayer][klocp1] >= PoreSpace[LastTaprootLayer])
        return make_tuple(NumLayersWithRoots);
//     Average soil resistance (avres) is computed at the root tip.
// avres = average value of RootGroFactor for the two soil cells at the tip of the taproot.
    double avres = 0.5 * (RootGroFactor[LastTaprootLayer][PlantRowColumn]
                          + RootGroFactor[LastTaprootLayer][klocp1]);
//     It is assumed that a linear empirical function of avres controls the rate of 
//  taproot elongation. The potential elongation rate of the taproot is also modified by
//  soil temperature (SoilTemOnRootGrowth function), soil resistance, and soil
//  moisture near the root tip.
//     Actual growth is added to the taproot length TapRootLength.
    double stday; // daily average soil temperature (C) at root tip.
    stday = 0.5 * (SoilTempDailyAvrg[LastTaprootLayer][PlantRowColumn]
                   + SoilTempDailyAvrg[LastTaprootLayer][klocp1]) - 273.161;
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
        return make_tuple(NumLayersWithRoots);
//     The following is executed when the taproot reaches a new soil layer.
    LastTaprootLayer++;
    DepthLastRootLayer += dl[LastTaprootLayer];
    if (LastTaprootLayer > NumLayersWithRoots - 1) {
        NumLayersWithRoots = LastTaprootLayer + 1;
        if (NumLayersWithRoots > nl)
            NumLayersWithRoots = nl;
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
        double tran; // root mass transferred to the cell below when the elongating taproot
        // reaches a new soil layer.
// first column
        tran = RootWeight[LastTaprootLayer - 1][PlantRowColumn][i] * 2 / dl[LastTaprootLayer - 1];
        if (tran > 0.5 * RootWeight[LastTaprootLayer - 1][PlantRowColumn][i])
            tran = 0.5 * RootWeight[LastTaprootLayer - 1][PlantRowColumn][i];
        RootWeight[LastTaprootLayer][PlantRowColumn][i] += tran;
        RootWeight[LastTaprootLayer - 1][PlantRowColumn][i] -= tran;
// second column
        tran = RootWeight[LastTaprootLayer - 1][klocp1][i] * 2 / dl[LastTaprootLayer - 1];
        if (tran > 0.5 * RootWeight[LastTaprootLayer - 1][klocp1][i])
            tran = 0.5 * RootWeight[LastTaprootLayer - 1][klocp1][i];
        RootWeight[LastTaprootLayer][klocp1][i] += tran;
        RootWeight[LastTaprootLayer - 1][klocp1][i] -= tran;
    }
    return make_tuple(NumLayersWithRoots);
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
    double sdl; // distance of a layer from tip of taproot, cm.
    sdl = TapRootLength - DepthLastRootLayer;
//     Loop on soil layers, from the lowest layer with roots upward:
    for (int l = LastTaprootLayer; l >= 0; l--) {
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
void LateralRootGrowthLeft(int l, const int &NumRootAgeGroups, double RootWeight[40][20][3], double RootAge[40][20])
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
    const double p1 = 0.10; // constant parameter.
    const double rlatr = 3.6; // potential growth rate of lateral roots, cm/day.
    const double rtran = 0.2; // the ratio of root mass transferred to a new soil
    // soil cell, when a lateral root grows into it.
//     On its initiation, lateral root length is assumed to be equal to the 
//  width of a soil column soil cell at the location of the taproot.
    if (rlat1[l] <= 0)
        rlat1[l] = wk[PlantRowColumn];
    double stday; // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][PlantRowColumn] - 273.161;
    double temprg; // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
//     Define the column with the tip of this lateral root (ktip)
    int ktip = 0; // column with the tips of the laterals to the left
    double sumwk = 0; // summation of columns width
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
//     If the lateral reaches a new soil soil cell: a proportion (tran) of 
//	mass of roots is transferred to the new soil cell.
        if (rlat1[l] > sumwk && ktip > 0) {
            int newktip; // column into which the tip of the lateral grows to left.
            newktip = ktip - 1;
            for (int i = 0; i < NumRootAgeGroups; i++) {
                double tran = RootWeight[l][ktip][i] * rtran;
                RootWeight[l][ktip][i] -= tran;
                RootWeight[l][newktip][i] += tran;
            }
//     RootAge is initialized for this soil cell.
//     RootColNumLeft of this layer is redefined.
            if (RootAge[l][newktip] == 0)
                RootAge[l][newktip] = 0.01;
            if (newktip < RootColNumLeft[l])
                RootColNumLeft[l] = newktip;
        }
    }
}

//////////////////////////////
void LateralRootGrowthRight(int l, const int &NumRootAgeGroups, double RootWeight[40][20][3], double RootAge[40][20])
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
    const double p1 = 0.10; // constant parameter.
    const double rlatr = 3.6; // potential growth rate of lateral roots, cm/day.
    const double rtran = 0.2; // the ratio of root mass transferred to a new soil
    // soil cell, when a lateral root grows into it.
//     On its initiation, lateral root length is assumed to be equal to the width 
//  of a soil column soil cell at the location of the taproot.
    int klocp1 = PlantRowColumn + 1;
    if (rlat2[l] <= 0)
        rlat2[l] = wk[klocp1];
    double stday; // daily average soil temperature (C) at root tip.
    stday = SoilTempDailyAvrg[l][klocp1] - 273.161;
    double temprg; // the effect of soil temperature on root growth.
    temprg = SoilTemOnRootGrowth(stday);
// define the column with the tip of this lateral root (ktip)
    int ktip = 0; // column with the tips of the laterals to the right
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
//     If the lateral reaches a new soil soil cell: a proportion (tran) of 
//	mass of roots is transferred to the new soil cell.
        if (rlat2[l] > sumwk && ktip < nk - 1) {
            int newktip; // column into which the tip of the lateral grows to left.
            newktip = ktip + 1; // column into which the tip of the lateral grows to left.
            for (int i = 0; i < NumRootAgeGroups; i++) {
                double tran = RootWeight[l][ktip][i] * rtran;
                RootWeight[l][ktip][i] -= tran;
                RootWeight[l][newktip][i] += tran;
            }
//     RootAge is initialized for this soil cell.
//     RootColNumLeft of this layer is redefined.
            if (RootAge[l][newktip] == 0)
                RootAge[l][newktip] = 0.01;
            if (newktip > RootColNumRight[l])
                RootColNumRight[l] = newktip;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void RootAging(int l, int k, double RootWeight[40][20][3], double RootAge[40][20])
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
    const double thtrn[2] = {20.0, 40.0}; // the time threshold, from the initial
// penetration of roots to a soil cell, after which some of the root
// mass of class i may be transferred into the next class (i+1).
    const double trn[2] = {0.0060, 0.0050}; //  the daily proportion of this transfer.
//
    double stday; // daily average soil temperature (c) of soil cell.
    stday = SoilTempDailyAvrg[l][k] - 273.161;
    RootAge[l][k] += SoilTemOnRootGrowth(stday);
//
    for (int i = 0; i < 2; i++) {
        if (RootAge[l][k] > thtrn[i]) {
            double xtr; // root mass transferred from one class to the next.
            xtr = trn[i] * RootWeight[l][k][i];
            RootWeight[l][k][i + 1] += xtr;
            RootWeight[l][k][i] -= xtr;
        }
    }
}

//////////////////////////////
double RootDeath(int l, int k, double DailyRootLoss, double RootWeight[40][20][3], const double RootAge[40][20]) {
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
    const double aa = 0.008; // a parameter in the equation for computing dthfac.
    const double dth[3] = {0.0001, 0.0002, 0.0001}; // the daily proportion of death of root tissue.
    const double dthmax = 0.10; // a parameter in the equation for computing dthfac.
    const double psi0 = -14.5; // a parameter in the equation for computing dthfac.
    const double thdth[3] = {30.0, 50.0, 100.0}; // the time threshold, from the initial
    // penetration of roots to a soil cell, after
    // which death of root tissue of class i may occur.
//
    for (int i = 0; i < 3; i++) {
        if (RootAge[l][k] > thdth[i]) {
            double dthfac; // the computed proportion of roots dying in each class.
            dthfac = dth[i];
            if (VolWaterContent[l][k] >= PoreSpace[l])
                dthfac = dthmax;
            else {
                if (i <= 1 && SoilPsi[l][k] <= psi0)
                    dthfac += aa * (psi0 - SoilPsi[l][k]);
                if (dthfac > dthmax)
                    dthfac = dthmax;
            }
            DailyRootLoss += RootWeight[l][k][i] * dthfac;
            RootWeight[l][k][i] -= RootWeight[l][k][i] * dthfac;
        }
    }
    return DailyRootLoss;
}

//////////////////////////////
double RootCultivation(int j, const int &NumRootAgeGroups, double DailyRootLoss, double RootWeight[40][20][3])
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
    int lcult = 0; // number of soil layers affected by cultivation.
    double sdpth = 0; // sum depth to the end of the layer.
    for (int l = 0; l < nl; l++) {
        sdpth += dl[l];
        if (sdpth >= CultivationDepth[j]) {
            lcult = l;
            break;
        }
    }
//
    double sumwk = 0; // sum of column widths from edge of slab to this column.
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
    return DailyRootLoss;
}

//////////////////////////////
void RootSummation(const string &ProfileName, const int &DayOfSimulation, const int &NumRootAgeGroups,
                   const int &NumLayersWithRoots, double RootWeight[40][20][3])
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
//     Compute the root weight (of all age classes) for all soil cells as
//  rootsv, and the total as roots.
    double roots = 0; // total weight of roots of all classes, g per slab.
    double rootsv[maxl][maxk]; // total dry weight of roots in a cell, g per soil cell.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++) {
            rootsv[l][k] = 0;
            for (int i = 0; i < 3; i++)
                rootsv[l][k] += RootWeight[l][k][i];
            roots += rootsv[l][k];
        }
//
    for (int l = 0; l < maxl; l++)
        for (int k = 0; k < maxk; k++)
            Scratch21[DayOfSimulation - 1].rootsv[l][k] = rootsv[l][k];
//     Convert total root weight from g per slab to g per plant.
    TotalRootWeight = roots * 100 * PerPlantArea / RowSpace;
//     If output flag is active, this flag also indicates the frequency of
//  output. If Kday is divisible by this frequency, compute sums of root
//  weight by classes (rwt), by layers (rootl), and by classes for each
//  layer (rootli).
    if (OutIndex[22] > 0) {
        int nni; // interval (days) between successive root outputs.
        nni = OutIndex[22];
        if ((Kday % nni) == 0 || Kday == 1) {
            double rootl[maxl]; // weight of all root classes in a soil layer, g.
            double rootli[maxl][3]; // weight of roots of class i in a soil layer, g.
            double rwt[3]; // total weight of roots of class i, g per slab.
            for (int i = 0; i < 3; i++)
                rwt[i] = 0;
            for (int l = 0; l < NumLayersWithRoots; l++) {
                rootl[l] = 0;
                for (int i = 0; i < 3; i++)
                    rootli[l][i] = 0;
                for (int k = 0; k < nk; k++)
                    for (int i = 0; i < 3; i++) {
                        rootl[l] += RootWeight[l][k][i];
                        rootli[l][i] += RootWeight[l][k][i];
                        rwt[i] += RootWeight[l][k][i];
                    }
            }
//    Write these computed results in file 34.
            ofstream File34(fs::path("output") / (ProfileName + ".RUT"), ios::app);
            File34 << endl;
            File34.unsetf(ios::left);
            File34 << " ROOT DATA FOR DAY Kday=";
            File34.width(5);
            File34 << Kday;
            File34 << endl << endl;
            File34 << " Layer    Total     Class=     1              2              3";
            File34 << endl;
//
            for (int l = 0; l < NumLayersWithRoots; l++) {
                File34.width(5);
                File34 << l + 1 << " ";
                File34.setf(ios::scientific);
                File34.precision(3);
                File34.width(14);
                File34 << rootl[l] << " ";
                for (int i = 0; i < NumRootAgeGroups; i++) {
                    File34.precision(3);
                    File34.width(14);
                    File34 << rootli[l][i] << " ";
                }
                File34 << endl;
            }
            File34 << endl;
//
            File34 << " Root Weight per Plant = ";
            File34.width(15);
            File34 << TotalRootWeight << endl;
//
            File34 << " Root Weight of each class per slab = " << endl;
            for (int i = 0; i < NumRootAgeGroups; i++) {
                File34.precision(3);
                File34.width(13);
                File34 << "  " << rwt[i];
            }
            File34 << endl << endl;
//    Compute and write (to file 34) the roots capable of uptake in mg
//  per cm**3 soil volume for each cell.
            File34 << "  Root weight capable of uptake (in mg / cm**3 soil volume) :" << endl;
            File34 << " Layer ";
            for (int i = 0; i < nk; i++) {
                File34.width(9);
                File34 << i + 1;
            }
            File34 << endl;

            for (int l = 0; l < NumLayersWithRoots; l++) {
                File34.width(5);
                File34 << l + 1 << " ";
                for (int k = 0; k < nk; k++) {
                    File34.setf(ios::scientific);
                    File34.precision(2);
                    File34.width(8);
                    File34 << 1000 * RootWtCapblUptake[l][k] / (dl[l] * wk[k]) << " ";
                }
                File34 << endl;
            }
        }
    }
}
