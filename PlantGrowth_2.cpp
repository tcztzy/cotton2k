// File  PlantGrowth_2.cpp
//
//   functions in this file:
// PotentialStemGrowth()
// PotentialLeafGrowth()
// TemperatureOnLeafGrowthRate()
// TemperatureOnFruitGrowthRate()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//////////////////////////////////////////////////////////////////
double PotentialStemGrowth(double stemnew)
//     This function computes and returns the potential stem growth of cotton
// plants. It is called from PlantGrowth().
//     The following argument is used:
//        stemnew - dry weight of active stem tissue.
//     The following global variables are referenced here:
//        DensityFactor, FruitingCode, Kday, VarPar.
{
    double potentialStemGrowth;
    //     There are two periods for computation of potential stem growth:
    //     (1) Before the appearance of a square on the third fruiting
    // branch. Potential stem growth is a functon of plant age (Kday, days
    // from emergence).
    if (FruitingCode[0][2][0] == 0)
        potentialStemGrowth = VarPar[12] * (VarPar[13] + VarPar[14] * Kday);
    //     (2) After the appearance of a square on the third fruiting
    // branch. It is assumed that all stem tissue that is more than 32
    // days old is not active. Potential stem growth is a function
    // of active stem tissue weight (stemnew), and plant density (denfac).
    else {
        double denfac;  // effect of plant density on stem growth rate.
        denfac = 1 - VarPar[15] * (1 - DensityFactor);
        if (denfac < 0.2) denfac = 0.2;
        potentialStemGrowth =
            denfac * VarPar[16] * (VarPar[17] + VarPar[18] * stemnew);
    }
    return potentialStemGrowth;
}
//////////////////////////
void PotentialLeafGrowth()
//     This function simulates the potential growth of leaves of cotton plants.
//  It is called from PlantGrowth(). It calls function
//  TemperatureOnLeafGrowthRate().
//
//     The following monomolecular growth function is used :
//       leaf area = smax * (1 - exp(-c * pow(t,p)))
//     where    smax = maximum leaf area.
//              t = time (leaf age).
//              c, p = constant parameters.
//     Note: p is constant for all leaves, whereas smax and c depend on leaf
//     position. The rate per day (the derivative of this function) is :
//           r = smax * c * p * exp(-c * pow(t,p)) * pow(t, (p-1))
//
//     The following global variables are referenced here:
//        AgeOfPreFruNode, AvrgDailyTemp, DensityFactor, LeafAge, LeafAreaNodes,
//        LeafAreaMainStem, LeafAreaPreFru, NumFruitBranches, NumNodes,
//        NumPreFruNodes, NumVegBranches, pixda, VarPar, WaterStress.
//     The following global variables are set here:
//        LeafWeightAreaRatio, PotGroAllLeaves, PotGroAllPetioles,
//        PotGroLeafAreaMainStem, PotGroLeafAreaNodes,  PotGroLeafAreaPreFru,
//        PotGroLeafWeightMainStem, PotGroLeafWeightNodes,
//        PotGroLeafWeightPreFru, PotGroPetioleWeightMainStem,
//        PotGroPetioleWeightNodes, PotGroPetioleWeightPreFru.
{
    //     The following constant parameters. are used in this function:
    const double p = 1.6;  // parameter of the leaf growth rate equation.
    const double vpotlf[14] = {3.0,      0.95,       1.2,   13.5,    -0.62143,
                               0.109365, 0.00137566, 0.025, 0.00005, 30.,
                               0.02,     0.001,      2.50,  0.18};
    //     Calculate water stress reduction factor for leaf growth rate
    // (wstrlf). This has been empirically calibrated in COTTON2K.
    double wstrlf =
        WaterStress * (1 + vpotlf[0] * (2 - WaterStress)) - vpotlf[0];
    if (wstrlf < 0.05) wstrlf = 0.05;
    //     Calculate wtfstrs, the effect of leaf water stress
    //     onLeafWeightAreaRatio (the ratio
    //  of leaf dry weight to leaf area). This has also been empirically
    //  calibrated in COTTON2K.
    double wtfstrs = vpotlf[1] + vpotlf[2] * (1 - wstrlf);
    //     Compute the ratio of leaf dry weight increment to leaf area increment
    //     (g per dm2),
    //  as a function of average daily temperature and water stress. Parameters
    //  for the effect of temperature are adapted from GOSSYM.
    double tdday =
        AvrgDailyTemp;  // limited value of today's average temperature.
    if (tdday < vpotlf[3]) tdday = vpotlf[3];
    LeafWeightAreaRatio =
        wtfstrs / (vpotlf[4] + tdday * (vpotlf[5] - tdday * vpotlf[6]));
    //     Assign zero to total potential growth of leaf and petiole.
    PotGroAllLeaves = 0;
    PotGroAllPetioles = 0;
    double c = 0;     // parameter of the leaf growth rate equation.
    double smax = 0;  // maximum possible leaf area, a parameter of the leaf
                      // growth rate equation.
    double rate;      // growth rate of area of a leaf.
    //     Compute the potential growth rate of prefruiting leaves.
    //  smax and c are functions of prefruiting node number.
    for (int j = 0; j < NumPreFruNodes; j++)  // loop by prefruiting node.
    {
        if (LeafAreaPreFru[j] <= 0) {
            PotGroLeafAreaPreFru[j] = 0;
            PotGroLeafWeightPreFru[j] = 0;
            PotGroPetioleWeightPreFru[j] = 0;
        } else {
            int jp1 = j + 1;
            smax = jp1 * (VarPar[2] - VarPar[3] * jp1);
            if (smax < VarPar[4]) smax = VarPar[4];
            c = vpotlf[7] + vpotlf[8] * jp1 * (jp1 - vpotlf[9]);
            rate = smax * c * p * exp(-c * pow(AgeOfPreFruNode[j], p)) *
                   pow(AgeOfPreFruNode[j], (p - 1));
            //     Growth rate is modified by water stress, pix and a function
            //     of average temperature. Compute potential growth of leaf
            //     area, leaf weight and petiole
            //  weight for leaf on node j. Add leaf weight potential growth to
            //  PotGroAllLeaves. Add potential growth of petiole weight to
            //  PotGroAllPetioles.
            if (rate >= 1e-12) {
                PotGroLeafAreaPreFru[j] =
                    rate * wstrlf * pixda *
                    TemperatureOnLeafGrowthRate(AvrgDailyTemp);
                PotGroLeafWeightPreFru[j] =
                    PotGroLeafAreaPreFru[j] * LeafWeightAreaRatio;
                PotGroPetioleWeightPreFru[j] =
                    PotGroLeafAreaPreFru[j] * LeafWeightAreaRatio * vpotlf[13];
                PotGroAllLeaves += PotGroLeafWeightPreFru[j];
                PotGroAllPetioles += PotGroPetioleWeightPreFru[j];
            }  // rate
        }      // LeafAreaPreFru
    }          // NumPreFruNodes
               //     denfac is the effect of plant density on leaf growth rate.
    double denfac = 1 - vpotlf[12] * (1 - DensityFactor);
    int nbrch;  // number of fruiting branches on a vegetative stem.
    int nnid;   // number of nodes on a fruiting branch.
                //    Loop for all vegetative stems.
    for (int k = 0; k < NumVegBranches; k++)  // loop of vegetative branches
    {
        nbrch = NumFruitBranches[k];
        for (int l = 0; l < nbrch; l++)  // loop of fruiting branches
        {
            //     smax and c are  functions of fruiting branch number.
            //     smax is modified by plant density, using the density factor
            //     denfac. Compute potential main stem leaf growth, assuming
            //     that the main
            //  stem leaf is initiated at the same time as leaf (k,l,0).
            if (LeafAreaMainStem[k][l] <= 0) {
                PotGroLeafAreaMainStem[k][l] = 0;
                PotGroLeafWeightMainStem[k][l] = 0;
                PotGroPetioleWeightMainStem[k][l] = 0;
            } else {
                int lp1 = l + 1;
                smax = VarPar[5] + VarPar[6] * lp1 * (VarPar[7] - lp1);
                smax = smax * denfac;
                if (smax < VarPar[4]) smax = VarPar[4];
                c = vpotlf[10] + lp1 * vpotlf[11];
                if (LeafAge[k][l][0] > 70)
                    rate = 0;
                else
                    rate = smax * c * p * exp(-c * pow(LeafAge[k][l][0], p)) *
                           pow(LeafAge[k][l][0], (p - 1));
                //     Add leaf and petiole weight potential growth to SPDWL and
                //     SPDWP.
                if (rate >= 1e-12) {
                    PotGroLeafAreaMainStem[k][l] =
                        rate * wstrlf * pixda *
                        TemperatureOnLeafGrowthRate(AvrgDailyTemp);
                    PotGroLeafWeightMainStem[k][l] =
                        PotGroLeafAreaMainStem[k][l] * LeafWeightAreaRatio;
                    PotGroPetioleWeightMainStem[k][l] =
                        PotGroLeafAreaMainStem[k][l] * LeafWeightAreaRatio *
                        vpotlf[13];
                    PotGroAllLeaves += PotGroLeafWeightMainStem[k][l];
                    PotGroAllPetioles += PotGroPetioleWeightMainStem[k][l];
                }  // rate
            }      // LeafAreaMainStem
            //     Assign smax value of this main stem leaf to smaxx, c to cc.
            //     Loop over the nodes of this fruiting branch.
            double smaxx =
                smax;  // value of smax for the corresponding main stem leaf.
            double cc = c;  // value of c for the corresponding main stem leaf.
            nnid = NumNodes[k][l];
            for (int m = 0; m < nnid;
                 m++)  // loop of nodes on a fruiting branch
            {
                if (LeafAreaNodes[k][l][m] <= 0) {
                    PotGroLeafAreaNodes[k][l][m] = 0;
                    PotGroLeafWeightNodes[k][l][m] = 0;
                    PotGroPetioleWeightNodes[k][l][m] = 0;
                }
                //     Compute potential growth of leaf area and leaf weight for
                //     leaf
                //  on fruiting branch node (k,l,m).
                //    Add leaf and petiole weight potential growth to spdwl and
                //    spdwp.
                else {
                    int mp1 = m + 1;
                    //     smax and c are reduced as a function of node number
                    //     on this fruiting branch.
                    smax = smaxx * (1 - VarPar[8] * mp1);
                    c = cc * (1 - VarPar[8] * mp1);
                    //     Compute potential growth for the leaves on fruiting
                    //     branches.
                    if (LeafAge[k][l][m] > 70)
                        rate = 0;
                    else
                        rate = smax * c * p *
                               exp(-c * pow(LeafAge[k][l][m], p)) *
                               pow(LeafAge[k][l][m], (p - 1));
                    if (rate >= 1e-12) {
                        //     Growth rate is modified by water stress and pix.
                        //     Potential growth
                        //  is computed as a function of average temperature.
                        PotGroLeafAreaNodes[k][l][m] =
                            rate * wstrlf * pixda *
                            TemperatureOnLeafGrowthRate(AvrgDailyTemp);
                        PotGroLeafWeightNodes[k][l][m] =
                            PotGroLeafAreaNodes[k][l][m] * LeafWeightAreaRatio;
                        PotGroPetioleWeightNodes[k][l][m] =
                            PotGroLeafAreaNodes[k][l][m] * LeafWeightAreaRatio *
                            vpotlf[13];
                        PotGroAllLeaves += PotGroLeafWeightNodes[k][l][m];
                        PotGroAllPetioles += PotGroPetioleWeightNodes[k][l][m];
                    }  // if rate
                }      // if LeafAreaNodes
            }          // loop nnid
        }              // loop nbrch
    }                  // loop NumVegBranches
}
/////////////////////////////////////////////////////////////////////////
double TemperatureOnLeafGrowthRate(double t) {
    //     This is the temperature function for leaf growth rate. It is called
    //     by function
    //  PotentialLeafGrowth(). It is based on the original code of GOSSYM, and
    //  the parameters are the same. The argument t is air temperature (C).
    //     ra is divided by the maximum value. Thus, the function returns values
    //     between 0. and 1.
    //  This will result :
    //     maximum value of TemperatureOnLeafGrowthRate = 1.00  at  t = 29.86747
    //         for t = 24    TemperatureOnLeafGrowthRate = 0.766
    //         for t = 27    TemperatureOnLeafGrowthRate = 0.953
    //         for t = 30    TemperatureOnLeafGrowthRate = 0.999
    //         for t = 36    TemperatureOnLeafGrowthRate = 0.737
    //         for t = 42    TemperatureOnLeafGrowthRate = 0.0
    //      and for t <= 24 :
    //         for t = 12    TemperatureOnLeafGrowthRate = 0.0
    //         for t = 16    TemperatureOnLeafGrowthRate = 0.269
    //         for t = 20    TemperatureOnLeafGrowthRate = 0.549
    //         for t = 24    TemperatureOnLeafGrowthRate = 0.768
    //
    //     Constant parameters used:
    const double par[8] = {24,        -1.14277,  0.0910026,   0.00152344,
                           -0.317136, 0.0300712, 0.000416356, 0.2162044};
    //
    double ra;  // intermediate value for computing temperatureOnLeafGrowthRate.
    if (t > par[0])
        ra = par[1] + t * (par[2] - t * par[3]);
    else
        ra = par[4] + t * (par[5] - t * par[6]);
    //
    if (ra < 0) ra = 0;
    double temperatureOnLeafGrowthRate = ra / par[7];
    return temperatureOnLeafGrowthRate;
}
////////////////////////////////////
double TemperatureOnFruitGrowthRate(double t)
//     This function computes the effect of air temperature (t) on growth
//  rate of bolls in cotton plants. It is called from PotentialFruitGrowth().
//     Some values computed by this function:
//   t (C)       tfr
//   12          0.
//   15          0.336
//   20          0.751
//   25          0.978
//   26          1.
//   28.5        1.024 (maximum)
//   30          1.016
//   35          0.866
//   40          0.527
//   45          0.
//
{
    const double p1 = -2.041;
    const double p2 = 0.215;
    const double p3 = 0.00377;
    //
    double tfr = p1 + t * (p2 - p3 * t);
    if (tfr < 0) tfr = 0;
    return tfr;
}
