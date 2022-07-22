//  RootGrowth_1.cpp
//
//   functions in this file:
// RootImpedance()
// SoilMechanicResistance()
// SoilAirOnRootGrowth()
// SoilNitrateOnRootGrowth()
// SoilWaterOnRootGrowth()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//////////////////////////
void RootImpedance()
//     This function calculates soil mechanical impedance to root growth,
//     rtimpd(l,k),
//  for all soil cells. It is called from PotentialRootGrowth(). The impedance
//  is a function of bulk density and water content in each soil soil cell. No
//  changes have been made in the original GOSSYM code.
//
//     The following global variables are referenced here:
//  BulkDensity, gh2oc, impede, inrim, ncurve, nk, nl,
//  SoilHorizonNum, tstbd, VolWaterContent.
//     The following global variables are set here:    RootImpede.
{
    for (int l = 0; l < nl; l++) {
        int j = SoilHorizonNum[l];
        double Bd = BulkDensity[j];  // bulk density for this layer
                                     //
        int jj;
        for (jj = 0; jj < inrim; jj++) {
            if (Bd <= tstbd[jj][0]) break;
        }
        int j1 = jj;
        if (j1 > inrim - 1) j1 = inrim - 1;
        int j0 = jj - 1;
        if (j0 < 0) j0 = 0;
        //
        for (int k = 0; k < nk; k++) {
            double Vh2o = VolWaterContent[l][k] / Bd;
            int ik;
            for (ik = 0; ik < ncurve; ik++) {
                if (Vh2o <= gh2oc[ik]) break;
            }
            int i1 = ik;
            if (i1 > ncurve - 1) i1 = ncurve - 1;
            int i0 = ik - 1;
            if (i0 < 0) i0 = 0;
            //
            if (j1 == 0) {
                if (i1 == 0 || Vh2o <= gh2oc[i1])
                    RootImpede[l][k] = impede[j1][i1];
                else
                    RootImpede[l][k] =
                        impede[j1][i0] - (impede[j1][i0] - impede[j1][i1]) *
                                             (Vh2o - gh2oc[i0]) /
                                             (gh2oc[i1] - gh2oc[i0]);
            } else {
                if (i1 == 0 || Vh2o <= gh2oc[i1])
                    RootImpede[l][k] =
                        impede[j0][i1] - (impede[j0][i1] - impede[j1][i1]) *
                                             (tstbd[j0][i1] - Bd) /
                                             (tstbd[j0][i1] - tstbd[j1][i1]);
                else {
                    double temp1 =
                        impede[j0][i1] - (impede[j0][i1] - impede[j1][i1]) *
                                             (tstbd[j0][i1] - Bd) /
                                             (tstbd[j0][i1] - tstbd[j1][i1]);
                    double temp2 =
                        impede[j0][i0] - (impede[j0][i0] - impede[j1][i1]) *
                                             (tstbd[j0][i0] - Bd) /
                                             (tstbd[j0][i0] - tstbd[j1][i0]);
                    RootImpede[l][k] = temp2 + (temp1 - temp2) *
                                                   (Vh2o - gh2oc[i0]) /
                                                   (gh2oc[i1] - gh2oc[i0]);
                }
            }
        }
    }
    //
}
//////////////////////////
double SoilMechanicResistance(int l, int k)
//     This function calculates soil mechanical resistance of cell l,k. It is
//     computed
//  on the basis of parameters read from the input and calculated in
//  RootImpedance().
//     It is called from PotentialRootGrowth().
//     The function has been adapted, without change, from the code of GOSSYM.
//     Soil mechanical
//  resistance is computed as an empirical function of bulk density and water
//  content. It should be noted, however, that this empirical function is based
//  on data for one type of soil only, and its applicability for other soil
//  types is questionable. The effect of soil moisture is only indirectly
//  reflected in this function. A new module (SoilWaterOnRootGrowth) has
//  therefore been added in COTTON2K to simulate an additional direct effect of
//  soil moisture on root growth.
//
//     The following global variables are referenced here:
//       nk, nl, rtimpd
//     The following arguments are used:
//       k, l - column and layer numbers pf this cell.
//
//     The minimum value of rtimpd of this and neighboring soil cells is used to
//     compute
//  rtpct. The code is based on a segment of RUTGRO in GOSSYM, and the values of
//  the p1 to p3 parameters are based on GOSSYM usage:
{
    const double p1 = 1.046;
    const double p2 = 0.034554;
    const double p3 = 0.5;
    //
    int lp1;          // layer below l.
    if (l == nl - 1)  // last layer
        lp1 = l;
    else
        lp1 = l + 1;
    //
    int km1, kp1;     // columns to the left and to the right of k.
    if (k == nk - 1)  // last column
        kp1 = k;
    else
        kp1 = k + 1;
    //
    if (k == 0)  // first column
        km1 = 0;
    else
        km1 = k - 1;
    //
    double rtimpd0 = RootImpede[l][k];
    double rtimpdkm1 = RootImpede[l][km1];
    double rtimpdkp1 = RootImpede[l][kp1];
    double rtimpdlp1 = RootImpede[lp1][k];
    //
    double rtimpdmin = rtimpd0;  // minimum value of rtimpd
    if (rtimpdkm1 < rtimpdmin) rtimpdmin = rtimpdkm1;
    if (rtimpdkp1 < rtimpdmin) rtimpdmin = rtimpdkp1;
    if (rtimpdlp1 < rtimpdmin) rtimpdmin = rtimpdlp1;
    //
    double rtpct;  // effect of soil mechanical resistance on root growth (the
                   // return value).
    rtpct = p1 - p2 * rtimpdmin;
    if (rtpct > 1) rtpct = 1;
    if (rtpct < p3) rtpct = p3;
    //
    return rtpct;
}
//////////////////////////
double SoilAirOnRootGrowth(double psislk, double poreSpace, double vh2oclk)
//     This function calculates the reduction of potential root growth rate in
//     cells with
//  low oxygen content (high water content). It is called from
//  PotentialRootGrowth().
//     It has been adapted from GOSSYM, but the critical value of soil moisture
//     potential
//  for root growth reduction (i.e., water logging conditions) has been changed.
//
//     The following input arguments are used:
//        poreSpace -  value of PoreSpace (v/v) for this layer.
//        psislk -  value of SoilPsi for this cell.
//        vh2oclk - water content (v/v) of this cell
{
    //     Constant parameters:
    double p1 = 0;
    double p2 = 1;
    double p3 = 0.1;
    //     The following is actually disabled by the choice of the calibration
    //     parameters. It
    //  may be redefined when more experimental data become available.
    double rtrdo;  // Effect of oxygen deficiency on root growth (the return
                   // value).
    if (psislk > p1)
        rtrdo = p2;
    else
        rtrdo = 1;
    //   Reduced root growth when water content is at pore - space saturation
    //  (below water table).
    if (vh2oclk >= poreSpace) rtrdo = p3;
    return rtrdo;
}
//////////////////////////
double SoilNitrateOnRootGrowth(double vno3clk) {
    //     This function calculates the reduction of potential root growth rate
    //     in cells with
    //  low nitrate content. It is called from PotentialRootGrowth().
    //     It has been adapted from GOSSYM. It is assumed that root growth is
    //     reduced when
    //  nitrate N content falls below a certain level.
    //
    //     The following argument is used:
    //        vno3clk - VolNo3NContent value for this cell
    //
    //     The following constant parameters are used:
    const double p1 = 0;
    const double p2 = 1;
    //     Note: This function actually does nothing. It is disabled by the
    //     choice of the constant
    //  parameters. It may be redefined when more experimental data become
    //  available.
    double rtrdn;  // effect of nitrate deficiency on root growth (the return
                   // value).
    if (vno3clk < p1)
        rtrdn = p2;
    else
        rtrdn = 1;
    return rtrdn;
}
//////////////////////////
double SoilWaterOnRootGrowth(double psislk)
//     This function returns the effect of soil  moisture in cell l,k on cotton
//     root potential
// growth rate. It is called from PotentialRootGrowth() and uses the matric
// potential of this cell.
//
//     The following argument is used:
//        psislk - soil water potential (bars) of this cell.
//
{
    //     The following constants are used:
    const double p1 = 20;
    const double p2 = 16;
    //     It is assumed that almost no root growth occurs when the soil is
    //  dryer than -p1 (-20 bars), and root growth rate is maximum at a matric
    //  potential of -4 bars (p2 - p1) or wetter.
    double smf;  // effect of soil moisture on root growth (the return value).
    //     smf is computed here as an empirical third degree function,
    //  with values between 0.02 and 1.
    smf = pow(((p1 + psislk) / p2), 3);
    if (smf < 0.02) smf = 0.02;
    if (smf > 1) smf = 1;
    return smf;
}
