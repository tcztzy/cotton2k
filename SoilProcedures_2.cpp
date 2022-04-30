// File SoilProcedures_2.cpp
//
//   functions in this file:
// CapillaryFlow()
// Drain()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//////////////////////////
void CapillaryFlow()
//     This function computes the capillary water flow between soil cells. It is
//     called by
//  SoilProcedures(), noitr times per day.  The number of iterations (noitr) has
//  been computed in SoilProcedures() as a function of the amount of water
//  applied. It is executed only once per day if no water is applied by rain or
//  irrigation.
//     It calls functions:   Drain(), NitrogenFlow(), psiq(), PsiOsmotic(),
//     WaterFlux().
//
//     The following global variables are referenced:
//       alpha, ,beta, Daynum, DayStart, dl, ElCondSatSoilToday, nk, nl,
//       PoreSpace, RowSpace, SoilHorizonNum, thad ,thts, WaterTableLayer, wk.
//     The following global variables are set:
//       CumWaterDrained, SoilPsi, VolNo3NContent, VolUreaNContent,
//       VolWaterContent.
{
    static long numiter;  // counter used for WaterFlux() calls.
    double wk1[40];       // dummy array for passing values of array wk.
                          //     Set initial values in first day.
    if (Daynum <= DayStart) {
        numiter = 0;
        for (int l = 0; l < nl; l++) wk1[l] = 0;
    }
    //     Increase the counter numiter, and compute the updated values of
    //     SoilPsi in each
    //  soil cell by calling functions psiq() and PsiOsmotic().
    numiter++;
    for (int l = 0; l < nl; l++) {
        int j = SoilHorizonNum[l];  //  the soil horizon number
        for (int k = 0; k < nk; k++)
            SoilPsi[l][k] =
                psiq(VolWaterContent[l][k], thad[l], thts[l], alpha[j],
                     beta[j]) -
                PsiOsmotic(VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
    }
    //
    int nlx = nl;  // The last layer without a water table.
    if (WaterTableLayer < nlx) nlx = WaterTableLayer - 1;
    int iv;  //  direction indicator: iv = 1 for vertical flow in each column;
             //    iv = 0 for horizontal flow in each layer.
    double q01[40];  // one dimensional array of a layer or a column of previous
                     // values of VolWaterContent.
    double q1[40];   // one dimensional array of a layer or a column of
                     // VolWaterContent.
    double
        psi1[40];    // one dimensional array of a layer or a column of SoilPsi.
    double nit[40];  // one dimensional array of a layer or a column of
                     // VolNo3NContent.
    double nur[40];  // one dimensional array of a layer or a column of
                     // VolUreaNContent.
                     //
    //     VERTICAL FLOW in each column. the direction indicator iv is set to 1.
    iv = 1;
    //     Loop over all columns. Temporary one-dimensional arrays are defined
    //     for each column:
    //  assign the VolWaterContent[] values to temporary one-dimensional arrays
    //  q1 and q01. Assign SoilPsi, VolNo3NContent and VolUreaNContent values to
    //  arrays psi1, nit and nur, respectively.
    for (int k = 0; k < nk; k++) {
        for (int l = 0; l < nlx; l++) {
            q1[l] = VolWaterContent[l][k];
            q01[l] = VolWaterContent[l][k];
            psi1[l] = SoilPsi[l][k] + PsiOsmotic(VolWaterContent[l][k], thts[l],
                                                 ElCondSatSoilToday);
            nit[l] = VolNo3NContent[l][k];
            nur[l] = VolUreaNContent[l][k];
        }  // end loop l
           //     Call the following functions: WaterFlux() calculates the water
           //     flow caused by potential
        //  gradients; NitrogenFlow() computes the movement of nitrates caused
        //  by the flow of water.
        WaterFlux(q1, psi1, dl, thad, thts, PoreSpace, nlx, iv, 0, numiter);
        NitrogenFlow(nl, q01, q1, dl, nit, nur);
        //     Reassign the updated values of q1, nit, nur and psi1 back to
        //  VolWaterContent, VolNo3NContent, VolUreaNContent and SoilPsi.
        for (int l = 0; l < nlx; l++) {
            VolWaterContent[l][k] = q1[l];
            VolNo3NContent[l][k] = nit[l];
            VolUreaNContent[l][k] = nur[l];
            SoilPsi[l][k] = psi1[l] - PsiOsmotic(VolWaterContent[l][k], thts[l],
                                                 ElCondSatSoilToday);
        }            // end loop l
    }                // end loop k
    double pp1[40];  // one dimensional array of a layer or a column of PP.
    double qr1[40];  // one dimensional array of a layer or a column of THAD.
    double qs1[40];  // one dimensional array of a layer or a column of THTS.
                     //
    //     HORIZONTAL FLUX in each layer. The direction indicator iv is set to
    //     0.
    iv = 0;
    //     Loop over all layers. Define the horizon number j for this layer.
    //     Temporary
    //  one-dimensional arrays are defined for each layer: assign the
    //  VolWaterContent values to  q1 and q01. Assign SoilPsi, VolNo3NContent,
    //  VolUreaNContent, thad and thts values of the soil cells to arrays psi1,
    //  nit, nur, qr1 and qs1, respectively.
    for (int l = 0; l < nlx; l++) {
        for (int k = 0; k < nk; k++) {
            q1[k] = VolWaterContent[l][k];
            q01[k] = VolWaterContent[l][k];
            psi1[k] = SoilPsi[l][k] + PsiOsmotic(VolWaterContent[l][k], thts[l],
                                                 ElCondSatSoilToday);
            qr1[k] = thad[l];
            qs1[k] = thts[l];
            pp1[k] = PoreSpace[l];
            nit[k] = VolNo3NContent[l][k];
            nur[k] = VolUreaNContent[l][k];
            wk1[k] = wk[k];
        }
        //     Call subroutines WaterFlux(), and NitrogenFlow() to compute water
        //     nitrate and
        //  urea transport in the layer.
        WaterFlux(q1, psi1, wk1, qr1, qs1, pp1, nk, iv, l, numiter);
        NitrogenFlow(nk, q01, q1, wk1, nit, nur);
        //     Reassign the updated values of q1, nit, nur and psi1 back to
        //  VolWaterContent, VolNo3NContent, VolUreaNContent and SoilPsi.
        for (int k = 0; k < nk; k++) {
            VolWaterContent[l][k] = q1[k];
            SoilPsi[l][k] = psi1[k] - PsiOsmotic(VolWaterContent[l][k], thts[l],
                                                 ElCondSatSoilToday);
            VolNo3NContent[l][k] = nit[k];
            VolUreaNContent[l][k] = nur[k];
        }  // snd loop k
    }      // end loop l
    //     Call Drain() to move excess water down in the column and compute
    //     drainage out
    //  of the column. Update cumulative drainage.
    double WaterDrainedOut = 0;  // water drained out of the slab, mm.
    WaterDrainedOut += Drain();
    if (WaterDrainedOut > 0) CumWaterDrained += 10 * WaterDrainedOut / RowSpace;
    //  Compute the soil water potential for all soil cells.
    for (int l = 0; l < nl; l++) {
        int j = SoilHorizonNum[l];
        for (int k = 0; k < nk; k++) {
            SoilPsi[l][k] =
                psiq(VolWaterContent[l][k], thad[l], thts[l], alpha[j],
                     beta[j]) -
                PsiOsmotic(VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////
double Drain()
//     This function computes the gravity flow of water in the slab, and returns
//     the
//  drainage of water out of the slab. It is called from GravityFlow() and
//  CapillaryFlow().
//
//     The following global variables are referenced:
//       dl, FieldCapacity, MaxWaterCapacity, nk, nl, NO3FlowFraction,
//       PoreSpace, RowSpace, WaterTableLayer, wk
//     The following global variables are set:
//       SoilNitrogenLoss, VolNo3NContent, VolUreaNContent, VolWaterContent,
{
    int nlx = nl;  // last soil layer for computing drainage.
    if (WaterTableLayer < nlx) nlx = WaterTableLayer;
    double oldvh2oc[20];  // stores previous values of VolWaterContent.
    double nitconc;       // nitrate N concentration in the soil solution.
    double nurconc;       // urea N concentration in the soil solution.
    //     The following is executed if this is not the bottom layer.
    for (int l = 0; l < nlx - 1; l++) {
        //     Compute the average water content (avwl) of layer l. Store the
        //  water content in array oldvh2oc.
        double avwl = 0;  // average water content in a soil layer
        for (int k = 0; k < nk; k++) {
            avwl += VolWaterContent[l][k] * wk[k] / RowSpace;
            oldvh2oc[k] = VolWaterContent[l][k];
        }
        //     Upper limit of water content in free drainage..
        double uplimit = MaxWaterCapacity[l];
        //
        //     Check if the average water content exceeds uplimit for this
        //     layer, and if it does,
        //  compute amount (wmov) to be moved to the next layer from each cell.
        double wmov;  // amount of water moving out of a cell.
        if (avwl > uplimit) {
            wmov = avwl - uplimit;
            wmov = wmov * dl[l] / dl[l + 1];
            for (int k = 0; k < nk; k++) {
                //     Water content of all soil cells in this layer will be
                //     uplimit. the amount (qmv)
                //  to be added to each cell of the next layer is computed
                //  (corrected for non uniform column widths). The water content
                //  in the next layer is computed.
                VolWaterContent[l][k] = uplimit;
                VolWaterContent[l + 1][k] += wmov * wk[k] * nk / RowSpace;
                //     The concentrations of nitrate and urea N in the soil
                //     solution are
                //  computed and their amounts in this layer and in the next one
                //  are updated.
                double qvout;  // amount of water moving out of a cell.
                qvout = (oldvh2oc[k] - uplimit);
                if (qvout > 0) {
                    nitconc = VolNo3NContent[l][k] / oldvh2oc[k];
                    if (nitconc < 1.e-30) nitconc = 0;
                    nurconc = VolUreaNContent[l][k] / oldvh2oc[k];
                    if (nurconc < 1.e-30) nurconc = 0;
                    VolNo3NContent[l][k] = VolWaterContent[l][k] * nitconc;
                    VolUreaNContent[l][k] = VolWaterContent[l][k] * nurconc;
                    //     Only a part ( NO3FlowFraction ) of N is moved with
                    //     water draining.
                    double
                        vno3mov;  // amount of nitrate N moving out of a cell.
                    double vnurmov;  // amount of urea N moving out of a cell.
                    vno3mov = qvout * nitconc;
                    VolNo3NContent[l + 1][k] +=
                        NO3FlowFraction[l] * vno3mov * dl[l] / dl[l + 1];
                    VolNo3NContent[l][k] += (1 - NO3FlowFraction[l]) * vno3mov;
                    vnurmov = qvout * nurconc;
                    VolUreaNContent[l + 1][k] +=
                        NO3FlowFraction[l] * vnurmov * dl[l] / dl[l + 1];
                    VolUreaNContent[l][k] += (1 - NO3FlowFraction[l]) * vnurmov;
                }
            }
        }
        //     If the average water content is not higher than uplimit, start
        //     another loop over columns.
        else {
            for (int k = 0; k < nk; k++) {
                //  Check each soil cell if water content exceeds uplimit,
                if (VolWaterContent[l][k] > uplimit) {
                    wmov = VolWaterContent[l][k] - uplimit;
                    VolWaterContent[l][k] = uplimit;
                    VolWaterContent[l + 1][k] += wmov * dl[l] / dl[l + 1];
                    nitconc = VolNo3NContent[l][k] / oldvh2oc[k];
                    if (nitconc < 1.e-30) nitconc = 0;
                    nurconc = VolUreaNContent[l][k] / oldvh2oc[k];
                    if (nurconc < 1.e-30) nurconc = 0;
                    VolNo3NContent[l][k] = VolWaterContent[l][k] * nitconc;
                    VolUreaNContent[l][k] = VolWaterContent[l][k] * nurconc;
                    //
                    VolNo3NContent[l + 1][k] +=
                        NO3FlowFraction[l] * wmov * nitconc * dl[l] / dl[l + 1];
                    VolUreaNContent[l + 1][k] +=
                        NO3FlowFraction[l] * wmov * nurconc * dl[l] / dl[l + 1];
                    VolNo3NContent[l][k] +=
                        (1 - NO3FlowFraction[l]) * wmov * nitconc;
                    VolUreaNContent[l][k] +=
                        (1 - NO3FlowFraction[l]) * wmov * nurconc;
                }  // end if Vol...
            }      // end loop k
        }          // end if avwl...
    }              // end loop l
                   //     For the lowermost soil layer, loop over columns:
    //     It is assumed that the maximum amount of water held at the lowest
    //     soil layer (nlx-1)
    //  of the slab is equal to FieldCapacity. If water content exceeds
    //  MaxWaterCapacity, compute the water drained out (Drainage), update
    //  water, nitrate and urea, compute nitrogen lost by drainage, and add it
    //  to the cumulative N loss SoilNitrogenLoss.
    double Drainage;  // drainage of water out of the slab, cm3 (return value)
    Drainage = 0;
    for (int k = 0; k < nk; k++) {
        if (VolWaterContent[nlx - 1][k] > MaxWaterCapacity[nlx - 1]) {
            Drainage +=
                (VolWaterContent[nlx - 1][k] - MaxWaterCapacity[nlx - 1]) *
                dl[nlx - 1] * wk[k];
            nitconc = VolNo3NContent[nlx - 1][k] / oldvh2oc[k];
            if (nitconc < 1.e-30) nitconc = 0;
            nurconc = VolUreaNContent[nlx - 1][k] / oldvh2oc[k];
            if (nurconc < 1.e-30) nurconc = 0;
            double saven;  //  intermediate variable for computing N loss.
            saven = (VolNo3NContent[nlx - 1][k] + VolUreaNContent[nlx - 1][k]) *
                    dl[nlx - 1] * wk[k];
            VolWaterContent[nlx - 1][k] = MaxWaterCapacity[nlx - 1];
            VolNo3NContent[nlx - 1][k] = nitconc * MaxWaterCapacity[nlx - 1];
            VolUreaNContent[nlx - 1][k] = nurconc * MaxWaterCapacity[nlx - 1];
            SoilNitrogenLoss += saven - (VolNo3NContent[nlx - 1][k] +
                                         VolUreaNContent[nlx - 1][k]) *
                                            dl[nlx - 1] * wk[k];
        }  // end if Vol...
    }      // end loop k
    return Drainage;
}
