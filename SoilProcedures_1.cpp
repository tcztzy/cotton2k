// File SoilProcedures_1.cpp
//
//   functions in this file:
// GetTargetStress()
// PredictDripIrrigation()
// PredictSurfaceIrrigation()
// AveragePsi()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
///////////////////////////////////////////////////////////////////////
double GetTargetStress()
//     This function computes and returns the target water stress factor.
//     A target water stress factor is defined for each growth stage.
//  The following growth stages are used: before first square; before
//  first bloom; first bloom to bloom + 20 days; bloom + 20 to bloom + 40;
//  bloom + 40 to first open boll; first open boll to 20% opening;
//  20% to 40% boll opening; 40% to 60% boll opening; 60% to 80% boll
//  opening; 80% to 90% boll opening.
//      Irrigation will be stopped after 90% boll opening.
//      Target stress parameters:
//    before 1st square:           strestgt[0]
//    before 1st bloom :           strestgt[1]
//    up to bloom + 20 days:       strestgt[2]
//    up to bloom + 40 days:       strestgt[3]
//    up to 1st open boll:         strestgt[4]
//    up to 20% open bolls:        strestgt[5]
//    up to 40% open bolls:        strestgt[6]
//    up to 60% open bolls:        strestgt[7]
//    up to 80% open bolls:        strestgt[8]
//    up to 90% open bolls:        strestgt[9]
//
//     This function is called from ComputeIrrigation().
//     The following global variables are referenced here:
//       Daynum, FirstBloom, FirstSquare, Kday, NumGreenBolls, NumOpenBolls.
//     The following global variable is set here:  DayStopPredIrrig,.
{
    const double strestgt[10] = {.70, .95, .99, .99, .99,
                                 .95, .90, .80, .60, .40};
    //    The target stress is defined according the current stage of plant
    //  development. Irrigation is stopped when 90% of the bolls are open,
    //  or if the target stress is zero.
    double tgtstr;  // the target water stress at this growth stage
    if (Kday > 0 && FirstSquare <= 0)
        tgtstr = strestgt[0];
    else if (FirstBloom <= 0)
        tgtstr = strestgt[1];
    else if (Daynum <= (FirstBloom + 20))
        tgtstr = strestgt[2];
    else if (Daynum <= (FirstBloom + 40))
        tgtstr = strestgt[3];
    else if (NumOpenBolls <= 0.01)
        tgtstr = strestgt[4];
    else if (NumOpenBolls < (0.25 * NumGreenBolls))
        tgtstr = strestgt[5];
    else if (NumOpenBolls < (0.667 * NumGreenBolls))
        tgtstr = strestgt[6];
    else if (NumOpenBolls < (1.5 * NumGreenBolls))
        tgtstr = strestgt[7];
    else if (NumOpenBolls < (4.0 * NumGreenBolls))
        tgtstr = strestgt[8];
    else if (NumOpenBolls < (9.0 * NumGreenBolls))
        tgtstr = strestgt[9];
    else {
        DayStopPredIrrig = Daynum;
        tgtstr = -9999;
    }
    if (tgtstr <= 0) {
        DayStopPredIrrig = Daynum;
        tgtstr = -9999;
    }
    return tgtstr;
}
//////////////////////////////////////////////////////////////////////////////////////
void PredictDripIrrigation(double TargetStress)
//     This function computes the amount of water (mm) needed for predicted drip
//  irrigation, considering the effects of water stress.
//     It is called from ComputeIrrigation(). It calls the function
//     GetFromClim(). The following global variables are referenced here:
//       ActualSoilEvaporation, ActualTranspiration, Daynum, DayStartPredIrrig,
//       Irrig, LastIrrigation, MaxIrrigation, MinDaysBetweenIrrig,
//       NumIrrigations, WaterStress.
//     The following global variable is set here:    AppliedWater.
//     The argument used is TargetStress.
//
{
    static bool irr1st;  // switch is set to true after first drip irrigation
    if (Daynum <= DayStartPredIrrig)
        irr1st = false;  // set to false on first day of irrigation prediction
                         //
    static double RequiredWater;  // the amount of water (mm) required for this
                                  // or next irrigation
    //     The amount of the first drip irrigation is set as 30 mm, or
    //     MaxIrrigation
    if (!irr1st) {
        //     Check if there is an irrigation defined by input, or rain, on
        //     this day.
        bool bIsIrr = false;
        for (int j = 0; j < NumIrrigations; j++) {
            if (Irrig[j].day == Daynum ||
                GetFromClim(RAIN, Daynum) > 1) {
                bIsIrr = true;  // there is an irrigation today
                break;
            }
        }
        //     The first predicted drip irrigation is applied if there is no
        //     scheduled irrigation
        //  or rain today, and water stress is not higher than 0.99
        if (!bIsIrr && WaterStress <= 0.99) {
            //     The first predicted drip irrigation is applied
            AppliedWater = fmin(30., MaxIrrigation);
            irr1st = true;
            RequiredWater = 0;
        }
        return;
    }
    //     The following is executed after the first drip irrigation has been
    //     applied. The "Required water" is computed by adding the amount needed
    //     to replace the water loss
    //  from the soil by evapotranspiration today.
    RequiredWater += ActualTranspiration + ActualSoilEvaporation -
                     GetFromClim(RAIN, Daynum);
    if (RequiredWater < 0) RequiredWater = 0;
    if ((Daynum - MinDaysBetweenIrrig) >= LastIrrigation) {
        //     If the minimum number of days (MinDaysBetweenIrrig) have passed
        //     after the last irrigation,
        //  A drip irrigation may be applied today.
        double facirr;  // factor used to modify irrigation amount.
        //     Compute the irrigation factor (facirr) from the ratio between the
        //     target water stress
        //  and the actual water stress. facirr can not be less than 0.80 or
        //  more than 1.25.
        //     Amount of irrigation is modified by the ratio of target stress to
        //     actual water stress.
        if (TargetStress > WaterStress)
            facirr = 1.20 * TargetStress / WaterStress;  // increase water
        else
            facirr = 0.90 * TargetStress / WaterStress;  // decrease water
        if (facirr < 0.80) facirr = 0.80;
        if (facirr > 1.25) facirr = 1.25;
        //     The amount of water to be applied (RequiredWater), computed from
        //     cumulative water
        //  loss since the last irrigation, is multiplied by facirr.
        if ((RequiredWater * facirr) > MaxIrrigation) {
            //     It is checked if it is greater than MaxIrrigation. in this
            //     case only the maximum possible
            //  amount is applied, and the difference will be added in the next
            //  irrigation.
            AppliedWater = MaxIrrigation;
            RequiredWater -= MaxIrrigation;
        } else {
            //     All the required water, modified by facirr, is applied.
            //     "Required Water" is set to zero.
            AppliedWater = RequiredWater * facirr;
            RequiredWater = 0;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
void PredictSurfaceIrrigation(double TargetStress)
//     This function computes the amount of predicted irrigation applied at soil
//     surface,
//  by sprinkler or furrow. It is called from ComputeIrrigation().
//
//     The following global variables are referenced here:
//       Daynum, DayStartPredIrrig, dl, IrrigationDepth, LastIrrigation,
//       MaxWaterCapacity, MinDaysBetweenIrrig, nl, RowSpace, VolWaterContent,
//       WaterStress, wk
//     The following global variable is set here:       AppliedWater.
//     The argument used: TargetStress
{
    //     Prediction of surface irrigation: loop over soil layers and check
    //  if you have reached the target depth of irrigation (IrrigationDepth).
    static int
        nDaysBelowTargetStress;  // number of days, since the last irrigation,
                                 // with water stress below the target stress.
    static int nIrrLayers = 0;   // number of soil layers to be irrigated.
    if (Daynum <= DayStartPredIrrig) {
        nDaysBelowTargetStress = 0;
        double sumdl =
            0;  // sum of thickness of all soil layers to be irrigated.
        for (int l = 0; l < nl; l++) {
            sumdl += dl[l];
            if (sumdl > IrrigationDepth) {
                nIrrLayers = l;
                break;
            }
        }
    }
    //     If the minimum number of days (MinDaysBetweenIrrig) minus 2 days have
    //     passed after
    //  the last irrigation, sum up days with water stress below target stress.
    if ((Daynum - MinDaysBetweenIrrig) >= (LastIrrigation - 2)) {
        //     Irrigation will be applied when water stress is less than the
        //     target stress for
        //  three days since the last irrigation.
        if (Daynum > DayStartPredIrrig && WaterStress < TargetStress) {
            nDaysBelowTargetStress++;
            if (nDaysBelowTargetStress >= 3) {
                double RequiredWater =
                    0;  // the amount of water required for irrigation
                for (int l = 0; l < nIrrLayers; l++)
                    for (int k = 0; k < nk; k++) {
                        //  as the deficit to MaxWaterCapacity down to
                        //  this depth. RequiredWater is converted from cm3 per
                        //  slab to mm.
                        double defcit;  // water content deficit to irrigation
                                        // depth
                        defcit =
                            MaxWaterCapacity[l] -
                            VolWaterContent[l][k];  // water content deficit
                        RequiredWater += dl[l] * wk[k] * defcit;
                    }
                AppliedWater = RequiredWater * 10 / RowSpace;
                //     The amount of water to be applied is checked not to
                //     exceed MaxIrrigation.
                if (AppliedWater > MaxIrrigation) AppliedWater = MaxIrrigation;
                nDaysBelowTargetStress = 0;
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
double AveragePsi()
//     This function computes and returns the average soil water potential of
//     the root zone of
//  the soil slab. This average is weighted by the amount of active roots (roots
//  capable of uptake) in each soil cell. Soil zones without roots are not
//  included.
//     It is called from SoilProcedures().
//     The functions psiq() and PsiOsmotic() are called.
//
//     The following global variables are referenced here:
//       airdr, alpha, beta, dl, ElCondSatSoilToday, NumLayersWithRoots,
//       RootColNumLeft, RootColNumRight, RootWtCapblUptake, SoilHorizonNum,
//       thetas, VolWaterContent, wk.
{
    //     Constants used:
    const double vrcumin = 0.1e-9;
    const double vrcumax = 0.025;
    //
    double psinum[9];  // sum of weighting coefficients for computing avgwat.
    double
        sumwat[9];  // sum of weighted soil water content for computing avgwat.
    double sumdl[9];  // sum of thickness of all soil layers containing roots.
                      //     Assign zero to some variables used for summation.
    for (int j = 0; j < 9; j++) {
        psinum[j] = 0;
        sumwat[j] = 0;
        sumdl[j] = 0;
    }
    //     Compute sum of dl as sumdl for each soil horizon.
    for (int l = 0; l < NumLayersWithRoots; l++) {
        int j = SoilHorizonNum[l];
        sumdl[j] += dl[l];
        for (int k = RootColNumLeft[l]; k <= RootColNumRight[l]; k++) {
            //     Check that RootWtCapblUptake in any cell is more than a
            //     minimum value vrcumin.
            if (RootWtCapblUptake[l][k] >= vrcumin) {
                //     Compute sumwat as the weighted sum of the water content,
                //     and psinum as the sum of
                //  these weights. Weighting is by root weight capable of
                //  uptake, or if it exceeds a maximum value (vrcumax) this
                //  maximum value is used for weighting.
                sumwat[j] += VolWaterContent[l][k] * dl[l] * wk[k] *
                             fmin(RootWtCapblUptake[l][k], vrcumax);
                psinum[j] +=
                    dl[l] * wk[k] * fmin(RootWtCapblUptake[l][k], vrcumax);
            }
        }
    }
    double sumpsi = 0;  // weighted sum of avgpsi
    double sumnum =
        0;  // sum of weighting coefficients for computing AverageSoilPsi.
    for (int j = 0; j < 9; j++) {
        if (psinum[j] > 0 && sumdl[j] > 0) {
            double avgwat;  // weighted average soil water content (V/V) in root
                            // zone
            //     Compute avgwat and the parameters to compute the soil water
            //     potential in each soil horizon
            avgwat = sumwat[j] / psinum[j];
            //     Soil water potential computed for a soil profile layer:
            double avgpsi =
                psiq(avgwat, airdr[j], thetas[j], alpha[j], beta[j]) -
                PsiOsmotic(avgwat, thetas[j], ElCondSatSoilToday);
            //     Use this to compute the average for the whole root zone.
            sumpsi += avgpsi * psinum[j];
            sumnum += psinum[j];
        }
    }
    double averagePsi;  // average soil water potential for the whole root zone
    if (sumnum > 0)
        averagePsi = sumpsi / sumnum;
    else
        averagePsi = 0;
    return averagePsi;
}
