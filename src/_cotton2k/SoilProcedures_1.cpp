// File SoilProcedures_1.cpp
//
//   functions in this file:
// SoilProcedures()
// RootsCapableOfUptake()
// ApplyFertilizer()
// GetTargetStress()
// AveragePsi()
// WaterTable()
//
#include <cmath>
#include "global.h"
#include "GeneralFunctions.h"
#include "Simulation.hpp"

void ComputeIrrigation(Simulation &, uint32_t);

double GetTargetStress(unsigned int, unsigned int, unsigned int, unsigned int, double, double);

/////////////////////////////////////////////////////////////////////////////////////////
void RootsCapableOfUptake(Soil &soil)
//     This function computes the weight of roots capable of uptake for all soil cells.
//  It is called from SoilProcedures().
//
//     The following global variables are referenced here:
//       nk, nl, RootColNumLeft, RootColNumRight, RootWeight.
//     The following global variable is set here:     RootWtCapblUptake
{
    const double cuind[3] = {1, 0.5, 0}; // the indices for the relative capability of uptake
    // (between 0 and 1) of water and nutrients by root age classes.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++)
            soil.cells[l][k].root.weight_capable_uptake = 0;
    //     Loop for all soil soil cells with roots. compute for each soil cell root-weight capable
    //  of uptake (RootWtCapblUptake) as the sum of products of root weight and capability of
    //  uptake index (cuind) for each root class in it.
    for (int l = 0; l < soil.number_of_layers_with_root; l++)
        for (int k = soil.layers[l].number_of_left_columns_with_root; k <= soil.layers[l].number_of_right_columns_with_root; k++)
            for (int i = 0; i < 3; i++)
                if (soil.cells[l][k].root.weight[i] > 1.e-15)
                    soil.cells[l][k].root.weight_capable_uptake += soil.cells[l][k].root.weight[i] * cuind[i];
}

/////////////////////////////////////////////////////////////////////////////////////////
void ApplyFertilizer(Simulation &sim, unsigned int u)
//     This function simulates the application of nitrogen fertilizer on each date
//  of application. It is called from SoilProcedures().
//
//     The following global variables are referenced here:
//       LocationColumnDrip, LocationLayerDrip,
//       NFertilizer, nk, nl, NumNitApps.
//     The following global variables are set here:
//       CumFertilizerN, VolNh4NContent, VolNo3NContent, VolUreaNContent.
{
    State &state = sim.states[u];
    const double ferc = 0.01; //  constant used to convert kgs per ha to mg cm-2
                              //     Loop over all fertilizer applications.
    for (int i = 0; i < NumNitApps; i++)
    {
        //     Check if fertilizer is to be applied on this day.
        if (sim.day_start + u == NFertilizer[i].day)
        {
            //     Compute the sum of mineral N applied as fertilizer (CumFertilizerN, in mg per slab).
            CumFertilizerN += ferc * sim.row_space * (NFertilizer[i].amtamm + NFertilizer[i].amtnit + NFertilizer[i].amtura);
            //     If this is a BROADCAST fertilizer application:
            if (NFertilizer[i].mthfrt == 0)
            {
                //     Compute the number of layers affected by broadcast fertilizer incorporation (lplow),
                //  assuming that the depth of incorporation is 20 cm.
                int lplow = 0;  // number of soil layers affected by cultivation
                double sdl = 0; // sum of depth of consecutive soil layers
                for (int l = 0; l < nl; l++)
                {
                    sdl += dl(l);
                    if (sdl >= 20)
                    {
                        lplow = l + 1;
                        break;
                    }
                }
                //     Calculate the actual depth of fertilizer incorporation in the soil (fertdp) as the sum of
                //  all soil layers affected by incorporation.
                double fertdp = 0; // depth of broadcast fertilizer incorporation, cm
                for (int l = 0; l < lplow; l++)
                    fertdp += dl(l);
                //     Update the nitrogen contents of all soil soil cells affected by this fertilizer application.
                for (int l = 0; l < lplow; l++)
                    for (int k = 0; k < nk; k++)
                    {
                        VolNh4NContent[l][k] += NFertilizer[i].amtamm * ferc / fertdp;
                        state.soil.cells[l][k].nitrate_nitrogen_content += NFertilizer[i].amtnit * ferc / fertdp;
                        VolUreaNContent[l][k] += NFertilizer[i].amtura * ferc / fertdp;
                    }
            }
            //     If this is a FOLIAR fertilizer application:
            else if (NFertilizer[i].mthfrt == 2)
            {
                //     It is assumed that 70% of the amount of ammonium or urea intercepted by the canopy
                //  is added to the leaf N content (state.leaf_nitrogen).
                state.leaf_nitrogen += 0.70 * state.light_interception * (NFertilizer[i].amtamm + NFertilizer[i].amtura) * 1000 / sim.plant_population;
                //     The amount not intercepted by the canopy is added to the soil. If the fertilizer is
                //  nitrate, it is assumed that all of it is added to the upper soil layer.
                //     Update nitrogen contents of the upper layer.
                for (int k = 0; k < nk; k++)
                {
                    VolNh4NContent[0][k] += NFertilizer[i].amtamm * (1 - 0.70 * state.light_interception) * ferc / dl(0);
                    state.soil.cells[0][k].nitrate_nitrogen_content += NFertilizer[i].amtnit * ferc / dl(0);
                    VolUreaNContent[0][k] += NFertilizer[i].amtura * (1 - 0.70 * state.light_interception) * ferc / dl(0);
                }
            }
            //     If this is a SIDE-DRESSING of N fertilizer:
            else if (NFertilizer[i].mthfrt == 1)
            {
                //     Define the soil column (ksdr) and the soil layer (lsdr) in
                //  which the side-dressed fertilizer is applied.
                int ksdr = NFertilizer[i].ksdr; // the column in which the side-dressed is applied
                int lsdr = NFertilizer[i].ksdr; // the layer in which the side-dressed is applied
                int n00 = 1;                    // number of soil soil cells in which side-dressed fertilizer is incorporated.
                                                //     If the volume of this soil cell is less than 100 cm3, it is assumed that the fertilizer
                                                //  is also incorporated in the soil cells below and to the sides of it.
                if ((dl(lsdr) * wk(ksdr, sim.row_space)) < 100)
                {
                    if (ksdr < nk - 1)
                        n00++;
                    if (ksdr > 0)
                        n00++;
                    if (lsdr < nl - 1)
                        n00++;
                }
                double addamm; // amount of ammonium N added to the soil by sidedressing (mg per cell)
                double addnit; // amount of nitrate N added to the soil by sidedressing (mg per cell)
                double addnur; // amount of urea N added to the soil by sidedressing (mg per cell)
                addamm = NFertilizer[i].amtamm * ferc * sim.row_space / n00;
                addnit = NFertilizer[i].amtnit * ferc * sim.row_space / n00;
                addnur = NFertilizer[i].amtura * ferc * sim.row_space / n00;
                //     Update the nitrogen contents of these soil cells.
                state.soil.cells[lsdr][ksdr].nitrate_nitrogen_content += addnit / (dl(lsdr) * wk(ksdr, sim.row_space));
                VolNh4NContent[lsdr][ksdr] += addamm / (dl(lsdr) * wk(ksdr, sim.row_space));
                VolUreaNContent[lsdr][ksdr] += addnur / (dl(lsdr) * wk(ksdr, sim.row_space));
                if ((dl(lsdr) * wk(ksdr, sim.row_space)) < 100)
                {
                    if (ksdr < nk - 1)
                    {
                        int kp1 = ksdr + 1; // column to the right of ksdr.
                        state.soil.cells[lsdr][kp1].nitrate_nitrogen_content += addnit / (dl(lsdr) * wk(kp1, sim.row_space));
                        VolNh4NContent[lsdr][kp1] += addamm / (dl(lsdr) * wk(kp1, sim.row_space));
                        VolUreaNContent[lsdr][kp1] += addnur / (dl(lsdr) * wk(kp1, sim.row_space));
                    }
                    if (ksdr > 0)
                    {
                        int km1 = ksdr - 1; // column to the left of ksdr.
                        state.soil.cells[lsdr][km1].nitrate_nitrogen_content += addnit / (dl(lsdr) * wk(km1, sim.row_space));
                        VolNh4NContent[lsdr][km1] += addamm / (dl(lsdr) * wk(km1, sim.row_space));
                        VolUreaNContent[lsdr][km1] += addnur / (dl(lsdr) * wk(km1, sim.row_space));
                    }
                    if (lsdr < nl - 1)
                    {
                        int lp1 = lsdr + 1;
                        state.soil.cells[lp1][ksdr].nitrate_nitrogen_content += addnit / (dl(lp1) * wk(ksdr, sim.row_space));
                        VolNh4NContent[lp1][ksdr] += addamm / (dl(lp1) * wk(ksdr, sim.row_space));
                        VolUreaNContent[lp1][ksdr] += addnur / (dl(lp1) * wk(ksdr, sim.row_space));
                    }
                }
            }
            //     If this is FERTIGATION (N fertilizer applied in drip irrigation):
            else if (NFertilizer[i].mthfrt == 3)
            {
                //      Convert amounts added to mg cm-3, and update the nitrogen content of the
                //  soil cell in which the drip outlet is situated.
                VolNh4NContent[LocationLayerDrip][LocationColumnDrip] +=
                    NFertilizer[i].amtamm * ferc * sim.row_space / (dl(LocationLayerDrip) * wk(LocationColumnDrip, sim.row_space));
                state.soil.cells[LocationLayerDrip][LocationColumnDrip].nitrate_nitrogen_content +=
                    NFertilizer[i].amtnit * ferc * sim.row_space / (dl(LocationLayerDrip) * wk(LocationColumnDrip, sim.row_space));
                VolUreaNContent[LocationLayerDrip][LocationColumnDrip] +=
                    NFertilizer[i].amtura * ferc * sim.row_space / (dl(LocationLayerDrip) * wk(LocationColumnDrip, sim.row_space));
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////
double GetTargetStress(unsigned int Daynum, unsigned int Kday, unsigned int FirstBloom, unsigned int FirstSquare, double NumGreenBolls, double NumOpenBolls)
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
//       FirstBloom, FirstSquare, Kday, NumGreenBolls, NumOpenBolls.
//     The following global variable is set here:  DayStopPredIrrig,.
{
    const double strestgt[10] = {.70, .95, .99, .99, .99, .95, .90, .80, .60, .40};
    //    The target stress is defined according the current stage of plant
    //  development. Irrigation is stopped when 90% of the bolls are open,
    //  or if the target stress is zero.
    double tgtstr; // the target water stress at this growth stage
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
    else
    {
        DayStopPredIrrig = Daynum;
        tgtstr = -9999;
    }
    if (tgtstr <= 0)
    {
        DayStopPredIrrig = Daynum;
        tgtstr = -9999;
    }
    return tgtstr;
}

////////////////////////////////////////////////////////////////////////////////
double AveragePsi(const State &state, double row_space)
//     This function computes and returns the average soil water potential of the root zone of
//  the soil slab. This average is weighted by the amount of active roots (roots capable of
//  uptake) in each soil cell. Soil zones without roots are not included.
//     It is called from SoilProcedures().
//     The functions psiq() and PsiOsmotic() are called.
//
//     The following global variables are referenced here:
//       airdr, alpha, vanGenuchtenBeta, dl, ElCondSatSoilToday, NumLayersWithRoots, RootColNumLeft,
//       RootColNumRight, RootWtCapblUptake, SoilHorizonNum, thetas.
{
    //     Constants used:
    const double vrcumin = 0.1e-9;
    const double vrcumax = 0.025;
    //
    double psinum[9]; // sum of weighting coefficients for computing avgwat.
    double sumwat[9]; // sum of weighted soil water content for computing avgwat.
    double sumdl[9];  // sum of thickness of all soil layers containing roots.
                      //     Assign zero to some variables used for summation.
    for (int j = 0; j < 9; j++)
    {
        psinum[j] = 0;
        sumwat[j] = 0;
        sumdl[j] = 0;
    }
    //     Compute sum of dl as sumdl for each soil horizon.
    for (int l = 0; l < state.soil.number_of_layers_with_root; l++)
    {
        int j = SoilHorizonNum[l];
        sumdl[j] += dl(l);
        for (int k = state.soil.layers[l].number_of_left_columns_with_root; k <= state.soil.layers[l].number_of_right_columns_with_root; k++)
        {
            //     Check that RootWtCapblUptake in any cell is more than a minimum value vrcumin.
            if (state.soil.cells[l][k].root.weight_capable_uptake >= vrcumin)
            {
                //     Compute sumwat as the weighted sum of the water content, and psinum as the sum of
                //  these weights. Weighting is by root weight capable of uptake, or if it exceeds a maximum
                //  value (vrcumax) this maximum value is used for weighting.
                sumwat[j] += state.soil.cells[l][k].water_content * dl(l) * wk(k, row_space) * std::min(state.soil.cells[l][k].root.weight_capable_uptake, vrcumax);
                psinum[j] += dl(l) * wk(k, row_space) * std::min(state.soil.cells[l][k].root.weight_capable_uptake, vrcumax);
            }
        }
    }
    double sumpsi = 0; // weighted sum of avgpsi
    double sumnum = 0; // sum of weighting coefficients for computing AverageSoilPsi.
    for (int j = 0; j < 9; j++)
    {
        if (psinum[j] > 0 && sumdl[j] > 0)
        {
            double avgwat; // weighted average soil water content (V/V) in root zone
                           //     Compute avgwat and the parameters to compute the soil water potential in each soil horizon
            avgwat = sumwat[j] / psinum[j];
            //     Soil water potential computed for a soil profile layer:
            double avgpsi = psiq(avgwat, airdr[j], thetas[j], alpha[j], vanGenuchtenBeta[j]) - PsiOsmotic(avgwat, thetas[j], ElCondSatSoilToday);
            //     Use this to compute the average for the whole root zone.
            sumpsi += avgpsi * psinum[j];
            sumnum += psinum[j];
        }
    }
    double averagePsi; // average soil water potential for the whole root zone
    if (sumnum > 0)
        averagePsi = sumpsi / sumnum;
    else
        averagePsi = 0;
    return averagePsi;
}

/////////////////////////////////////////////////////////////////////////////
void WaterTable(Simulation &sim, unsigned int u)
//     This function sets the water saturation of the soil layers below the water
//  table, if it has been defined in the input. It is called from  SoilProcedures()
//  if water table data have been input.
//
//     The following global variables are referenced here:
//       dl, DayWaterTableInput, ElCondSatSoil, LevelsOfWaterTable, MaxWaterCapacity,
//       nk, nl, NumWaterTableData, PoreSpace, MaxWaterCapacity, wk.
//
//     The following global variables are set here:
//       ElCondSatSoilToday, WaterTableLayer.
//
{
    if (NumWaterTableData <= 0)
        return;
    //     Find the depth of water table for this day.
    double lwtable = 201; // level of water table on this day, cm
    for (int i = 0; i < NumWaterTableData; i++)
    {
        if (sim.day_start + u >= DayWaterTableInput[i])
        {
            lwtable = LevelsOfWaterTable[i];
            ElCondSatSoilToday = ElCondSatSoil[i];
        }
    }
    //     Find the number of the uppermost layer of water table
    if (lwtable > 200)
        WaterTableLayer = 1000;
    else
    {
        double sumdl = 0; // sum of depth of consecutive soil layers
        for (int l = 0; l < nl; l++)
        {
            sumdl += dl(l);
            if (sumdl > lwtable)
            {
                WaterTableLayer = l;
                break;
            }
        }
    }
    // check the water balance in the soil.
    for (int l = 0; l < nl; l++)
        if (l >= WaterTableLayer)
            for (int k = 0; k < nk; k++)
                sim.states[u].soil.cells[l][k].water_content = PoreSpace[l];
        else
            //     Make sure that (in case water table was lowered) water content is not
            //  higher than MaxWaterCapacity.
            for (int k = 0; k < nk; k++)
                if (sim.states[u].soil.cells[l][k].water_content > MaxWaterCapacity[l])
                    sim.states[u].soil.cells[l][k].water_content = MaxWaterCapacity[l];
}
