// File SoilProcedures_1.cpp
//
//   functions in this file:
// SoilProcedures()
// RootsCapableOfUptake()
// ApplyFertilizer()
// ComputeIrrigation()
// GetTargetStress()
// PredictDripIrrigation()
// PredictSurfaceIrrigation()
// OutputPredictedIrrigation()
// AveragePsi()
// WaterTable()
//
#include "global.h"
#include "GeneralFunctions.h"

void RootsCapableOfUptake(Soil &);

void ApplyFertilizer(Simulation &, unsigned int);

void ComputeIrrigation(Simulation &, uint32_t);

double GetTargetStress(const int &, const int &, const int &, double, double);

void PredictDripIrrigation(Simulation &, uint32_t, double, const double &);

void PredictSurfaceIrrigation(Simulation &, unsigned int, double);

void OutputPredictedIrrigation(double, double, const string &, const int &, const int &, const double &);

double AveragePsi(const State &, double);

void WaterTable(Simulation &, unsigned int); // WATERTBL
// SoilProcedure_2
void CapillaryFlow(Simulation &, unsigned int);

void DripFlow(SoilCell[40][20], double, double);

// SoilProcedures_3
void WaterUptake(Simulation &, unsigned int); // UPTAKE
void GravityFlow(SoilCell[40][20], double, double);

//////////////////////////
void SoilProcedures(Simulation &sim, uint32_t u)
//     This function manages all the soil related processes, and is executed once each
//  day. It is called from SimulateThisDay() and it calls the following functions:
//  ApplyFertilizer(), AveragePsi(), CapillaryFlow(), ComputeIrrigation(), DripFlow(),
//  GravityFlow(), RootsCapableOfUptake(), WaterUptake(), WaterTable()
//
//     The following global variables are referenced here:
//       ActualTranspiration, Clim, DayStartPredIrrig, DayStopPredIrrig,
//       dl, Irrig, IrrigMethod, isw, Kday, MaxIrrigation, nk, nl, NumIrrigations,
//       NumWaterTableData, OutIndex, PerPlantArea, SoilPsi, SupplyNH4N, SupplyNO3N,
//       VolWaterContent, wk.
//     The following global variables are set here:
//       AverageSoilPsi, CumNitrogenUptake, CumTranspiration, CumWaterAdded, LocationColumnDrip,
//       LocationLayerDrip, noitr.
{
    State &state = sim.states[u];
    //     The following constant parameters are used:
    const double cpardrip = 0.2;
    const double cparelse = 0.4;
    //     Call function ApplyFertilizer() for nitrogen fertilizer application.
    ApplyFertilizer(sim, u);
    double DripWaterAmount = 0; // amount of water applied by drip irrigation
    double WaterToApply;        // amount of water applied by non-drip irrigation or rainfall
                                //     Check if there is rain on this day
    WaterToApply = sim.climate[u].Rain;
    //     If irrigation is to be predicted for this day, call ComputeIrrigation()
    //  to compute the actual amount of irrigation.
    if (MaxIrrigation > 0)
        if (state.daynum >= DayStartPredIrrig && state.daynum < DayStopPredIrrig)
        {
            ComputeIrrigation(sim, u);
            if (IrrigMethod == 2)
                DripWaterAmount = state.applied_water;
            else
                WaterToApply += state.applied_water;
            state.applied_water = 0;
        }
    //     When water is added by an irrigation defined in the input: update the amount
    //  of applied water.
    for (int i = 0; i < NumIrrigations; i++)
    {
        if (state.daynum == sim.irrigation[i].day)
        {
            if (sim.irrigation[i].method == 2)
            {
                DripWaterAmount += sim.irrigation[i].amount;
                LocationColumnDrip = sim.irrigation[i].LocationColumnDrip;
                LocationLayerDrip = sim.irrigation[i].LocationLayerDrip;
            }
            else
                WaterToApply += sim.irrigation[i].amount;
            break;
        }
    }
    CumWaterAdded += WaterToApply + DripWaterAmount;
    if (Kday > 0)
        Scratch21[u].amitri = WaterToApply + DripWaterAmount;
    //     The following will be executed only after plant emergence
    if (state.daynum >= sim.day_emerge && isw > 0)
    {
        RootsCapableOfUptake(state.soil);   // function computes roots capable of uptake for each soil cell
        AverageSoilPsi = AveragePsi(state, sim.row_space); // function computes the average matric soil water
                                            //                      potential in the root zone, weighted by the roots-capable-of-uptake.
        WaterUptake(sim, u);                // function  computes water and nitrogen uptake by plants.
                                            //     Update the cumulative sums of actual transpiration (CumTranspiration, mm) and total uptake
                                            //  of nitrogen (CumNitrogenUptake, mg N per slab, converted from total N supply, g per plant).
        state.cumulative_transpiration += state.actual_transpiration;
        CumNitrogenUptake += (SupplyNO3N + SupplyNH4N) * 10 * sim.row_space / PerPlantArea;
    }
    //     Call function WaterTable() for saturating soil below water table.
    if (NumWaterTableData > 0)
        WaterTable(sim, u);
    if (WaterToApply > 0)
    {
        //     For rain or surface irrigation.
        //     The number of iterations is computed from the thickness of the first soil layer.
        noitr = (int)(cparelse * WaterToApply / (dl(0) + 2) + 1);
        double applywat; // the amount of water applied, mm per iteration.
        applywat = WaterToApply / noitr;
        //     The following subroutines are called noitr times per day:
        //  If water is applied, GravityFlow() is called when the method of irrigation is not by
        //  drippers, followed by CapillaryFlow().
        for (int iter = 0; iter < noitr; iter++)
        {
            GravityFlow(sim.states[u].soil.cells, applywat, sim.row_space);
            CapillaryFlow(sim, u);
        }
    }
    if (DripWaterAmount > 0)
    {
        //     For drip irrigation.
        //     The number of iterations is computed from the volume of the soil cell in which
        //  the water is applied.
        noitr = (int)(cpardrip * DripWaterAmount / (dl(LocationLayerDrip) * wk(LocationColumnDrip, sim.row_space)) + 1);
        double applywat; // the amount of water applied, mm per iteration.
        applywat = DripWaterAmount / noitr;
        // If water is applied, DripFlow() is called followed by CapillaryFlow().
        for (int iter = 0; iter < noitr; iter++)
        {
            DripFlow(state.soil.cells, applywat, sim.row_space);
            CapillaryFlow(sim, u);
        }
    }
    //     When no water is added, there is only one iteration in this day.
    if (WaterToApply + DripWaterAmount <= 0)
    {
        noitr = 1;
        CapillaryFlow(sim, u);
    }
    //     If the output flag OutIndex(17) is non-zero, write to output file *.WAT.
    //  This flag is also the interval in days between outputs. This is used for checking only.
    if (OutIndex[17] > 0)
    {
        int kkk = state.daynum / OutIndex[17];
        if (kkk * OutIndex[17] == state.daynum)
        {
            ofstream File42(fs::path("output") / (string(sim.profile_name) + ".WAT"), ios::app);
            File42 << endl
                   << " Average Values by Layers on Day of Year ";
            File42 << state.daynum << endl;
            File42 << " Layer          VolWaterContent          SoilPsi " << endl;
            for (int l = 0; l < nl; l++)
            {
                //     Compute for each layer the average water content and the average matric water potential.
                double psiavx = 0; // average of SoilPsi for a soil layer, bars.
                double vavg = 0;   // average water content in a soil layer, cm3 cm-3.
                for (int k = 0; k < nk; k++)
                {
                    vavg += VolWaterContent[l][k];
                    psiavx += SoilPsi[l][k];
                }
                vavg = vavg / nk;
                psiavx = psiavx / nk;
                File42.width(6);
                File42 << l + 1;
                File42.setf(ios::fixed);
                File42.precision(3);
                File42.width(15);
                File42 << vavg;
                File42.width(15);
                File42 << psiavx;
                File42 << endl;
            }
        }
    }
}

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
//       dl, LightIntercept, LocationColumnDrip, LocationLayerDrip,
//       NFertilizer, nk, nl, NumNitApps, wk.
//     The following global variables are set here:
//       CumFertilizerN, LeafNitrogen, VolNh4NContent, VolNo3NContent, VolUreaNContent.
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
                //  is added to the leaf N content (LeafNitrogen).
                LeafNitrogen += 0.70 * LightIntercept * (NFertilizer[i].amtamm + NFertilizer[i].amtura) * 1000 / PlantPopulation;
                //     The amount not intercepted by the canopy is added to the soil. If the fertilizer is
                //  nitrate, it is assumed that all of it is added to the upper soil layer.
                //     Update nitrogen contents of the upper layer.
                for (int k = 0; k < nk; k++)
                {
                    VolNh4NContent[0][k] += NFertilizer[i].amtamm * (1 - 0.70 * LightIntercept) * ferc / dl(0);
                    state.soil.cells[0][k].nitrate_nitrogen_content += NFertilizer[i].amtnit * ferc / dl(0);
                    VolUreaNContent[0][k] += NFertilizer[i].amtura * (1 - 0.70 * LightIntercept) * ferc / dl(0);
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

////////////////////////////////////////////////////////////////////////////////
void ComputeIrrigation(Simulation &sim, uint32_t u)
//     This function computes the amount of water (mm) applied by a predicted
//  irrigation. It is called from SoilProcedures().
//     It calls GetTargetStress(), PredictDripIrrigation(), PredictSurfaceIrrigation(),
//  OutputPredictedIrrigation().
//     The following global variables are referenced here:
//       AppliedWater, IrrigMethod.
//     The following global variable is set here:       LastIrrigation.
{
    double TargetStress = GetTargetStress(sim.day_start + u, sim.first_bloom, sim.first_square, sim.states[u].number_of_green_bolls, sim.states[u].number_of_open_bolls);
    if (TargetStress == -9999)
        return;
    //
    if (IrrigMethod == 2)
        PredictDripIrrigation(sim, u, TargetStress, sim.states[u].water_stress);
    else
        PredictSurfaceIrrigation(sim, u, TargetStress);
    //     If the amount of water to be applied (AppliedWater) is non zero update the date of
    //  last irrigation, and write report in output file *.B01.
    if (sim.states[u].applied_water > 0.00001)
    {
        LastIrrigation = sim.day_start + u;
        OutputPredictedIrrigation(sim.states[u].applied_water, TargetStress, sim.profile_name, sim.day_start + u, sim.year, sim.states[u].water_stress);
    }
}

///////////////////////////////////////////////////////////////////////
double GetTargetStress(const int &Daynum, const int &FirstBloom, const int &FirstSquare, double NumGreenBolls, double NumOpenBolls)
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

//////////////////////////////////////////////////////////////////////////////////////
void PredictDripIrrigation(Simulation &sim, uint32_t u, double TargetStress, const double &WaterStress)
//     This function computes the amount of water (mm) needed for predicted drip
//  irrigation, considering the effects of water stress.
//     It is called from ComputeIrrigation().
//     The following global variables are referenced here:
//       ActualSoilEvaporation, ActualTranspiration, DayStartPredIrrig, Irrig,
//       LastIrrigation, MaxIrrigation, MinDaysBetweenIrrig, NumIrrigations, WaterStress.
//     The following global variable is set here:    AppliedWater.
//     The argument used is TargetStress.
//
{
    State &state = sim.states[u];
    static bool irr1st; // switch is set to true after first drip irrigation
    if (sim.day_start + u <= DayStartPredIrrig)
        irr1st = false;          // set to false on first day of irrigation prediction
                                 //
    static double RequiredWater; // the amount of water (mm) required for this or next irrigation
                                 //     The amount of the first drip irrigation is set as 30 mm, or MaxIrrigation
    if (!irr1st)
    {
        //     Check if there is an irrigation defined by input, or rain, on this day.
        bool bIsIrr = false;
        for (int j = 0; j < NumIrrigations; j++)
        {
            if (sim.irrigation[j].day == sim.day_start + u || sim.climate[u].Rain > 1)
            {
                bIsIrr = true; // there is an irrigation today
                break;
            }
        }
        //     The first predicted drip irrigation is applied if there is no scheduled irrigation
        //  or rain today, and water stress is not higher than 0.99
        if (!bIsIrr && WaterStress <= 0.99)
        {
            //     The first predicted drip irrigation is applied
            state.applied_water = min(30., MaxIrrigation);
            irr1st = true;
            RequiredWater = 0;
        }
        return;
    }
    //     The following is executed after the first drip irrigation has been applied.
    //     The "Required water" is computed by adding the amount needed to replace the water loss
    //  from the soil by evapotranspiration today.
    RequiredWater += state.actual_transpiration + state.actual_soil_evaporation - sim.climate[u].Rain;
    if (RequiredWater < 0)
        RequiredWater = 0;
    if ((sim.day_start + u - MinDaysBetweenIrrig) >= LastIrrigation)
    {
        //     If the minimum number of days (MinDaysBetweenIrrig) have passed after the last irrigation,
        //  A drip irrigation may be applied today.
        double facirr; // factor used to modify irrigation amount.
                       //     Compute the irrigation factor (facirr) from the ratio between the target water stress
                       //  and the actual water stress. facirr can not be less than 0.80 or more than 1.25.
                       //     Amount of irrigation is modified by the ratio of target stress to actual water stress.
        if (TargetStress > WaterStress)
            facirr = 1.20 * TargetStress / WaterStress; // increase water
        else
            facirr = 0.90 * TargetStress / WaterStress; // decrease water
        if (facirr < 0.80)
            facirr = 0.80;
        if (facirr > 1.25)
            facirr = 1.25;
        //     The amount of water to be applied (RequiredWater), computed from cumulative water
        //  loss since the last irrigation, is multiplied by facirr.
        if ((RequiredWater * facirr) > MaxIrrigation)
        {
            //     It is checked if it is greater than MaxIrrigation. in this case only the maximum possible
            //  amount is applied, and the difference will be added in the next irrigation.
            state.applied_water = MaxIrrigation;
            RequiredWater -= MaxIrrigation;
        }
        else
        {
            //     All the required water, modified by facirr, is applied. "Required Water" is set to zero.
            state.applied_water = RequiredWater * facirr;
            RequiredWater = 0;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void PredictSurfaceIrrigation(Simulation &sim, unsigned int u, double TargetStress)
//     This function computes the amount of predicted irrigation applied at soil surface,
//  by sprinkler or furrow. It is called from ComputeIrrigation().
//
//     The following global variables are referenced here:
//       DayStartPredIrrig, dl, IrrigationDepth, LastIrrigation, MaxWaterCapacity,
//       MinDaysBetweenIrrig, nl, VolWaterContent, WaterStress, wk
//     The following global variable is set here:       AppliedWater.
//     The argument used: TargetStress
{
    State &state = sim.states[u];
    //     Prediction of surface irrigation: loop over soil layers and check
    //  if you have reached the target depth of irrigation (IrrigationDepth).
    static int nDaysBelowTargetStress; // number of days, since the last irrigation,
    // with water stress below the target stress.
    static int nIrrLayers = 0; // number of soil layers to be irrigated.
    if (sim.day_start + u <= DayStartPredIrrig)
    {
        nDaysBelowTargetStress = 0;
        double sumdl = 0; // sum of thickness of all soil layers to be irrigated.
        for (int l = 0; l < nl; l++)
        {
            sumdl += dl(l);
            if (sumdl > IrrigationDepth)
            {
                nIrrLayers = l;
                break;
            }
        }
    }
    //     If the minimum number of days (MinDaysBetweenIrrig) minus 2 days have passed after
    //  the last irrigation, sum up days with water stress below target stress.
    if ((sim.day_start + u - MinDaysBetweenIrrig) >= (LastIrrigation - 2))
    {
        //     Irrigation will be applied when water stress is less than the target stress for
        //  three days since the last irrigation.
        if (sim.day_start + u > DayStartPredIrrig && state.water_stress < TargetStress)
        {
            nDaysBelowTargetStress++;
            if (nDaysBelowTargetStress >= 3)
            {
                double RequiredWater = 0; // the amount of water required for irrigation
                for (int l = 0; l < nIrrLayers; l++)
                    for (int k = 0; k < nk; k++)
                    {
                        //  as the deficit to MaxWaterCapacity down to
                        //  this depth. RequiredWater is converted from cm3 per slab to mm.
                        double defcit;                                        // water content deficit to irrigation depth
                        defcit = MaxWaterCapacity[l] - VolWaterContent[l][k]; // water content deficit
                        RequiredWater += dl(l) * wk(k, sim.row_space) * defcit;
                    }
                state.applied_water = RequiredWater * 10 / sim.row_space;
                //     The amount of water to be applied is checked not to exceed MaxIrrigation.
                if (state.applied_water > MaxIrrigation)
                    state.applied_water = MaxIrrigation;
                nDaysBelowTargetStress = 0;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////
void OutputPredictedIrrigation(double AppliedWater, double TargetStress, const string &ProfileName, const int &Daynum, const int &year, const double &WaterStress)
//      This function is called from ComputeIrrigation().
//      It writes output ofapplication of predicted irrigation to file *.B01
//      Function DoyToDate() is used.
//      Arguments used: AppliedWater, TargetStress.
//      Global variables referenced: iyear, WaterStress.
{
    ofstream File20(fs::path("output") / (ProfileName + ".B01"), ios::app);
    File20 << " Predicted irrigation on " << DoyToDate(Daynum, year) << " - ";
    if (OutIndex[1] == 0)
        File20 << AppliedWater << " mm. ";
    else
        File20 << AppliedWater / 25.4 << " inches. ";
    File20 << " Water Stress: Actual = ";
    File20.setf(ios::fixed);
    File20.precision(3);
    File20.width(6);
    File20 << WaterStress;
    File20 << " Target = ";
    File20.width(6);
    File20 << TargetStress;
    File20 << endl;
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
//       RootColNumRight, RootWtCapblUptake, SoilHorizonNum, thetas, VolWaterContent, wk.
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
                sumwat[j] += VolWaterContent[l][k] * dl(l) * wk(k, row_space) * min(state.soil.cells[l][k].root.weight_capable_uptake, vrcumax);
                psinum[j] += dl(l) * wk(k, row_space) * min(state.soil.cells[l][k].root.weight_capable_uptake, vrcumax);
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
//       addwtbl, ElCondSatSoilToday, WaterTableLayer, VolWaterContent.
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
    //     The total water entering the soil slab (addwtbl) is computed. It is used to check
    //  the water balance in the soil.
    double vh2ocx; // previous water content of a cell
    for (int l = 0; l < nl; l++)
    {
        if (l >= WaterTableLayer)
        {
            for (int k = 0; k < nk; k++)
            {
                vh2ocx = VolWaterContent[l][k];
                VolWaterContent[l][k] = PoreSpace[l];
                addwtbl += 10 * (VolWaterContent[l][k] - vh2ocx) * dl(l) * wk(k, sim.row_space) / sim.row_space;
            }
        }
        else
        {
            //     Make sure that (in case water table was lowered) water content is not
            //  higher than MaxWaterCapacity and adjust addwtbl.
            for (int k = 0; k < nk; k++)
            {
                if (VolWaterContent[l][k] > MaxWaterCapacity[l])
                {
                    vh2ocx = VolWaterContent[l][k];
                    VolWaterContent[l][k] = MaxWaterCapacity[l];
                    addwtbl += 10 * (VolWaterContent[l][k] - vh2ocx) * dl(l) * wk(k, sim.row_space) / sim.row_space;
                }
            }
        }
    }
}
