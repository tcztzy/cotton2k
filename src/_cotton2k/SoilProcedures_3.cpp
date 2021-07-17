//  SoilProcedures_3.cpp
//
//   functions in this file:
// GravityFlow()
// WaterUptake()
// PsiOnTranspiration()
// NitrogenUptake()
// WaterFlux()
// WaterBalance()
// NitrogenFlow()
// SoilSum()
//
#include <cmath>
#include "global.h"
#include "GeneralFunctions.h"
#include "Simulation.hpp"

extern "C"{
    double PsiOnTranspiration(double);
}

void NitrogenUptake(State &, SoilCell &, int, int, double, double, double);

void WaterBalance(double[], double[], double [], int);

// SoilProcedures_2
double Drain(SoilCell[40][20], double);

////////////////////////////////////////////////////////////////////////////
void GravityFlow(SoilCell soil_cells[40][20], double applywat, double row_space)
//     This function computes the water redistribution in the soil or surface irrigation
//  (by flooding or sprinklers). It is called by SoilProcedures(). It calls function Drain().
//     The following argument is used:          ApplyWat = amount of water applied, mm.
//     The following global variables are referenced:
//       dl, nk.
//     The following global variables are set:
//       CumWaterDrained.
{
//     Add the applied amount of water to the top soil cell of each column.
    for (int k = 0; k < nk; k++)
        soil_cells[0][k].water_content += 0.10 * applywat / dl(0);
//     Call function Drain() to compute downflow of water.
    double WaterDrainedOut; // water drained out of the slab, mm.
    WaterDrainedOut = Drain(soil_cells, row_space);
//     If there is drainage out of the slab, transform it to mm,
//  and update the cumulative drainage (CumWaterDrained)
    if (WaterDrainedOut > 0)
        CumWaterDrained += 10 * WaterDrainedOut / row_space;
}

///////////////////////////////////////////////////////////////////////////////////
void WaterUptake(Simulation &sim, unsigned int u)
//     This function computes the uptake of water by plant roots from the soil
//  (i.e., actual transpiration rate). It is called from SoilProcedures().
//     It calls PsiOnTranspiration(), psiq(), PsiOsmotic().
//
//     The following global variables are referenced:
//  nk, nl, NumLayersWithRoots, ReferenceTransp, RootColNumLeft,
//  RootColNumRight, SoilHorizonNum, thetar, TotalRequiredN.
//     The following global variables are set:
//  ActualTranspiration, SoilPsi.
{
    State &state = sim.states[u];
    // Compute the modified light interception factor (LightInter1) for use in computing transpiration rate.
    double LightInter1; // modified light interception factor by canopy
    LightInter1 = state.light_interception * 1.55 - 0.32;
    if (LightInter1 < state.light_interception)
        LightInter1 = state.light_interception;
    if (LightInter1 > 1)
        LightInter1 = 1;

    // The potential transpiration is the product of the daytime Penman equation and LightInter1.
    double PotentialTranspiration = sim.states[u].evapotranspiration * LightInter1;
    Scratch21[u].ep = PotentialTranspiration;
    double upf[40][20];  // uptake factor, computed as a ratio, for each soil cell
    double uptk[40][20]; // actual transpiration from each soil cell, cm3 per day
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++) {
            upf[l][k] = 0;
            uptk[l][k] = 0;
        }
    double sumep = 0; // sum of actual transpiration from all soil soil cells, cm3 per day.

    // Compute the reduction due to soil moisture supply by function PsiOnTranspiration().
    double Transp; // the actual transpiration converted to cm3 per slab units.
    Transp = .10 * sim.row_space * PotentialTranspiration * PsiOnTranspiration(AverageSoilPsi);
    double difupt = 0; // the cumulative difference between computed transpiration and actual transpiration, in cm3, due to limitation of PWP.
    do {
        double supf = 0; // sum of upf for all soil cells
        for (int l = 0; l < sim.states[u].soil.number_of_layers_with_root; l++) {
            int j = SoilHorizonNum[l];
            // Compute, for each layer, the lower and upper water content limits for the transpiration function. These are set from limiting soil water potentials (-15 to -1 bars).
            double vh2lo; // lower limit of water content for the transpiration function
            double vh2hi; // upper limit of water content for the transpiration function
            vh2lo = qpsi(-15, thad[l], thts[l], alpha[j], vanGenuchtenBeta[j]);
            vh2hi = qpsi(-1, thad[l], thts[l], alpha[j], vanGenuchtenBeta[j]);
            for (int k = state.soil.layers[l].number_of_left_columns_with_root; k <= state.soil.layers[l].number_of_right_columns_with_root; k++) {
                double redfac; // reduction factor for water uptake, caused by low levels of soil
                // water, as a linear function of cell.water_content, between vh2lo and vh2hi.
                redfac = (state.soil.cells[l][k].water_content - vh2lo) / (vh2hi - vh2lo);
                if (redfac < 0)
                    redfac = 0;
                if (redfac > 1)
                    redfac = 1;
                // The computed 'uptake factor' (upf) for each soil cell is the product of 'root weight capable of uptake' and redfac.
                upf[l][k] = sim.states[u].soil.cells[l][k].root.weight_capable_uptake * redfac;
                supf += upf[l][k];
            }
        }

        difupt = 0;
        for (int l = 0; l < state.soil.number_of_layers_with_root; l++)
            for (int k = state.soil.layers[l].number_of_left_columns_with_root; k <= state.soil.layers[l].number_of_right_columns_with_root; k++)
                if (upf[l][k] > 0 && state.soil.cells[l][k].water_content > thetar[l]) {
                    // The amount of water extracted from each cell is proportional to its 'uptake factor'.
                    double upth2o;  // transpiration from a soil cell, cm3 per day
                    upth2o = Transp * upf[l][k] / supf;
                    // Update cell.water_content, storing its previous value as vh2ocx.
                    double vh2ocx; // previous value of water_content of this cell
                    vh2ocx = state.soil.cells[l][k].water_content;
                    state.soil.cells[l][k].water_content -= upth2o / (dl(l) * wk(k, sim.row_space));
                    // If the new value of cell.water_content is less than the permanent wilting point, modify the value of upth2o so that water_content will be equal to it.
                    if (state.soil.cells[l][k].water_content < thetar[l]) {
                        state.soil.cells[l][k].water_content = thetar[l];
                        double xupt;  // intermediate computation of upth2o

                        // Compute the difference due to this correction and add it to difupt.
                        xupt = (vh2ocx - thetar[l]) * dl(l) * wk(k, sim.row_space);
                        difupt += upth2o - xupt;
                        upth2o = xupt;
                    }
                    if (upth2o < 0)
                        upth2o = 0;

                    // Compute sumep as the sum of the actual amount of water extracted from all soil cells. Recalculate uptk of this soil cell as cumulative upth2o.
                    sumep += upth2o;
                    uptk[l][k] += upth2o;
                }

        // If difupt is greater than zero, redefine the variable Transp as difuptfor use in next loop.
        if (difupt > 0)
            Transp = difupt;
    } while (difupt > 0);

    // recompute SoilPsi for all soil cells with roots by calling function PSIQ,
    for (int l = 0; l < state.soil.number_of_layers_with_root; l++) {
        int j = SoilHorizonNum[l];
        for (int k = state.soil.layers[l].number_of_left_columns_with_root; k <= state.soil.layers[l].number_of_right_columns_with_root; k++)
            SoilPsi[l][k] = psiq(state.soil.cells[l][k].water_content, thad[l], thts[l], alpha[j], vanGenuchtenBeta[j])
                            - PsiOsmotic(state.soil.cells[l][k].water_content, thts[l], ElCondSatSoilToday);
    } // end l & k loops

    // compute ActualTranspiration as actual water transpired, in mm.
    state.actual_transpiration = sumep * 10 / sim.row_space;

    // Zeroize the amounts of NH4 and NO3 nitrogen taken up from the soil.
    state.supplied_nitrate_nitrogen = 0;
    state.supplied_ammonium_nitrogen = 0;

    // Compute the proportional N requirement from each soil cell with roots, and call function NitrogenUptake() to compute nitrogen uptake.
    if (sumep > 0 && state.total_required_nitrogen > 0) {
        for (int l = 0; l < state.soil.number_of_layers_with_root; l++)
            for (int k = state.soil.layers[l].number_of_left_columns_with_root; k <= state.soil.layers[l].number_of_right_columns_with_root; k++) {
                if (uptk[l][k] > 0) {
                    double reqnc; // proportional allocation of TotalRequiredN to each cell
                    reqnc = state.total_required_nitrogen * uptk[l][k] / sumep;
                    NitrogenUptake(state, state.soil.cells[l][k], l, k, reqnc, sim.row_space, sim.per_plant_area);
                } // end if uptk
            } // end loop k & l
    } // end if sumep
}

///////////////////////////////////////////////////////////////////////////
void NitrogenUptake(State &state, SoilCell &soil_cell, int l, int k, double reqnc, double row_space, double per_plant_area)
//     The function NitrogenUptake() computes the uptake of nitrate and ammonium N
//  from a soil cell. It is called by WaterUptake().
//     The arguments of this function are:
//        k, l - column and layer index of this cell.
//        reqnc - maximum N uptake (proportional to total N
//                required for plant growth), g N per plant.
//     The following global variables are set here:
//       VolNh4NContent, VolNo3NContent
{
//     Constant parameters:
    const double halfn = 0.08;// the N concentration in soil water (mg cm-3) at which
    // uptake is half of the possible rate.
    const double cparupmax = 0.5;  // constant parameter for computing upmax.
    const double p1 = 100, p2 = 5; // constant parameters for computing AmmonNDissolved.
//
    double coeff; // coefficient used to convert g per plant to mg cm-3 units.
    coeff = 10 * row_space / (per_plant_area * dl(l) * wk(k, row_space));
//     A Michaelis-Menten procedure is used to compute the rate of nitrate uptake from
//  each cell. The maximum possible amount of uptake is reqnc (g N per plant), and
//  the half of this rate occurs when the nitrate concentration in the soil solution is
//  halfn (mg N per cm3 of soil water).
//     Compute the uptake of nitrate from this soil cell, upno3c in g N per plant units.
//  Define the maximum possible uptake, upmax, as a fraction of VolNo3NContent.
    if (soil_cell.nitrate_nitrogen_content > 0) {
        double upno3c; // uptake rate of nitrate, g N per plant per day
        upno3c = reqnc * soil_cell.nitrate_nitrogen_content
                 / (halfn * soil_cell.water_content + soil_cell.nitrate_nitrogen_content);
        double upmax; // maximum possible uptake rate, mg N per soil cell per day
        upmax = cparupmax * soil_cell.nitrate_nitrogen_content;
//     Make sure that uptake will not exceed upmax and update VolNo3NContent and upno3c.
        if ((coeff * upno3c) < upmax)
            soil_cell.nitrate_nitrogen_content -= coeff * upno3c;
        else {
            soil_cell.nitrate_nitrogen_content -= upmax;
            upno3c = upmax / coeff;
        }
//     upno3c is added to the total uptake by the plant (supplied_nitrate_nitrogen).
        state.supplied_nitrate_nitrogen += upno3c;
    }
//     Ammonium in the soil is in a dynamic equilibrium between the adsorbed and the soluble
//  fractions. The parameters p1 and p2 are used to compute the dissolved concentration,
//  AmmonNDissolved, of ammoniumnitrogen. bb, cc, ee are intermediate values for computing.
    if (VolNh4NContent[l][k] > 0) {
        double bb = p1 + p2 * soil_cell.water_content - VolNh4NContent[l][k];
        double cc = p2 * soil_cell.water_content * VolNh4NContent[l][k];
        double ee = bb * bb + 4 * cc;
        if (ee < 0)
            ee = 0;
        double AmmonNDissolved; // ammonium N dissolved in soil water, mg cm-3 of soil
        AmmonNDissolved = 0.5 * (sqrt(ee) - bb);
//     Uptake of ammonium N is now computed from AmmonNDissolved , using the Michaelis-Menten
//  method, as for nitrate. upnh4c is added to the total uptake supplied_ammonium_nitrogen.
        if (AmmonNDissolved > 0) {
            double upnh4c; // uptake rate of ammonium, g N per plant per day
            upnh4c = reqnc * AmmonNDissolved / (halfn * soil_cell.water_content + AmmonNDissolved);
            double upmax; // maximum possible uptake rate, mg N per soil cell per day
            upmax = cparupmax * VolNh4NContent[l][k];
            if ((coeff * upnh4c) < upmax)
                VolNh4NContent[l][k] -= coeff * upnh4c;
            else {
                VolNh4NContent[l][k] -= upmax;
                upnh4c = upmax / coeff;
            }
            state.supplied_ammonium_nitrogen += upnh4c;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void WaterFlux(double q1[], double psi1[], double dd[], double qr1[],
               double qs1[], double pp1[], int nn, int iv, int ll, long numiter)
//     This function computes the movement of water in the soil, caused by potential differences
//  between cells in a soil column or in a soil layer. It is called by function
//  CapillaryFlow(). It calls functions WaterBalance(), psiq(), qpsi() and wcond().
//
//     The following arguments are used:
//       q1 = array of volumetric water content, v/v.
//       psi1 = array of matric soil water potential, bars.
//       dd1 = array of widths of soil cells in the direction of flow, cm.
//       qr1 = array of residual volumetric water content.
//       qs1 = array of saturated volumetric water content.
//       pp1 = array of pore space v/v.
//       nn = number of cells in the array.
//       iv = indicator if flow direction, iv = 1 for vertical iv = 0 for horizontal.
//       ll = layer number if the flow is horizontal.
//       numiter = counter for the number of iterations.
//
//     Global variables referenced:
//       alpha, vanGenuchtenBeta, RatioImplicit, SaturatedHydCond, SoilHorizonNum.
//
{
    double delt = 1 / (double) noitr; // the time step of this iteration (fraction of day)
    double cond[40]; // values of hydraulic conductivity
    double kx[40]; // non-dimensional conductivity to the lower layer or to the column on the right
    double ky[40]; // non-dimensional conductivity to the upper layer or to the column on the left
//     Loop over all soil cells. if this is a vertical flow, define the profile index j
//  for each soil cell. compute the hydraulic conductivity of each soil cell, using the
//  function wcond(). Zero the arrays kx and ky.
    int j = SoilHorizonNum[ll]; // for horizontal flow (iv = 0)
    for (int i = 0; i < nn; i++) {
        if (iv == 1)
            j = SoilHorizonNum[i]; // for vertical flow
        cond[i] = wcond(q1[i], qr1[i], qs1[i], vanGenuchtenBeta[j], SaturatedHydCond[j], pp1[i]);
        kx[i] = 0;
        ky[i] = 0;
    }
//     Loop from the second soil cell. compute the array dy (distances
//  between the midpoints of adjacent cells).
//     Compute the average conductivity avcond[i] between soil cells
//  i and (i-1). for low values of conductivity in both cells,(( an
//  arithmetic mean is computed)). for higher values the harmonic mean is
//  used, but if in one of the cells the conductivity is too low (less
//  than a minimum value of condmin ), replace it with condmin.
//
    double dy[40]; // the distance between the midpoint of a layer (or a column) and the midpoint
    // of the layer above it (or the column to the left of it)
    double avcond[40]; // average hydraulic conductivity of two adjacent soil cells
    double condmin = 0.000006;  // minimum value of conductivity, used for computing averages
    for (int i = 1; i < nn; i++) {
        dy[i] = 0.5 * (dd[i - 1] + dd[i]);
        if (cond[i - 1] <= condmin && cond[i] <= condmin)
            avcond[i] = condmin;
        else if (cond[i - 1] <= condmin && cond[i] > condmin)
            avcond[i] = 2 * condmin * cond[i] / (condmin + cond[i]);
        else if (cond[i] <= condmin && cond[i - 1] > condmin)
            avcond[i] = 2 * condmin * cond[i - 1] / (condmin + cond[i - 1]);
        else
            avcond[i] = 2 * cond[i - 1] * cond[i] / (cond[i - 1] + cond[i]);
    }
//     The numerical solution of the flow equation is a combination of the implicit method
//  (weighted by RatioImplicit) and the explicit method (weighted by 1-RatioImplicit).
//     Compute the explicit part of the solution, weighted by (1-RatioImplicit).
//  store water content values, before changing them, in array qx.
    double qx[40]; // previous value of q1.
    double addq[40]; // water added to qx
    double sumaddq = 0; // sum of addq
    for (int i = 0; i < nn; i++)
        qx[i] = q1[i];
//     Loop from the second to the last but one soil cells.
    for (int i = 1; i < nn - 1; i++) {
//     Compute the difference in soil water potential between adjacent cells (deltpsi).
//  This difference is not allowed to be greater than 1000 bars, in order to prevent computational
//  overflow in cells with low water content.
        double deltpsi = psi1[i - 1] - psi1[i]; // difference of soil water potentials (in bars)
        // between adjacent soil soil cells
        if (deltpsi > 1000)
            deltpsi = 1000;
        if (deltpsi < -1000)
            deltpsi = -1000;
//     If this is a vertical flux, add the gravity component of water potential.
        if (iv == 1)
            deltpsi += 0.001 * dy[i];
//     Compute dumm1 (the hydraulic conductivity redimensioned to cm), and check that it will
//  not exceed conmax multiplied by the distance between soil cells, in order to prevent
//  overflow errors.
        double dumm1; // redimensioned hydraulic conductivity components between adjacent cells.
        dumm1 = 1000 * avcond[i] * delt / dy[i];
        if (dumm1 > conmax * dy[i])
            dumm1 = conmax * dy[i];
//     Water entering soil cell i is now computed, weighted by (1 - RatioImplicit).
//  It is not allowed to be greater than 25% of the difference between the cells.
//     Compute water movement from soil cell i-1 to i:
        double dqq1; // water added to cell i from cell (i-1)
        dqq1 = (1 - RatioImplicit) * deltpsi * dumm1;
        double deltq; // difference of soil water content (v/v) between adjacent cells.
        deltq = qx[i - 1] - qx[i];
        if (fabs(dqq1) > fabs(0.25 * deltq)) {
            if (deltq > 0 && dqq1 < 0)
                dqq1 = 0;
            else if (deltq < 0 && dqq1 > 0)
                dqq1 = 0;
            else
                dqq1 = 0.25 * deltq;
        }
//     This is now repeated for water movement from i+1 to i.
        deltpsi = psi1[i + 1] - psi1[i];
        deltq = qx[i + 1] - qx[i];
        if (deltpsi > 1000)
            deltpsi = 1000;
        if (deltpsi < -1000)
            deltpsi = -1000;
        if (iv == 1)
            deltpsi -= 0.001 * dy[i + 1];
        dumm1 = 1000 * avcond[i + 1] * delt / dy[i + 1];
        if (dumm1 > (conmax * dy[i + 1]))
            dumm1 = conmax * dy[i + 1];
        double dqq2 = (1 - RatioImplicit) * deltpsi * dumm1; // water added to cell i from cell (i+1)
        if (fabs(dqq2) > fabs(0.25 * deltq)) {
            if (deltq > 0 && dqq2 < 0)
                dqq2 = 0;
            else if (deltq < 0 && dqq2 > 0)
                dqq2 = 0;
            else
                dqq2 = 0.25 * deltq;
        }
        addq[i] = (dqq1 + dqq2) / dd[i];
        sumaddq += dqq1 + dqq2;
//     Water content of the first and last soil cells is
//  updated to account for flow to or from their adjacent soil cells.
        if (i == 1) {
            addq[0] = -dqq1 / dd[0];
            sumaddq -= dqq1;
        }
        if (i == nn - 2) {
            addq[nn - 1] = -dqq2 / dd[nn - 1];
            sumaddq -= dqq2;
        }
    }
//     Water content q1[i] and soil water potential psi1[i] are updated.
    for (int i = 0; i < nn; i++) {
        q1[i] = qx[i] + addq[i];
        if (iv == 1)
            j = SoilHorizonNum[i];
        psi1[i] = psiq(q1[i], qr1[i], qs1[i], alpha[j], vanGenuchtenBeta[j]);
    }
//     Compute the implicit part of the solution, weighted by RatioImplicit, starting
//  loop from the second cell.
    for (int i = 1; i < nn; i++) {
//     Mean conductivity (avcond) between adjacent cells is made "dimensionless" (ky) by
//  multiplying it by the time step (delt)and dividing it by cell length (dd) and by dy.
//  It is also multiplied by 1000 for converting the potential differences from bars to cm.
        ky[i] = 1000 * avcond[i] * delt / (dy[i] * dd[i]);
//     Very low values of ky are converted to zero, to prevent underflow computer errors, and
//  very high values are converted to maximum limit (conmax), to prevent overflow errors.
        if (ky[i] < 0.0000001)
            ky[i] = 0;
        if (ky[i] > conmax)
            ky[i] = conmax;
    }
//     ky[i] is the conductivity between soil cells i and i-1, whereas kx[i] is between i and i+1.
//  Another loop, until the last but one soil cell, computes kx in a similar manner.
    for (int i = 0; i < nn - 1; i++) {
        kx[i] = 1000 * avcond[i + 1] * delt / (dy[i + 1] * dd[i]);
        if (kx[i] < 0.0000001) kx[i] = 0;
        if (kx[i] > conmax) kx[i] = conmax;
    }
//     Arrays used for the implicit numeric solution:
    double a1[40], b1[40], cau[40], cc1[40], d1[40], dau[40];
    for (int i = 0; i < nn; i++) {
//     Arrays a1, b1, and cc1 are computed for the implicit part of
//  the solution, weighted by RatioImplicit.
        a1[i] = -kx[i] * RatioImplicit;
        b1[i] = 1 + RatioImplicit * (kx[i] + ky[i]);
        cc1[i] = -ky[i] * RatioImplicit;
        if (iv == 1) {
            j = SoilHorizonNum[i];
            a1[i] = a1[i] - 0.001 * kx[i] * RatioImplicit;
            cc1[i] = cc1[i] + 0.001 * ky[i] * RatioImplicit;
        }
//     The water content of each soil cell is converted to water
//  potential by function psiq and stored in array d1 (in bar units).
        d1[i] = psiq(q1[i], qr1[i], qs1[i], alpha[j], vanGenuchtenBeta[j]);
    }
//     The solution of the simultaneous equations in the implicit method alternates between
//  the two directions along the arrays. The reason for this is because the direction of the
//  solution may cause some cumulative bias. The counter numiter determines the direction
//  of the solution.
//     The solution in this section starts from the last soil cell (nn).
    if ((numiter % 2) == 0) {
//     Intermediate arrays dau and cau are computed.
        cau[nn - 1] = psi1[nn - 1];
        dau[nn - 1] = 0;
        for (int i = nn - 2; i > 0; i--) {
            double p = a1[i] * dau[i + 1] + b1[i]; // temporary
            dau[i] = -cc1[i] / p;
            cau[i] = (d1[i] - a1[i] * cau[i + 1]) / p;
        }
        if (iv == 1)
            j = SoilHorizonNum[0];
        psi1[0] = psiq(q1[0], qr1[0], qs1[0], alpha[j], vanGenuchtenBeta[j]);
//     psi1 is now computed for soil cells 1 to nn-2. q1 is
//  computed from psi1 by function qpsi.
        for (int i = 1; i < nn - 1; i++) {
            if (iv == 1)
                j = SoilHorizonNum[i];
            psi1[i] = dau[i] * psi1[i - 1] + cau[i];
            q1[i] = qpsi(psi1[i], qr1[i], qs1[i], alpha[j], vanGenuchtenBeta[j]);
        }
    }
//     The alternative direction of solution is executed here. the
//  solution in this section starts from the first soil cell.
    else {
//     Intermediate arrays dau and cau are computed, and the computations
//  described previously are repeated in the opposite direction.
        cau[0] = psi1[0];
        dau[0] = 0;
        for (int i = 1; i < nn - 1; i++) {
            double p = a1[i] * dau[i - 1] + b1[i]; // temporary
            dau[i] = -cc1[i] / p;
            cau[i] = (d1[i] - a1[i] * cau[i - 1]) / p;
        }
        if (iv == 1)
            j = SoilHorizonNum[nn - 1];
        psi1[nn - 1] = psiq(q1[nn - 1], qr1[nn - 1], qs1[nn - 1], alpha[j], vanGenuchtenBeta[j]);
        for (int i = nn - 2; i > 0; i--) {
            if (iv == 1)
                j = SoilHorizonNum[i];
            psi1[i] = dau[i] * psi1[i + 1] + cau[i];
            q1[i] = qpsi(psi1[i], qr1[i], qs1[i], alpha[j], vanGenuchtenBeta[j]);
        }
    }
//     The limits of water content are now checked and corrected, and
//  function WaterBalance() is called to correct water amounts.
    for (int i = 0; i < nn; i++) {
        if (q1[i] < qr1[i])
            q1[i] = qr1[i];
        if (q1[i] > qs1[i])
            q1[i] = qs1[i];
        if (q1[i] > pp1[i])
            q1[i] = pp1[i];
    }
    WaterBalance(q1, qx, dd, nn);
}

////////////////////////
void WaterBalance(double q1[], double qx[], double dd[], int nn)
//     This function checks and corrects the water balance in the soil cells
//  within a soil column or a soil layer. It is called by WaterFlux().
//     The implicit part of the solution may cause some deviation in the total amount of water
//  to occur. This module corrects the water balance if the sum of deviations is not zero, so
//  that the total amount of water in the array will not change. The correction is proportional
//  to the difference between the previous and present water amounts in each soil cell.
//
//     The following arguments are used here:
//        dd[] - one dimensional array of layer or column widths.
//        nn - the number of cells in this layer or column.
//        qx[] - one dimensional array of a layer or a column of the
//               previous values of cell.water_content.
//        q1[] - one dimensional array of a layer or a column of cell.water_content.
//
{
    double dev = 0;  // Sum of differences of water amount in soil
    double dabs = 0; // Sum of absolute value of differences in water content in
    // the array between beginning and end of this time step.
    for (int i = 0; i < nn; i++) {
        dev += dd[i] * (q1[i] - qx[i]);
        dabs += fabs(q1[i] - qx[i]);
    }
    if (dabs > 0)
        for (int i = 0; i < nn; i++)
            q1[i] = q1[i] - fabs(q1[i] - qx[i]) * dev / (dabs * dd[i]);

}

////////////////////////////////////////////////////////////////////////////////
void NitrogenFlow(int nn, double q01[], double q1[], double dd[], double nit[], double nur[])
//     This function computes the movement of nitrate and urea between the soil cells,
//  within a soil column or within a soil layer, as a result of water flux.
//     It is called by function CapillaryFlow().
//     It is assumed that there is only a passive movement of nitrate and urea
//  (i.e., with the movement of water).
//     The following arguments are used here:
//       dd[] - one dimensional array of layer or column widths.
//       nit[] - one dimensional array of a layer or a column of VolNo3NContent.
//       nn - the number of cells in this layer or column.
//       nur[] - one dimensional array of a layer or a column of VolUreaNContent.
//       q01[] - one dimensional array of a layer or a column of the previous values of cell.water_content.
//       q1[] - one dimensional array of a layer or a column of cell.water_content.
//
{
//     Zeroise very small values to prevent underflow.
    for (int i = 0; i < nn; i++) {
        if (nur[i] < 1e-20)
            nur[i] = 0;
        if (nit[i] < 1e-20)
            nit[i] = 0;
    }
//     Declare and zeroise arrays.
    double qdn[40] = {40 * 0}; // amount of nitrate N moving to the previous cell.
    double qup[40] = {40 * 0}; // amount of nitrate N moving to the following cell.
    double udn[40] = {40 * 0}; // amount of urea N moving to the previous cell.
    double uup[40] = {40 * 0}; // amount of urea N moving to the following cell.
    for (int i = 0; i < nn; i++) {
//     The amout of water in each soil cell before (aq0) and after (aq1) water movement is
//  computed from the previous values of water content (q01), the present values (q1),
//  and layer thickness. The associated transfer of soluble nitrate N (qup and qdn) and urea N
//  (uup and udn) is now computed. qup and uup are upward movement (from cell i+1 to i),
//  qdn and udn are downward movement (from cell i-1 to i).
        double aq0 = q01[i] * dd[i]; // previous amount of water in cell i
        double aq1 = q1[i] * dd[i];  // amount of water in cell i now
        if (i == 0) {
            qup[i] = 0;
            uup[i] = 0;
        } else {
            qup[i] = -qdn[i - 1];
            uup[i] = -udn[i - 1];
        }
//
        if (i == nn - 1) {
            qdn[i] = 0;
            udn[i] = 0;
        } else {
            qdn[i] = (aq1 - aq0) * nit[i + 1] / q01[i + 1];
            if (qdn[i] < (-0.2 * nit[i] * dd[i]))
                qdn[i] = -0.2 * nit[i] * dd[i];
            if (qdn[i] > (0.2 * nit[i + 1] * dd[i + 1]))
                qdn[i] = 0.2 * nit[i + 1] * dd[i + 1];
            udn[i] = (aq1 - aq0) * nur[i + 1] / q01[i + 1];
            if (udn[i] < (-0.2 * nur[i] * dd[i]))
                udn[i] = -0.2 * nur[i] * dd[i];
            if (udn[i] > (0.2 * nur[i + 1] * dd[i + 1]))
                udn[i] = 0.2 * nur[i + 1] * dd[i + 1];
        }
    }
//     Loop over all cells to update nit and nur arrays.
    for (int i = 0; i < nn; i++) {
        nit[i] += (qdn[i] + qup[i]) / dd[i];
        nur[i] += (udn[i] + uup[i]) / dd[i];
    }
}

///////////////////////////////////////////////////////////////////////////
void SoilSum(State &state, double row_space)
//     This function computes sums of some soil variables. It is called by SimulateThisDay()
//     It computes the total amounts per slab of water in mm (TotalSoilWater). The sums of
//  nitrogen are in mg N per slab: nitrate (TotalSoilNo3N), ammonium (TotalSoilNh4N),
//  urea (TotalSoilUreaN) and the sum of all mineral nitrogen forms (TotalSoilNitrogen).
//
//     The following global variables are referenced here:
//       dl, nk, nl, RowSpace, VolNh4NContent, VolNo3NContent, VolUreaNContent, wk.
//     The following global variables are set here:
//       TotalSoilNitrogen, TotalSoilWater, TotalSoilNh4N, TotalSoilNo3N, TotalSoilUreaN.
{
    TotalSoilWater = 0;
    TotalSoilNo3N = 0;
    TotalSoilNh4N = 0;
    TotalSoilUreaN = 0;
//
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++) {
            TotalSoilNo3N += state.soil.cells[l][k].nitrate_nitrogen_content * dl(l) * wk(k, row_space);
            TotalSoilNh4N += VolNh4NContent[l][k] * dl(l) * wk(k, row_space);
            TotalSoilUreaN += VolUreaNContent[l][k] * dl(l) * wk(k, row_space);
            TotalSoilWater += state.soil.cells[l][k].water_content * dl(l) * wk(k, row_space);
        }
//
    TotalSoilNitrogen = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
    TotalSoilWater = TotalSoilWater * 10 / row_space;
}
