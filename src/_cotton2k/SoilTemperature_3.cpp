//  File SoilTemperature_3.cpp
//
//   List of functions in this file:
//       CanopyBalance()
//       SoilHeatFlux()
//       ThermalCondSoil()
//       HeatBalance()
//
#include <cmath>
#include "Simulation.hpp"
#include "global.h"
#include "exceptions.h"
#include "GeneralFunctions.h"
#include "Input.h"

using namespace std;

double ThermalCondSoil(double, double, int);

void HeatBalance(int nn);

//  arrays with file scope:
double dz[maxl];   // equal to the dl array in a columnn, or wk in a row.
double ts1[maxl];  // array of soil temperatures.
double ts0[maxl];  // array of previous soil temperatures.
double hcap[maxl]; // heat capacity of soil layer (cal cm-3 oC-1).
/////////////////////////////////////////
void CanopyBalance(int ihr, int k, double etp1, double rlzero, double rsv,
                   double c2, double sf, double so, double thet, double &tv, const int &Daynum)
//     This function solves the energy balance equations at the foliage / air interface, and
//  computes the resulting temperature of the foliage. It is Called from EnergyBalance().
//     Units for all energy fluxes are: cal cm-2 sec-1.
//
//     The following argument is set in this function:
//       tv - temperature of plant canopy.
//     The following arguments are referenced in this function:
//       c2 - multiplier for sensible heat transfer at plant surface.
//       etp1 - transpiration (mm / sec).
//       ihr - the time in hours.
//       k - soil column number.
//       rlzero - incoming long wave radiation (ly / sec).
//       rsv - global radiation absorbed by the vegetation
//       sf - fraction of shaded soil area
//       so - temperature of soil surface (k).
//       thet - air temperature (k).
{
    //     Constants:
    const double ef = 0.95;                      // emissivity of the foliage surface
    const double eg = 0.95;                      // emissivity of the soil surface
    const double stefa1 = 1.38e-12;              // stefan-boltsman constant.
                                                 //
    double rlv1;                                 // long wave radiation reaching the canopy
    rlv1 = sf * ef * rlzero                      // from sky
           + sf * ef * eg * stefa1 * pow(so, 4); // from soil
                                                 //     rlv4 is the multiplier of tv**4 for emitted long wave radiation from vegetation, corrected
                                                 //  for the amount reflected back from soil surface and absorbed by foliage. This is two-sided
                                                 //  ( note that when eg = ef = 1, the coefficient corr will be 2)
    double corr = 1 + eg / (ef + eg - ef * eg);  // coefficient
    double rlv4 = stefa1 * sf * ef * corr;
    double dsenfheat = c2 * (1 - 0.6 * sf); // derivative of senfheat
                                            //     Start iterations for tv:
    int mot = 0;                            // count of iterations for plant surface
    while (mot < 50)
    {
        //     Latent heat flux (hvlat) is computed from the transpiration rate.
        double hvlat = (75.5255 - 0.05752 * tv) * etp1;
        double dhvlat = -0.05752 * etp1; // derivative of hvlat
                                         //     Emitted long wave radiation from vegetation (cclwe)
        double cclwe = rlv4 * pow(tv, 4);
        //     Sensible heat transfer from vegetation
        double tvex = tv; // previous value of tv
        double tafk;      // average air temperature above soil surface (k) in canopy
        tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv);
        double senfheat; // sensible heat transfer from foliage
        senfheat = c2 * (tv - tafk);
        //     Compute the energy balance at the plant surface (cc), and
        //  if it is small enough end the computation.
        double cc = cclwe       // (a) long wave emission from vegetation
                    - rlv1      // long wave radiation reaching the vegetation
                    + hvlat     // (b) latent heat transfer
                    - rsv       // global radiation on vegetation
                    + senfheat; // (c) sensible heat transfer from vegetation to air
        if (fabs(cc) < 10e-6)
        {
            tv = 0.5 * (tvex + tv);
            return; // end iterations for tv
        }
        //     If cc is not small enough, compute its derivative by tv (ccp).
        double dcclwe = 4 * rlv4 * pow(tv, 3); // derivative of cclwe
                                               //     ccp is the derivative of energy balance at the plant surface (by tv)
        double ccp = dcclwe                    // (a)
                     + dhvlat                  // (b)
                     + dsenfheat;              // (c)
                                               //     Correct the canopy temperature by  the ratio of cc to ccp.
        double ccadjust = cc / ccp;            // adjustment of tv before next iteration
                                               //     If adjustment is small enough, no more iterations are needed.
        if (fabs(ccadjust) < 0.002)
            return;
        //     If ccadjust is not the same sign as ccadx, reduce fluctuations
        double ccadx; // previous value of ccadjust
        if (mot >= 2)
        {
            if (fabs(ccadjust - ccadx) > fabs(ccadjust + ccadx))
            {
                ccadjust = (ccadjust + ccadx) * 0.5;
                tv = (tv + tvex) * 0.5;
            }
        }
        if (ccadjust > 10)
            ccadjust = 10;
        if (ccadjust < -10)
            ccadjust = -10;
        tv = tv - ccadjust;
        tvex = tv;
        ccadx = ccadjust;
        mot++;
    }
    //     If reached 50 iterations there must be an error somewhere!
    string msg = " Infinite loop in CanopyBalance(). Abnormal stop!! \n";
    char C1[12];
    sprintf(C1, "%3d %3d %3d", Daynum, ihr, k);
    msg += " Daynum, ihr, k = " + (string)C1 + "\n";
    sprintf(C1, "%10.3g", so);
    msg += " so      = " + (string)C1 + "\n";
    sprintf(C1, "%10.3g", tv);
    msg += " tv = " + (string)C1 + "\n";
    throw SimulationEnd(msg);
}

///////////////////////////////////////////////////////////////////
void SoilHeatFlux(State &state, double dlt, int iv, int nn, int layer, int n0, double row_space)
//     This function computes heat flux in one direction between soil cells.
//  It is called from SoilTemperature(), and calls ThermalCondSoil() and HeatBalance().
//     Note: the units are: thermal conductivity = cal cm-1 s-1 oC-1;
//                          heat capacity = cal cm-3 oC-1;
//                      	thermal diffusivity = cm2 s-1;
//                          ckx and cky are dimensionless;
//
//     The following global variables are referenced here:
//       dl, HeatCapacitySoilSolid, maxl, PoreSpace.
//     The following global variable is set here:    SoilTemp.
//     The following input arguments are used in this function:
//       dlt - time (seconds) of one iteration.
//       iv -  = 1 for vertical flux, = 0 for horizontal flux.
//       layer - soil layer number.
//       n0 - number of layer or column of this array
//       nn - number of soil cells in the array.
{
    //     Constant parameters:
    const double beta1 = 0.90; // weighting factor for the implicit method of computation.
    const double ca = 0.0003;  // heat capacity of air (cal cm-3 oC-1).
                               //     Set soil layer number l (needed to define HeatCapacitySoilSolid, PoreSpace, ThermalCondSoil).
                               //  Compute for each soil cell the heat capacity and heat diffusivity.
    static long numiter = 0;   // number of this iteration.
    int l = layer;             // soil layer number.
    double q1[maxl];           // array of water content.
    double asoi[maxl];         // array of thermal diffusivity of soil cells (cm2 s-1).
    for (int i = 0; i < nn; i++)
    {
        if (iv == 1)
        {
            l = i;
            q1[i] = state.soil.cells[i][n0].water_content;
            ts1[i] = SoilTemp[i][n0];
            dz[i] = dl(i);
        }
        else
        {
            q1[i] = state.soil.cells[n0][i].water_content;
            ts1[i] = SoilTemp[n0][i];
            dz[i] = wk(i, row_space);
        }
        hcap[i] = HeatCapacitySoilSolid[l] + q1[i] + (PoreSpace[l] - q1[i]) * ca;
        asoi[i] = ThermalCondSoil(q1[i], ts1[i], l) / hcap[i];
    }
    //     The numerical solution of the flow equation is a combination of the implicit method
    //  (weighted by beta1) and the explicit method (weighted by 1-beta1).
    double dltt;         // computed time step required.
    double avdif[maxl];  // average thermal diffusivity between adjacent cells.
    double dy[maxl];     // array of distances between centers of adjacent cells (cm).
    double dltmin = dlt; // minimum time step for the explicit solution.
    avdif[0] = 0;
    dy[0] = 0;
    for (int i = 1; i < nn; i++)
    {
        //     Compute average diffusivities avdif between layer i and the previous (i-1),
        //  and dy(i), distance (cm) between centers of layer i and the previous (i-1)
        avdif[i] = (asoi[i] + asoi[i - 1]) / 2;
        dy[i] = (dz[i - 1] + dz[i]) / 2;
        //     Determine the minimum time step required for the explicit solution.
        dltt = 0.2 * dy[i] * dz[i] / avdif[i] / (1 - beta1);
        if (dltt < dltmin)
            dltmin = dltt;
    }
    //     Use time step of dlt1 seconds, for iterx iterations
    int iterx = (int)(dlt / dltmin); // computed number of iterations.
    if (dltmin < dlt)
        iterx++;
    double dlt1 = dlt / iterx; // computed time (seconds) of an iteration.
                               // start iterations. Store temperature data in array ts0. count iterations.
    for (int ii = 0; ii < iterx; ii++)
    {
        for (int i = 0; i < nn; i++)
        {
            ts0[i] = ts1[i];
            if (iv == 1)
                l = i;
            asoi[i] = ThermalCondSoil(q1[i], ts1[i], l) / hcap[i];
            if (i > 0)
                avdif[i] = (asoi[i] + asoi[i - 1]) / 2;
        }
        numiter++;
        //     The solution of the simultaneous equations in the implicit method alternates between
        //  the two directions along the arrays. The reason for this is because the direction of the
        //  solution may cause some cumulative bias. The counter numiter determines the direction
        //  of the solution.
        double cau[maxl], dau[maxl]; // arrays used for the implicit numerical solution.
        double ckx, cky;             // nondimensional diffusivities to next and previous layers.
        double vara, varb;           // used for computing the implicit solution.
        if ((numiter % 2) == 0)
        {
            //     1st direction of computation, for an even iteration number:
            dau[0] = 0;
            cau[0] = ts1[0];
            //     Loop from the second to the last but one soil cells. Compute
            //  nondimensional diffusivities to next and previous layers.
            for (int i = 1; i < nn - 1; i++)
            {
                ckx = avdif[i + 1] * dlt1 / (dz[i] * dy[i + 1]);
                cky = avdif[i] * dlt1 / (dz[i] * dy[i]);
                //     Correct value of layer 1 for explicit heat movement to/from layer 2
                if (i == 1)
                    cau[0] = ts1[0] - (1 - beta1) * (ts1[0] - ts1[1]) * cky * dz[1] / dz[0];
                vara = 1 + beta1 * (ckx + cky) - beta1 * ckx * dau[i - 1];
                dau[i] = beta1 * cky / vara;
                varb = ts1[i] + (1 - beta1) * (cky * ts1[i - 1] + ckx * ts1[i + 1] - (cky + ckx) * ts1[i]);
                cau[i] = (varb + beta1 * ckx * cau[i - 1]) / vara;
            }
            //     Correct value of last layer (nn-1) for explicit heat movement to/from layer nn-2
            ts1[nn - 1] = ts1[nn - 1] - (1 - beta1) * (ts1[nn - 1] - ts1[nn - 2]) * ckx * dz[nn - 2] / dz[nn - 1];
            //     Continue with the implicit solution
            for (int i = nn - 2; i >= 0; i--)
                ts1[i] = dau[i] * ts1[i + 1] + cau[i];
        }
        else
        {
            //     Alternate direction of computation for odd iteration number
            dau[nn - 1] = 0;
            cau[nn - 1] = ts1[nn - 1];
            for (int i = nn - 2; i > 0; i--)
            {
                ckx = avdif[i + 1] * dlt1 / (dz[i] * dy[i + 1]);
                cky = avdif[i] * dlt1 / (dz[i] * dy[i]);
                if (i == nn - 2)
                    cau[nn - 1] =
                        ts1[nn - 1] - (1 - beta1) * (ts1[nn - 1] - ts1[nn - 2]) * ckx * dz[nn - 2] / dz[nn - 1];
                vara = 1 + beta1 * (ckx + cky) - beta1 * cky * dau[i + 1];
                dau[i] = beta1 * ckx / vara;
                varb = ts1[i] + (1 - beta1) * (ckx * ts1[i + 1] + cky * ts1[i - 1] - (cky + ckx) * ts1[i]);
                cau[i] = (varb + beta1 * cky * cau[i + 1]) / vara;
            }
            ts1[0] = ts1[0] - (1 - beta1) * (ts1[0] - ts1[1]) * cky * dz[1] / dz[0];
            for (int i = 1; i < nn; i++)
                ts1[i] = dau[i] * ts1[i - 1] + cau[i];
        } // end if numiter
          //     Call HeatBalance to correct quantitative deviations caused by the
          //  imlicit part of the solution.
        HeatBalance(nn);
    } // end of iterx loop
      //     Set values of SoiTemp
    for (int i = 0; i < nn; i++)
    {
        if (iv == 1)
            SoilTemp[i][n0] = ts1[i];
        else
            SoilTemp[n0][i] = ts1[i];
    }
}

/////////////////////////////////////////
double ThermalCondSoil(double q0, double t0, int l0)
//     This function computes and returns the thermal conductivity of the soil
//  (cal cm-1 s-1 oC-1). It is based on the work of De Vries(1963).
//
//     The following global variables are referenced here:
//       ClayVolumeFraction, FieldCapacity, HeatCondDrySoil, MarginalWaterContent, PoreSpace,
//       SandVolumeFraction.
//     The following arguments are used in this function:
//       l0 - soil layer.
//       q0 - volumetric soil moisture content.
//       t0 - soil temperature (K).
{
    //     Constant parameters:
    const double bclay = 7.0;   // heat conductivity of clay (= 7 mcal cm-1 s-1 oc-1).
    const double bsand = 20.0;  // heat conductivity of sand (= 20 mcal cm-1 s-1 oc-1).
    const double cka = 0.0615;  // heat conductivity of air (= 0.0615 mcal cm-1 s-1 oc-1).
    const double ckw = 1.45;    // heat conductivity of water (= 1.45 mcal cm-1 s-1 oc-1).
                                //     Convert soil temperature to degrees C.
    double tcel = t0 - 273.161; // soil temperature, in C.
                                //     Compute cpn, the apparent heat conductivity of air in soil pore spaces, when saturated with
                                //  water vapor, using a function of soil temperature, which changes linearly between 36 and 40 C.
    double bb;                  // effect of temperature on heat conductivity of air saturated with water vapor.
    if (tcel <= 36)
        bb = 0.06188;
    else if (tcel > 36 && tcel <= 40)
        bb = 0.06188 + (tcel - 36) * (0.05790 - 0.06188) / (40 - 36);
    else
        bb = 0.05790;
    double cpn; // apparent heat conductivity of air in soil pore spaces, when it is saturated with water vapor.
    cpn = cka + 0.05 * exp(bb * tcel);
    //     Compute xair, air content of soil per volume, from soil porosity and moisture content.
    //     Compute thermal conductivity (a) for wet soil (soil moisture higher than field capacity),
    //                                  (b) for less wet soil.
    //     In each case compute first ga, and then dair.
    double xair; // air content of soil, per volume.
    xair = PoreSpace[l0] - q0;
    if (xair < 0)
        xair = 0;
    double dair;  // aggregation factor for air in soil pore spaces.
    double ga;    // shape factor for air in pore spaces.
    double hcond; // computed heat conductivity of soil, mcal cm-1 s-1 oc-1.
    if (q0 >= FieldCapacity[l0])
    {
        //     (a) Heat conductivity of soil wetter than field capacity.
        ga = 0.333 - 0.061 * xair / PoreSpace[l0];
        dair = form(cpn, ckw, ga);
        hcond = (q0 * ckw + dsand * bsand * SandVolumeFraction[l0] + dclay * bclay * ClayVolumeFraction[l0] + dair * cpn * xair) / (q0 + dsand * SandVolumeFraction[l0] + dclay * ClayVolumeFraction[l0] + dair * xair);
    }
    else
    {
        //     (b) For soil less wet than field capacity, compute also ckn (heat conductivity
        // of air in the soil pores).
        double qq;  // soil water content for computing ckn and ga.
        double ckn; // heat conductivity of air in pores in soil.
        qq = q0;
        if (qq < MarginalWaterContent[l0])
            qq = MarginalWaterContent[l0];
        ckn = cka + (cpn - cka) * qq / FieldCapacity[l0];
        ga = 0.041 + 0.244 * (qq - MarginalWaterContent[l0]) / (FieldCapacity[l0] - MarginalWaterContent[l0]);
        dair = form(ckn, ckw, ga);
        hcond = (qq * ckw + dsand * bsand * SandVolumeFraction[l0] + dclay * bclay * ClayVolumeFraction[l0] + dair * ckn * xair) / (qq + dsand * SandVolumeFraction[l0] + dclay * ClayVolumeFraction[l0] + dair * xair);
        //     When soil moisture content is less than the limiting value MarginalWaterContent,
        //  modify the value of hcond.
        if (qq <= MarginalWaterContent[l0])
            hcond = (hcond - HeatCondDrySoil[l0]) * q0 / MarginalWaterContent[l0] + HeatCondDrySoil[l0];
    } // q0
      //     The result is hcond converted from mcal to cal.
    double result = hcond / 1000;
    return result;
}

////////////////////////////////////////////////////////////////////////////////////
void HeatBalance(int nn)
//     This function checks and corrects the heat balance in the soil soil cells, within
//  a soil layer. It is called by function SoilHeatFlux() only for horizontal flux.
//     The implicit part of the solution may cause some deviation in the total heat sum
//  to occur. This module corrects the heat balance if the sum of absolute deviations is
//  not zero, so that the total amount of heat in the array does not change. The correction
//  is proportional to the difference between the previous and present heat amounts.
//
//     The following arguments are referenced here:
//       nn - the number of soil cells in this layer or column.
//     The following global or file scope variables are referenced:
//       dz, hcap, ts0
//     The following file scope variable is set:    ts1
{
    double dabs = 0; //Sum of absolute value of differences in heat content in
    // the array between beginning and end of this time step.
    double dev = 0; // Sum of differences of heat amount in soil.
    for (int i = 0; i < nn; i++)
    {
        dev += dz[i] * hcap[i] * (ts1[i] - ts0[i]);
        dabs += fabs(ts1[i] - ts0[i]);
    }
    if (dabs > 0)
    {
        for (int i = 0; i < nn; i++)
            ts1[i] = ts1[i] - fabs(ts1[i] - ts0[i]) * dev / (dabs * dz[i] * hcap[i]);
    }
}
