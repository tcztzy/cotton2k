//  File SoilTemperature_2.cpp
//
//   List of functions in this file:
//       SensibleHeatTransfer()
//       SoilSurfaceBalance()
//
#include <math.h>
#include <stdio.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
/////////////////////////////////////////////////////////////////////////////
double SensibleHeatTransfer(double tsf, double tenviron, double height,
                            double wndcanp)
//     This function computes the sensible heat transfer coefficient, using the
//     friction
//  potential (shear) temperature (thstar), and the surface friction (shear)
//  velocity (ustar) at the atmospheric boundary. It is called three times from
//  EnergyBalance(): for canopy, soil surface, or mulch boundaries with their
//  environment.
//
//     The following arguments are referenced in this function:
//       tenviron - temperature (K) of the environment - air at 200 cm height
//                  for columns with no canopy, or tafk when canopy is present .
//       tsf - surface temperature, K, of soil, mulch or canopy.
//       wndcanp - wind speed in the canopy (if present), cm s-1.
//       height - canopy height, cm, or zero for soil surface and mulch
//       boundaries.
//     The return value computed in this function:
//       sensibleHeatTransfer - raw sensible heat transfer coefficient
//     Global variable set:   bEnd
{
    //	    Constant values used:
    const double grav = 980;     // acceleration due to gravity (980 cm sec-2).
    const double s40 = 0.13;     // calibration constant.
    const double s42 = 0.63;     // calibration constant.
    const double stmin = 5;      // minimal value of ustar.
    const double vonkar = 0.40;  // Von-Karman constant (0.40).
    const double zalit1 = 0.0962;  // parameter .....
    //     Wind velocity not allowed to be less than 100 cm s-1.
    double u = wndcanp;  // wind speed at 200 cm height, cm / s.
    if (u < 100) u = 100;
    //     Assign initial values to z0 and gtop, and set dt.
    double z0 = s40 * height;  // surface roughness parameter, cm.
    if (z0 < 1) z0 = 1;
    double gtop = log((200 - s42 * height) /
                      z0);  // logarithm of ratio of height of measurement to
                            // surface roughness parameter.
    double dt = tsf - tenviron;  // temperature difference.
    //     Set approximate initial values for ustar and thstar (to reduce
    //     iterations).
    double thstar;  // friction potential (shear) temperature.
    double ustar;   // Surface friction (shear) velocity (cm sec-1).
    if (dt >= 0) {
        ustar = 1.873 + 0.570172 * dt + .07438568 * u;
        thstar = -0.05573 * dt + 1.987 / u - 6.657 * dt / u;
    } else {
        ustar = -4.4017 + 1.067 * dt + 0.25957 * u - 0.001683 * dt * u;
        if (ustar < 5) ustar = 5;
        thstar = -0.0096 - 0.1149 * dt + 0.0000377 * u + 0.0002367 * dt * u;
        if (thstar < 0.03) thstar = 0.03;
    }
    double tbot1 = tsf;  // surface temperature corrected for friction (shear)
                         // potential temperature.
    long mtest = 0;      // count of iterations.
    double g1;           // temporary derived variable
    double ug1chk;       // previous value of ug1.
    double ug1;          // ratio of ustar to g1.
    double ug1res;       // residual value of ug1.
                         //     Start iterations.
    do {
        double thekz;  // previous value of thstar.
        double uchek;  // previous value of ustar.
        ug1chk = 0;    // previous value of UG1.
        if (mtest > 0) {
            //     Assign values to tbot1, uchek, thekz, and ug1chk.
            tbot1 =
                tsf + zalit1 * thstar * pow((ustar * z0 / 15), 0.45) / vonkar;
            uchek = ustar;
            thekz = thstar;
            if (g1 != 0) ug1chk = ustar / g1;
        }
        //     Compute air temperature at 1 cm,, and compute zl and lstar.
        double zl;  // nondimensional height.
        if (fabs(thstar) < 1e-30)
            zl = 0;
        else {
            double thetmn = (tenviron + tbot1) *
                            0.5;  // mean temperature (K) of air and surface.
            double lstar = (thetmn * ustar * ustar) / (vonkar * grav * thstar);
            zl = (200 - s42 * height) / lstar;
            //     Added to decrease fluctuations in zl:
            if (zl < -5) zl = -5;
            if (zl > 0.5) zl = 0.5;
        }
        //     Compute g1u, and g2.
        double g1u, g2;  // temporary derived variables.
        if (zl > 0) {
            g1u = -4.7 * zl;
            g2 = -6.35135 * zl;
            if (g2 < -1) g2 = -1;
        } else {
            double tmp1 = pow((1 - 15 * zl), 0.25);  // intermediate variable.
            g1u = 2 * log((1 + tmp1) / 2) + log((1 + tmp1 * tmp1) / 2) -
                  2 * atan(tmp1 + 1.5708);
            g2 = 2 * log((1 + sqrt(1 - 9 * zl)) / 2);
        }
        if (g2 > gtop) g2 = gtop;
        //     Compute ustar and check for minimum value.
        ustar = vonkar * u / (gtop - g1u);
        if (ustar < stmin) ustar = stmin;
        //     Compute g1 and thstar.
        g1 = 0.74 * (gtop - g2) + zalit1 * pow((ustar * z0 / 0.15), 0.45);
        thstar = -dt * vonkar / g1;
        //     If more than 30 iterations, reduce fluctuations
        if (mtest > 30) {
            thstar = (thstar + thekz) * 0.5;
            ustar = (ustar + uchek) * 0.5;
        }
        mtest++;
        //     Stop simulation if no convergence after 100 iterations.
        if (mtest > 100) {
            fprintf(
                stderr,
                " Infinite loop in SensibleHeatTransfer(). Abnormal stop!! \n");
            fprintf(stderr, " tenviron = %10.3g\n", tenviron);
            fprintf(stderr, " tsf      = %10.3g\n", tsf);
            fprintf(stderr, " PlantHeight = %10.3g\n", height);
            fprintf(stderr, " u = %10.3g\n", u);
            bEnd = true;
            return 0;
        }
        //     Compute ug1 and  ug1res to check convergence
        ug1 = ustar / g1;
        if (fabs(ug1chk) <= 1.e-30)
            ug1res = fabs(ug1);
        else
            ug1res = fabs((ug1chk - ug1) / ug1chk);
        //     If ug1 did not converge, go to next iteration.
    } while (fabs(ug1 - ug1chk) > 0.05 && ug1res > 0.01);
    //     No more iterations. Compute sensibleHeatTransfer.
    double sensibleHeatTransfer = ustar * vonkar / g1;
    return sensibleHeatTransfer;
}
/////////////////////////////////////////
void SoilSurfaceBalance(int ihr, int k, double ess, double rlzero, double rss,
                        double sf, double hsg, double &so, double &so2,
                        double &so3, double thet, double tm, double tv)
//     This function is called from EnergyBalance(). It calls function
//     ThermalCondSoil(). It solves the energy balance equations at the soil
//     surface, and
//  computes the resulting temperature of the soil surface.
//     Units for all energy fluxes are: cal cm-2 sec-1.
//
//     The following global variables are referenced here:
//       Daynum, dl, MulchTranLW, VolWaterContent.
//     The following global variable is set here:      bEnd.
//     The following arguments are set in this function:
//       so - temperature of soil surface.
//       so2 - temperature of soil 2nd layer
//       so3 - temperature of soil 3rd layer
//     The following arguments are referenced in this function:
//       ess - evaporation from soil surface (mm / sec).
//       hsg - multiplier for computing sensible heat transfer from soil to air
//             or from soil to mulch if this is a mulched column.
//       ihr - the time in hours.
//       k - soil column number.
//       rlzero - incoming long wave radiation (ly / sec).
//       rss - global radiation absorbed by soil surface
//       sf - fraction of shaded soil area
//       thet - air temperature (K).
//       tm - temperature of plastic mulch (K). When tm = 0 there is no plastic
//       mulch. tv - temperature of plant canopy (K).
//
{
    // Constants:
    const double ef = 0.95;          // emissivity of the foliage surface
    const double eg = 0.95;          // emissivity of the soil surface
    const double stefa1 = 1.38e-12;  // Stefan-Boltsman constant.
    //     Long wave radiation reaching the soil:
    double rls1;                       // long wave energy reaching soil surface
    if (sf >= 0.05)                    // shaded column
        rls1 = (1 - sf) * eg * rlzero  // from sky on unshaded soil surface
               + sf * eg * ef * stefa1 *
                     pow(tv, 4);  // from foliage on shaded soil surface
    else
        rls1 = eg * rlzero;  // from sky in unshaded column
    if (tm > 0)  // modify by mulch transmissivity and add lw emitted from mulch
        rls1 =
            rls1 * MulchTranLW + eg * (1 - MulchTranLW) * stefa1 * pow(tm, 4);
    //     rls4 is the multiplier of so**4 for emitted long wave radiation from
    //     soil
    double rls4 = eg * stefa1;
    double bbex;       // previous value of bbadjust.
    int mon = 0;       // count of iterations for soil surface energy balance
    double soex = so;  //  previous value of so.
                       //     Start itrations for soil surface enegy balance.
    while (mon < 50) {
        //     Compute latent heat flux from soil evaporation: convert from mm
        //     sec-1 to
        //  cal cm-2 sec-1. Compute derivative of hlat
        //     hlat is the energy used for evaporation from soil surface (latent
        //     heat)
        double hlat = (75.5255 - 0.05752 * so) * ess;
        double dhlat = -0.05752 * ess;  // derivative of hlat
        //     Compute the thermal conductivity of layers 1 to 3 by function
        //     ThermalCondSoil().
        double rosoil1;  // heat conductivity of 1st soil layer in cal / (cm sec
                         // deg).
        double rosoil2;  // heat conductivity of 2nd soil layer in cal / (cm sec
                         // deg).
        double rosoil3;  // heat conductivity of 3rd soil layer in cal / (cm sec
                         // deg).
        rosoil1 = ThermalCondSoil(VolWaterContent[0][k], so, 1);
        rosoil2 = ThermalCondSoil(VolWaterContent[1][k], so2, 2);
        rosoil3 = ThermalCondSoil(VolWaterContent[2][k], so3, 3);
        //     Compute average rosoil between layers 1 to 3,and heat transfer
        //     from
        //  soil surface to 3rd soil layer.
        double rosoil;  // multiplier for heat flux between 1st and 3rd soil
                        // layers.
        rosoil = (rosoil1 * dl[0] + rosoil2 * dl[1] + rosoil3 * dl[2]) /
                 (dl[0] + dl[1] + dl[2]) / (.5 * dl[0] + dl[1] + .5 * dl[2]);
        //     bbsoil is the heat energy transfer by conductance from soil
        //     surface to soil
        double bbsoil = rosoil * (so - so3);
        //     emtlw is emitted long wave radiation from soil surface
        double emtlw = rls4 * pow(so, 4);
        //     Sensible heat transfer and its derivative
        double dsenheat;  // derivative of senheat
        double senheat;   // sensible heat transfer from soil surface
        if (tm > 0) {
            senheat = hsg * (so - tm);
            dsenheat = hsg;
        } else {
            double tafk;  // average air temperature above soil surface (K)
            tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv);
            senheat = hsg * (so - tafk);
            dsenheat = hsg * (1 - sf * 0.1);
        }
        //     Compute the energy balance bb. (positive direction is upward)
        double bb =
            emtlw     // (a) long wave radiation emitted from soil surface
            - rls1    // long wave radiation reaching the soil surface
            + bbsoil  // (b) heat transfer from soil surface to next soil layer
            + hlat    // (c) latent heat transfer
            - rss     // global radiation reaching the soil surface
            + senheat;  // (d) heat transfer from soil surface to air
                        //
        if (fabs(bb) < 10e-6) return;  // end computation for so
        //     If bb is not small enough, compute its derivative by so.
        double demtlw;  // The derivative of emitted long wave radiation (emtlw)
        demtlw = 4 * rls4 * pow(so, 3);
        //     Compute derivative of bbsoil
        double sop001 = so + 0.001;  // soil surface temperature plus 0.001
        double rosoil1p;  // heat conductivity of 1st soil layer for so+0.001
        rosoil1p = ThermalCondSoil(VolWaterContent[0][k], sop001, 1);
        double rosoilp;  // rosoil for so+0.001
        rosoilp = (rosoil1p * dl[0] + rosoil2 * dl[1] + rosoil3 * dl[2]) /
                  (dl[0] + dl[1] + dl[2]) / (.5 * dl[0] + dl[1] + .5 * dl[2]);
        double drosoil = (rosoilp - rosoil) / 0.001;     // derivative of rosoil
        double dbbsoil = rosoil + drosoil * (so - so3);  // derivative of bbsoil
        //     The derivative of the energy balance function
        double bbp = demtlw       // (a)
                     + dbbsoil    // (b)
                     + dhlat      // (c)
                     + dsenheat;  // (d)
        //     Correct the upper soil temperature by the ratio of bb to bbp.
        double bbadjust;  // the adjustment of soil surface temperature before
                          // next iteration
        bbadjust = bb / bbp;
        //     If adjustment is small enough, no more iterations are needed.
        if (fabs(bbadjust) < 0.002) return;
        //     If bbadjust is not the same sign as bbex, reduce fluctuations
        if (mon <= 1)
            bbex = 0;
        else if (mon >= 2) {
            if (fabs(bbadjust + bbex) < fabs(bbadjust - bbex)) {
                bbadjust = (bbadjust + bbex) / 2;
                so = (so + soex) / 2;
            }
        }
        //
        if (bbadjust > 10) bbadjust = 10;
        if (bbadjust < -10) bbadjust = -10;
        //
        so = so - bbadjust;
        so2 = so2 + (so - soex) / 2;
        so3 = so3 + (so - soex) / 3;
        soex = so;
        bbex = bbadjust;
        mon++;
    }  // end while loop
       //     If (mon >= 50) send message on error and end simulation.
    fprintf(stderr,
            " Infinite loop in SoilSurfaceBalance(). Abnormal stop!! \n");
    fprintf(stderr, "Daynum, ihr, k = %3d %3d %3d\n", Daynum, ihr, k);
    fprintf(stderr, " so      = %10.g\n", so);
    fprintf(stderr, " so2 = %10.g\n", so2);
    fprintf(stderr, " so3 = %10.g\n", so3);
    bEnd = true;
}
