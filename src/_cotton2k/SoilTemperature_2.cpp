//  File SoilTemperature_2.cpp
//
//   List of functions in this file:
//       SensibleHeatTransfer()
//
#include <cmath>
#include "Simulation.hpp"
#include "global.h"
#include "exceptions.h"

using namespace std;

// SoilTemperature_3

/////////////////////////////////////////////////////////////////////////////
double SensibleHeatTransfer(double tsf, double tenviron, double height, double wndcanp)
//     This function computes the sensible heat transfer coefficient, using the friction
//  potential (shear) temperature (thstar), and the surface friction (shear) velocity (ustar)
//  at the atmospheric boundary. It is called three times from EnergyBalance(): for canopy or
//  soil surface with their environment.
//
//     The following arguments are referenced in this function:
//       tenviron - temperature (K) of the environment - air at 200 cm height
//                  for columns with no canopy, or tafk when canopy is present .
//       tsf - surface temperature, K, of soil or canopy.
//       wndcanp - wind speed in the canopy (if present), cm s-1.
//       height - canopy height, cm, or zero for soil surface.
//     The return value computed in this function:
//       sensibleHeatTransfer - raw sensible heat transfer coefficient
{
    //	    Constant values used:
    const double grav = 980;      // acceleration due to gravity (980 cm sec-2).
    const double s40 = 0.13;      // calibration constant.
    const double s42 = 0.63;      // calibration constant.
    const double stmin = 5;       // minimal value of ustar.
    const double vonkar = 0.40;   // Von-Karman constant (0.40).
    const double zalit1 = 0.0962; // parameter .....
                                  //     Wind velocity not allowed to be less than 100 cm s-1.
    double u = wndcanp;           // wind speed at 200 cm height, cm / s.
    if (u < 100)
        u = 100;
    //     Assign initial values to z0 and gtop, and set dt.
    double z0 = s40 * height; // surface roughness parameter, cm.
    if (z0 < 1)
        z0 = 1;
    double gtop = log(
        (200 - s42 * height) / z0); // logarithm of ratio of height of measurement to surface roughness parameter.
    double dt = tsf - tenviron;     // temperature difference.
                                    //     Set approximate initial values for ustar and thstar (to reduce iterations).
    double thstar;                  // friction potential (shear) temperature.
    double ustar;                   // Surface friction (shear) velocity (cm sec-1).
    if (dt >= 0)
    {
        ustar = 1.873 + 0.570172 * dt + .07438568 * u;
        thstar = -0.05573 * dt + 1.987 / u - 6.657 * dt / u;
    }
    else
    {
        ustar = -4.4017 + 1.067 * dt + 0.25957 * u - 0.001683 * dt * u;
        if (ustar < 5)
            ustar = 5;
        thstar = -0.0096 - 0.1149 * dt + 0.0000377 * u + 0.0002367 * dt * u;
        if (thstar < 0.03)
            thstar = 0.03;
    }
    double tbot1 = tsf; // surface temperature corrected for friction (shear) potential temperature.
    long mtest = 0;     // count of iterations.
    double g1;          // temporary derived variable
    double ug1chk;      // previous value of ug1.
    double ug1;         // ratio of ustar to g1.
    double ug1res;      // residual value of ug1.
                        //     Start iterations.
    do
    {
        double thekz; // previous value of thstar.
        double uchek; // previous value of ustar.
        ug1chk = 0;   // previous value of UG1.
        if (mtest > 0)
        {
            //     Assign values to tbot1, uchek, thekz, and ug1chk.
            tbot1 = tsf + zalit1 * thstar * pow((ustar * z0 / 15), 0.45) / vonkar;
            uchek = ustar;
            thekz = thstar;
            if (g1 != 0)
                ug1chk = ustar / g1;
        }
        //     Compute air temperature at 1 cm,, and compute zl and lstar.
        double zl; // nondimensional height.
        if (fabs(thstar) < 1e-30)
            zl = 0;
        else
        {
            double thetmn = (tenviron + tbot1) * 0.5; // mean temperature (K) of air and surface.
            double lstar = (thetmn * ustar * ustar) / (vonkar * grav * thstar);
            zl = (200 - s42 * height) / lstar;
            //     Added to decrease fluctuations in zl:
            if (zl < -5)
                zl = -5;
            if (zl > 0.5)
                zl = 0.5;
        }
        //     Compute g1u, and g2.
        double g1u, g2; // temporary derived variables.
        if (zl > 0)
        {
            g1u = -4.7 * zl;
            g2 = -6.35135 * zl;
            if (g2 < -1)
                g2 = -1;
        }
        else
        {
            double tmp1 = pow((1 - 15 * zl), 0.25); // intermediate variable.
            g1u = 2 * log((1 + tmp1) / 2) + log((1 + tmp1 * tmp1) / 2) - 2 * atan(tmp1 + 1.5708);
            g2 = 2 * log((1 + sqrt(1 - 9 * zl)) / 2);
        }
        if (g2 > gtop)
            g2 = gtop;
        //     Compute ustar and check for minimum value.
        ustar = vonkar * u / (gtop - g1u);
        if (ustar < stmin)
            ustar = stmin;
        //     Compute g1 and thstar.
        g1 = 0.74 * (gtop - g2) + zalit1 * pow((ustar * z0 / 0.15), 0.45);
        thstar = -dt * vonkar / g1;
        //     If more than 30 iterations, reduce fluctuations
        if (mtest > 30)
        {
            thstar = (thstar + thekz) * 0.5;
            ustar = (ustar + uchek) * 0.5;
        }
        mtest++;
        //     Stop simulation if no convergence after 100 iterations.
        if (mtest > 100)
        {
            string msg = " Infinite loop in SensibleHeatTransfer(). Abnormal stop!! \n";
            char C1[12];
            sprintf(C1, "%10.3g", tenviron);
            msg += " tenviron = " + (string)C1 + "\n";
            sprintf(C1, "%10.3g", tsf);
            msg += " tsf      = " + (string)C1 + "\n";
            sprintf(C1, "%10.3g", height);
            msg += " PlantHeight = " + (string)C1 + "\n";
            sprintf(C1, "%10.3g", u);
            msg += " u = " + (string)C1 + "\n";
            throw SimulationEnd(msg);
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
