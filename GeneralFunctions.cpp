//     File GeneralFunctions.cpp contains the following functions, which are
//     used in
//  several places in the model:
//  Date conversions -
//      LeapYear()
//  Soil functions -
//      psiq()
//      qpsi()
//      wcond()
//      PsiOsmotic()
//  Extracting climate data -
//      GetFromClim()
//
#include "GeneralFunctions.h"

#include <math.h>

#include <boost/algorithm/string.hpp>
#include <iostream>

#include "global.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
/////////////////////////////////////////////////////////////////
int LeapYear(int nYear)
//    This function returns 1 if this is a leap year, or 0 if not.
//    Argument input:
//           nYear = the year (4 digit integer).
//    Return value:
//           nLeap =  1 if this is a leap year, 0 if not.
//
{
    int nLeap = 0;
    if (nYear % 4 == 0) nLeap = 1;
    if (nYear % 100 == 0) nLeap = 0;
    if (nYear % 400 == 0) nLeap = 1;

    return nLeap;
}
//////////////////////////////////////////////////////////////////
double psiq(double q, double qr, double qsat, double alpha, double beta)
//     This function computes soil water matric potential (in bars)
//  for a given value of soil water content, using the Van-Genuchten equation.
//
//     The following arguments are used:
//        alpha, beta  - parameters of the van-genuchten equation.
//        q - soil water content, cm3 cm-3.
//        qr - residual water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//
{
    //      For very low values of water content (near the residual water
    //  content) psiq is -500000 bars, and for saturated or higher water
    //  content psiq is -0.00001 bars.
    if ((q - qr) < 0.00001)
        return -500000;
    else if (q >= qsat)
        return -0.00001;
    //     The following equation is used (FORTRAN notation):
    //      PSIX = (((QSAT-QR) / (Q-QR))**(1/GAMA) - 1) **(1/BETA) / ALPHA
    double gama = 1 - 1 / beta;
    double gaminv = 1 / gama;
    double term = (qsat - qr) / (q - qr);  //  intermediate variable
    term = pow(term, gaminv);
    double psix = pow((term - 1), (1 / beta)) / alpha;
    if (psix < 0.01) psix = 0.01;
    //      psix (in cm) is converted to bars (negative value).
    psix = (0.01 - psix) * 0.001;
    if (psix < -500000) psix = -500000;
    if (psix > -0.00001) psix = -0.00001;
    return psix;
}
////////////////////////////////////////////////////////////////////
double qpsi(double psi, double qr, double qsat, double alpha, double beta)
//     This function computes soil water content (cm3 cm-3) for
//  a given value of matric potential, using the Van-Genuchten equation.
//
//     The following arguments are used:
//        alpha, beta  - parameters of the van-genuchten equation.
//        psi - soil water matric potential (bars).
//        qr - residual water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//
{
    //     For very high values of PSI, saturated water content is assumed.
    //     For very low values of PSI, air-dry water content is assumed.
    if (psi >= -0.00001)
        return qsat;
    else if (psi <= -500000)
        return qr;
    //     The soil water matric potential is transformed from bars (psi)
    //  to cm in positive value (psix).
    double psix = 1000 * fabs(psi + 0.00001);
    //     The following equation is used (in FORTRAN notation):
    //      QPSI = QR + (QSAT-QR) / (1 + (ALPHA*PSIX)**BETA)**(1-1/BETA)
    double gama = 1 - 1 / beta;
    double term = 1 + pow((alpha * psix), beta);  //  intermediate variable
    double swfun =
        qr + (qsat - qr) / pow(term, gama);  //  computed water content
    if (swfun < (qr + 0.0001)) swfun = qr + 0.0001;
    return swfun;
}  ////////////////////////////////////////////////////////////////////////
double wcond(double q, double qr, double qsat, double beta,
             double SaturatedHydCond, double PoreSpace)
//     This function computes soil water hydraulic conductivity
//  for a given value of soil water content, using the Van-Genuchten
//  equation. The units of the computed conductivity are the same as the given
//  saturated conductivity (SaturatedHydCond).
//
//     The following arguments are used:
//        beta  - parameter of the van-genuchten equation.
//        SaturatedHydCond - saturated hydraulic conductivity (at qsat).
//        PoreSpace - pore space volume.
//        q - soil water content, cm3 cm-3.
//        qr - residual water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//
{
    //
    //     For very low values of water content (near the residual water
    //  content) wcond is 0.
    if ((q - qr) < 0.0001) return 0;
    //     Water content for saturated conductivity is minimum of PoreSpace and
    //     qsat. For very high values of water content (exceeding the saturated
    //  water content or pore space) conductivity is SaturatedHydCond.
    double xsat = min(qsat, PoreSpace);
    if (q >= xsat) return SaturatedHydCond;
    //      The following equation is used (in FORTRAN notation):
    //      WCOND = CONDSAT * ((Q-QR)/(XSAT-QR))**0.5
    //             * (1-(1-((Q-QR)/(XSAT-QR))**(1/GAMA))**GAMA)**2
    double gama = 1 - 1 / beta;
    double gaminv = 1 / gama;
    double sweff =
        (q - qr) /
        (xsat - qr);  // intermediate variable (effective water content).
    double acoeff =
        pow((1 - pow(sweff, gaminv)), gama);  // intermediate variable
    double bcoeff = pow((1 - acoeff), 2);     // intermediate variable
    double conductivity = pow(sweff, 0.5) * bcoeff * SaturatedHydCond;
    return conductivity;
}
///////////////////////////////////////////////////////////////////////////////
double PsiOsmotic(double q, double qsat, double ec)
//      This function computes soil water osmotic potential (in bars, positive
//      value).
//
//     The following arguments are used:
//        q - soil water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//        ec - electrical conductivity of saturated extract (mmho/cm)
//
{
    double ReturnValue;
    if (ec > 0) {
        ReturnValue = 0.36 * ec * qsat / q;
        if (ReturnValue > 6) ReturnValue = 6;
        return ReturnValue;
    } else
        return 0;
}
///////////////////////////////////////////////////////////////////////////////
double GetFromClim(std::string item, int Doy)
//     This function extracts daily climate values for day of year Doy
//  from the structure Clim.
//     Input arguments:
//       item - string defining which item to extract.
//       Doy -  defines day of year to extract.
//
{
    int i;
    for (i = 0; i < 400; i++) {
        if (Clim[i].nDay == Doy) break;
    }
    //
    if (i > 399) i = 399;
    if (Doy < Clim[0].nDay) i = 0;
    if (item == "tmin")
        return Clim[i].Tmin;
    else if (item == "tmax")
        return Clim[i].Tmax;
    else if (item == "rad")
        return Clim[i].Rad;
    else if (item == "rain")
        return Clim[i].Rain;
    else if (item == "wind")
        return Clim[i].Wind;
    else if (item == "tdew")
        return Clim[i].Tdew;
    else
        return -99;
}
