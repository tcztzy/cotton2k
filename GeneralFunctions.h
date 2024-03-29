//  General auxiliary functions
//
#include "global.h"  // Necessary for file I/O
//
//  definition of functions
//  =======================
//  soil water functions:
double qpsi(double psi, double qr, double qsat, double alpha, double beta);
double wcond(double q, double qr, double qsat, double beta,
             double SaturatedHydCond, double PoreSpace);
double psiq(double q, double qr, double qsat, double alpha, double beta);
double PsiOsmotic(double q, double qsat, double ec);
double GetFromClim(CLIMATE_METRIC item, int Doy);
inline double drop_leaf_age(double lai) {
    return 140. - 1. * lai;
}
