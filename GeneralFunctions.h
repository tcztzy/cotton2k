//  General auxiliary functions
//
#include <fstream>  // Necessary for file I/O
using namespace std;
//
//  definition of functions
//  =======================
// Date conversion functions:
int LeapYear(int nYear);
//  soil water functions:
double qpsi(double psi, double qr, double qsat, double alpha, double beta);
double wcond(double q, double qr, double qsat, double beta,
             double SaturatedHydCond, double PoreSpace);
double psiq(double q, double qr, double qsat, double alpha, double beta);
double PsiOsmotic(double q, double qsat, double ec);
double GetFromClim(std::string item, int Doy);
inline double drop_leaf_age(double lai) {
    return 140. - 1. * lai;
}
