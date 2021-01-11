//  General auxiliary functions
//
#include <fstream> // Necessary for file I/O
#include "cotton2k.h"

using namespace std;

//
//  definition of functions
//  =======================
string GetLineData(ifstream &DataFile);

//  soil water functions:
double GetFromClim(const Climstruct Clim[400], const string &item, const int &Doy);

extern "C"
{
    uint64_t LeapYear(uint64_t);
    int64_t DateToDoy(const char *, int32_t);
    const char* DoyToDate(int Doy, int m_YearStart);
    double PsiOsmotic(double q, double qsat, double ec);
    double qpsi(double psi, double qr, double qsat, double alpha, double beta);
    double psiq(double q, double qr, double qsat, double alpha, double beta);
    double wcond(double q, double qr, double qsat, double beta, double SaturatedHydCond, double PoreSpace);
}
