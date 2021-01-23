//  General auxiliary functions
#include <fstream> // Necessary for file I/O
#include "cotton2k.h"

using namespace std;

string GetLineData(ifstream &DataFile);

extern "C"
{
    int64_t DateToDoy(const char *, int32_t);
    const char* DoyToDate(int Doy, int m_YearStart);
    double PsiOsmotic(double q, double qsat, double ec);
    double qpsi(double psi, double qr, double qsat, double alpha, double beta);
    double psiq(double q, double qr, double qsat, double alpha, double beta);
    double wcond(double q, double qr, double qsat, double beta, double SaturatedHydCond, double PoreSpace);
}
