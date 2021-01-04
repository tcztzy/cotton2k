//  General auxiliary functions
//
#include <fstream>  // Necessary for file I/O
#include "cotton2k.h"

using namespace std;

//
//  definition of functions
//  =======================
string GetLineData(ifstream &DataFile);

// Date conversion functions:
int DateToDoy(string Date, int m_YearStart);

string DoyToDate(int Doy, int m_YearStart);

int LeapYear(int nYear);

//  soil water functions:
double qpsi(double psi, double qr, double qsat, double alpha, double beta);

double wcond(double q, double qr, double qsat, double beta, double SaturatedHydCond, double PoreSpace);

double psiq(double q, double qr, double qsat, double alpha, double beta);

double PsiOsmotic(double q, double qsat, double ec);

double GetFromClim(const Climstruct Clim[400], const string& item, const int& Doy);
