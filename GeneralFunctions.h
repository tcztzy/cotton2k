//  General auxiliary functions
//
#include <fstream>  // Necessary for file I/O
using namespace std;
//
//  definition of functions
//  =======================
    CString GetLineData(ifstream &DataFile);
// Date conversion functions:
    int DateToDoy(CString Date, int m_YearStart);
    CString DoyToDate(int Doy, int m_YearStart);
    int LeapYear(int nYear);
//  soil water functions:
   double qpsi ( double psi, double qr, double qsat, double alpha, double beta );
   double wcond ( double q, double qr, double qsat, double beta, double SaturatedHydCond, double PoreSpace );
   double psiq ( double q, double qr, double qsat, double alpha, double beta );
   double PsiOsmotic ( double q, double qsat, double ec );
   double GetFromClim(CString item, int Doy);
