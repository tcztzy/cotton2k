//     File GeneralFunctions.cpp contains the following functions, which are used in
//  several places in the model:
//  String handling -
//      GetLineData()
//  Date conversions - 
//      DateToDoy()
//      DoyToDate()
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

///////////////////////////////////////////////////////////////////////
string GetLineData(ifstream &DataFile)
//     This function safely reads a line from input. All data are read into a temporary character
//  string with the "getline" command. It stops at end of line, and converts it to string.
//     Input        : the ifstream file which is being read.
//     return value : string with the trailing blanks removed.
//
{
    char tmp[91];   // Temporarily stores a character array.
    DataFile.getline(tmp, 90);
    string LineString = tmp;
    LineString.erase(LineString.find_last_not_of(" \r\n\t\f\v") + 1);
    return LineString;
}

/////////////////////////////////////////////////////////////////////////
// Date conversion functions:
//
string DoyToDate(int Doy, int m_YearStart)
//     This function converts day of year (sometimes called 'Julian date')
// to calendar date, allowing for leap years and for days in the following year.
//     Arguments input:
//            Doy = Day of year.
//            m_YearStart = year number (4 digit) at start of simulation.
//     Return value:
//         Calender date as an 11-character string 'dd-MON-yyyy'
//
//     NOTE: This function does not check for correct input. Make sure that Doy 
// is an integer between 1 and 730 (or 731 if m_YearStart or m_YearStart+1 is a leap year),
// and m_YearStart is a 4-digit number).
//
{
    static string MonthName[] =
            {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
    static int mday[] =
            {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (Doy <= 0)
        return "           ";

//  Adjust number of days in February for leap years.
    mday[1] = 28 + LeapYear(m_YearStart);

    int iday = Doy;
    int iy = m_YearStart;
//  If this is the following year
    int jadd = 365 + LeapYear(m_YearStart);
    if (Doy > jadd) {
        iday = Doy - jadd;
        iy++;
        mday[1] = 28 + LeapYear(iy);
    }
//   Compute month and day from Julian date.
    int month = 1;
    for (int i = 0; i < 12; i++) {
        if (iday <= mday[i])
            break;
        month++;
        iday = iday - mday[i];
    }
//   Construct date string.
    string Date1 = "  ";
    char Cday[3];
    sprintf(Cday, "%2d", iday);
    if (Cday[0] == ' ')
        Cday[0] = '0';
    Date1.replace(0, 1, 1, Cday[0]);
    Date1.replace(1, 1, 1, Cday[1]);

    string Date2 = "    ";
    char Cyear[5];
    sprintf(Cyear, "%4d", iy);
    Date2.replace(0, 1, 1, Cyear[0]);
    Date2.replace(1, 1, 1, Cyear[1]);
    Date2.replace(2, 1, 1, Cyear[2]);
    Date2.replace(3, 1, 1, Cyear[3]);

    string DateOut;
    DateOut = Date1 + "-" + MonthName[month - 1] + "-" + Date2;
    return DateOut;
}

///////////////////////////////////////////////////////////////////////////////
double GetFromClim(const Climstruct Clim[400], const string& item, const int& Doy)
//     This function extracts daily climate values for day of year Doy
//  from the structure Clim.
//     Input arguments:
//       item - string defining which item to extract.
//       Doy -  defines day of year to extract.
//
{
    int i;
    for (i = 0; i < 400; i++) {
        if (Clim[i].nDay == Doy)
            break;
    }
//
    if (i > 399)
        i = 399;
    if (Doy < Clim[0].nDay)
        i = 0;
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
