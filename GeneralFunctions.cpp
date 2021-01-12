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
