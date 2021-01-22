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
