//  File GettingInput_3.cpp
//    Functions in this file:
// ReadPlantMapInput()
//
#include <boost/algorithm/string.hpp>
#include <iostream>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadPlantMapInput()
//     This sunbroutine opens and reads an ascii file with input of
//  observed plant map adjustment data. It is used to adjust the simulation.
//     It is called by ReadInput().
//     The following global variables are referenced:    iyear, PlantmapFileName
//     The following global variables are set:
//       MapDataGreenBollNum[], MapDataDate[], MapDataMainStemNodes[],
//       MapDataPlantHeight[], MapDataSquareNum[], MapDataAllSiteNum[];
//
{
    if (PlantmapFileName.length() <= 0) return;
    std::string Exten = ".MAP";
    //     Check file extension
    if (PlantmapFileName.substr(PlantmapFileName.length() - 4) != Exten) {
        int strlen = PlantmapFileName.length();
        int newlen = strlen;
        for (int ii = 4; ii > 0; ii--) {
            char Dummy = PlantmapFileName[strlen - ii];
            if (Dummy == '.') {
                newlen = strlen - ii;
                break;
            }
        }
        PlantmapFileName = PlantmapFileName.substr(0, newlen) + Exten;
    }
    std::string m_FilePath = "PLANTMAP\\" + PlantmapFileName;
    //
    ifstream DataFile(m_FilePath, ios::in);
    if (DataFile.fail()) {
        DataFile.close();
        return;
    }
    //     Line #1: Read file description.
    std::string Dummy = GetLineData(DataFile);
    std::string m_PmapDesc;  // Description of the Profile file
    if (Dummy.length() > 20) {
        m_PmapDesc = Dummy.substr(20, 55);
        boost::algorithm::trim_right(m_PmapDesc);
    } else
        m_PmapDesc = "";
    //     Read other lines -
    int i = 0;
    while (DataFile.eof() == 0) {
        std::string StrTemp = GetLineData(DataFile);
        if (StrTemp.length() <= 0) break;
        //
        MapDataDate[i] =
            DateToDoy(StrTemp.substr(0, 11), iyear);  // day of year
        MapDataPlantHeight[i] =
            (double)stof(StrTemp.substr(11, 9));  // Plant height, cm
        MapDataMainStemNodes[i] =
            (double)stof(StrTemp.substr(20, 10));  // Number of mainstem nodes
        MapDataSquareNum[i] = (double)stof(
            StrTemp.substr(30, 10));  // Number of squares per plant
        MapDataGreenBollNum[i] = (double)stof(
            StrTemp.substr(40, 10));  // Number of green bolls per plant
        if (StrTemp.length() >= 61)  // old style *.MAP file from older versions
            MapDataAllSiteNum[i] = (double)stof(
                StrTemp.substr(60, 10));  // Number of total sites per plant
        else
            MapDataAllSiteNum[i] = (double)stof(
                StrTemp.substr(50, 10));  // Number of total sites per plant
                                          //
        i++;
        if (i >= 20) return;
    }
}
