//  GettingInput_1.cpp
//
//   functions in this file:
// ReadInput()
// ReadProfileFile()
//
#include <math.h>

#include <boost/algorithm/string.hpp>
#include <iostream>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

/////////////////////////////////////////////////////////////
void ReadInput()
//     This is the main function for reading input. It is called from
//     RunTheModel(). The following global variables are set here:
//        PlantWeightAtStart , SoilNitrogenAtStart
//     The following global variables are referenced here:
//        ReserveC, TotalLeafWeight, TotalRootWeight, TotalSoilNh4N,
//        TotalSoilNo3N, TotalSoilUreaN, TotalStemWeight.
//
{
    //     The following functions are called to read initial values of some
    //     variables from
    //  input files, or initialize them otherwise.
    ReadProfileFile();
    ReadPlantMapInput();
    InitializeRootData();
    //     initialize some variables at the start of simulation.
    SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
    PlantWeightAtStart =
        TotalRootWeight + TotalStemWeight + TotalLeafWeight() + ReserveC;
}
/////////////////////////////////////////////////////////////////////////////
void ReadProfileFile()
//     This function opens and reads the profile file. It is called from
//     ReadInput().
//  It calls GetLineData(), DateToDoy() and OpenOutputFiles().
//     The following global or file-scope variables are set here:
//  AgrInputFileName, bEnd, bLat, bLong, CO2EnrichmentFactor,
//  DayEmerge, DayEndCO2, DayFinish, DayPlant, DayStart, DayStartCO2,
//  Elevation, isw, iyear, Latitude, Longitude, m_mulchdata, MulchIndicator,
//  PlantmapFileName, PlantsPerM,
//  RowSpace, SkipRowWidth, SoilHydFileName, SoilInitFileName.
//     The following global variable is referenced here:   ProfileName
//
{
    std::string strFileName =
        "PROFILES\\" + ProfileName + ".PRO";  // file name with path
    ifstream DataFile(strFileName, ios::in);
    if (DataFile.fail()) {
        std::cerr << "Error opening " + strFileName + "." << std::endl;
        DataFile.close();
        return;
    }
    //     Line #1: Read file description.
    std::string Dummy = GetLineData(DataFile);
    std::string m_fileDesc;  // Description of the Profile file
    //     Line #2: Read dates of emergence, start and end of simulation, and
    //     planting date.
    std::string DateEmerge, DateSimStart, DateSimEnd, DatePlant;
    Dummy = GetLineData(DataFile);
    int nLength = Dummy.length();
    //     For advanced users only: if there is CO2 enrichment, read also CO2
    //     factor, DOY dates
    //	for start and stop of enrichment (these are left blank if there is no
    // CO2 enrichment).
    //     Line #3: Names of weather files: actual and predicted.
    Dummy = GetLineData(DataFile);
    //     For advanced users only: If soil mulch is used, read relevant
    //     parameters.
    //     Line #4: Names of files for soil hydraulic data, soil initial
    //  conditions, agricultural input, and plant map adjustment.
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    if (nLength > 1) {
        SoilHydFileName = Dummy.substr(0, 20);
        boost::algorithm::erase_all(SoilHydFileName, " ");
    }
    if (nLength > 20) {
        SoilInitFileName = Dummy.substr(20, 20);
        boost::algorithm::erase_all(SoilInitFileName, " ");
    }
    if (nLength > 40) {
        AgrInputFileName = Dummy.substr(40, 20);
        boost::algorithm::erase_all(AgrInputFileName, " ");
    }
    if (nLength > 60) {
        PlantmapFileName = Dummy.substr(60);
        boost::algorithm::erase_all(PlantmapFileName, " ");
    } else
        PlantmapFileName = "";
    //     Line #5: Latitude and longitude of this site, elevation (in m
    //  above sea level), and the index number for this geographic site.
    Dummy = GetLineData(DataFile);
    //     Line #6: Row spacing in cm, skip-row spacing in cm (blank or 0
    //  for no skip rows), number of plants per meter of row, and index
    //  number for the cultivar.
    Dummy = GetLineData(DataFile);
    //     Line #7: Frequency in days for output of soil maps, and dates
    //  for start and stop of this output (blank or 0 if no such output is
    //  required. Same is repeated for output of plant maps.
    Dummy = GetLineData(DataFile);
    //     Line #8: 23 output flags.
    // - Line 8 consists of zeros and ones.  "1" tells the simulator to
    //   produce a particular report and "0" indicates that no report
    //   should be produced.
    Dummy = GetLineData(DataFile);
    DataFile.close();
    //     Calendar dates of emergence, planting, start and stop of simulation,
    //     start and stop of
    // output of soil slab and plant maps are converted to DOY dates by calling
    // function DateToDoy.
    //     If the output frequency indicators are zero, they are set to 999.
}