//  GettingInput_1.cpp
//
//   functions in this file:
// ReadInput()
// ReadProfileFile()
// WriteInitialInputData()
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
    ReadSoilImpedance();
    WriteInitialInputData();
    InitSoil();
    ReadPlantMapInput();
    InitializeSoilData();
    InitializeSoilTemperature();
    InitializeRootData();
    //     initialize some variables at the start of simulation.
    SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
    PlantWeightAtStart =
        TotalRootWeight + TotalStemWeight + TotalLeafWeight + ReserveC;
}
/////////////////////////////////////////////////////////////////////////////
void ReadProfileFile()
//     This function opens and reads the profile file. It is called from
//     ReadInput().
//  It calls GetLineData(), DateToDoy() and OpenOutputFiles().
//     The following global or file-scope variables are set here:
//  AgrInputFileName, bEnd, bLat, bLong, CO2EnrichmentFactor,
//  DayEmerge, DayEndCO2, DayFinish, DayPlant, DayStart, DayStartCO2,
//  DayStartPlantMaps, DayStartSoilMaps, DayStopPlantMaps, DayStopSoilMaps,
//  Elevation, isw, iyear, Latitude, Longitude, m_mulchdata, MulchIndicator,
//  OutIndex, PlantmapFileName, PlantMapFreq, PlantsPerM,
//  RowSpace, SkipRowWidth, SoilHydFileName, SoilInitFileName,
//  SoilMapFreq.
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
    //CO2 enrichment).
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
    std::string SoilMapStartDate, SoilMapStopDate, PlantMapStartDate,
        PlantMapStopDate;
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    if (nLength > 9) SoilMapFreq = stoi(Dummy.substr(0, 10));
    if (nLength >= 16) {
        SoilMapStartDate = Dummy.substr(14, 11);
        boost::algorithm::erase_all(SoilMapStartDate, " ");
    }
    if (nLength >= 31) {
        SoilMapStopDate = Dummy.substr(29, 11);
        boost::algorithm::erase_all(SoilMapStopDate, " ");
    }
    if (nLength > 41) PlantMapFreq = stoi(Dummy.substr(40, 10));
    if (nLength >= 56) {
        PlantMapStartDate = Dummy.substr(54, 11);
        boost::algorithm::erase_all(PlantMapStartDate, " ");
    }
    if (nLength >= 71) {
        PlantMapStopDate = Dummy.substr(69, 11);
        boost::algorithm::erase_all(PlantMapStopDate, " ");
    }
    //     Line #8: 23 output flags.
    // - Line 8 consists of zeros and ones.  "1" tells the simulator to
    //   produce a particular report and "0" indicates that no report
    //   should be produced.
    for (int n = 0; n < 24; n++) OutIndex[n] = 0;
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    for (int n = 0; n < 23 && 3 * n < nLength; n++) {
        int n1 = 3 * n;
        OutIndex[n + 1] = stoi(Dummy.substr(n1, 3));
    }
    DataFile.close();
    //     Calendar dates of emergence, planting, start and stop of simulation,
    //     start and stop of
    // output of soil slab and plant maps are converted to DOY dates by calling
    // function DateToDoy.
    DayStartSoilMaps = DateToDoy(SoilMapStartDate, iyear);
    DayStopSoilMaps = DateToDoy(SoilMapStopDate, iyear);
    DayStartPlantMaps = DateToDoy(PlantMapStartDate, iyear);
    DayStopPlantMaps = DateToDoy(PlantMapStopDate, iyear);
    //     If the output frequency indicators are zero, they are set to 999.
    if (SoilMapFreq <= 0) SoilMapFreq = 999;
    if (PlantMapFreq <= 0) PlantMapFreq = 999;
    //     Call function OpenOutputFiles() to open the output files.
    OpenOutputFiles(m_fileDesc);
}
//////////////////////////////////////////////////////////
void WriteInitialInputData()
//     This function writes the input data to File20 (*.B01). It is executed
//     once
//  at the beginning of the simulation. It is called by ReadInput().
//
//     The following global or file-scope variables are set here:
//  DayEndMulch, DayStartMulch, MulchTranLW, MulchTranSW.
//     The following global or file-scope variables are referenced here:
//  AgrInputFileName, CO2EnrichmentFactor, DayEmerge, DayEndCO2,
//  DayFinish, DayPlant, DayStart, DayStartCO2, Elevation, iyear, Latitude,
//  Longitude, m_mulchdata, maxk, MulchIndicator, OutIndex, PerPlantArea,
//  PlantmapFileName, PlantsPerM, ProfileName, RowSpace,
//  SiteName, SkipRowWidth, SoilHydFileName, SoilInitFileName, VarName.
{
    ofstream File20("Output\\" + ProfileName + ".B01", ios::app);
    File20 << "    Latitude.. ";
    File20.setf(ios::fixed);
    File20.width(8);
    File20.precision(2);
    File20 << fabs(Latitude);
    File20.width(7);
        File20 << "  North";
    File20 << "         Longitude.. ";
    File20.width(8);
    File20.precision(2);
    File20 << fabs(Longitude);
    File20.width(7);
        File20 << "   East";
    File20 << endl;
    //
    if (OutIndex[1] == 0)  // write profile data in metric units
    {
        File20 << "    Elevation (m).... ";
        File20.width(8);
        File20.precision(2);
        File20 << Elevation;
    } else {
        File20 << "    Elevation (ft).... ";
        File20.width(8);
        File20.precision(2);
        File20 << Elevation * 3.28;
    }
    File20 << endl;
    //
    File20 << "    Start Simulation " << DoyToDate(DayStart, iyear);
    File20 << "       Stop Simulation.... " << DoyToDate(DayFinish, iyear)
           << endl;
    //
    File20 << "    Planting date...." << DoyToDate(DayPlant, iyear);
    if (DayEmerge <= 0)
        File20 << "       Emergence date is simulated   " << endl;
    else
        File20 << "       Emergence date...   " << DoyToDate(DayEmerge, iyear)
               << endl;
    //
    if (OutIndex[1] == 0)  // write profile data in metric units
    {
        File20 << "    Row Spacing (cm) ";
        File20.width(8);
        File20.precision(2);
        File20 << RowSpace;
        File20 << "          Plants Per Row-m... ";
        File20.width(8);
        File20 << 0 << endl;
        File20 << "    Skip Width (cm). ";
        File20.width(8);
        File20 << 0;
        File20 << "          Plants Per Ha...... ";
        File20.width(8);
        File20.precision(1);
        File20 << PlantPopulation << endl;
    } else {
        File20 << "    Row Spacing (in) ";
        File20.width(8);
        File20.precision(2);
        File20 << RowSpace / 2.54;
        File20 << "          Plants Per Row-ft.. ";
        File20.width(8);
        File20 << 0 * 0.305 << endl;
        File20 << "    Skip Width (in). ";
        File20.width(8);
        File20 << 0 / 2.54;
        File20 << "          Plants Per Acre.... ";
        File20.width(8);
        File20.precision(1);
        File20 << PlantPopulation * 0.405 << endl;
    }
    //
    if (CO2EnrichmentFactor > 1)  //  write CO2 enrichment input data.
    {
        File20 << "          CO2 enrichment factor              ...  ";
        File20.width(8);
        File20.precision(4);
        File20 << CO2EnrichmentFactor << endl;
        File20 << "    from ...... " << DoyToDate(DayStartCO2, iyear);
        File20 << "    to ...... " << DoyToDate(DayEndCO2, iyear) << endl;
    }
    //
    if (false && MulchIndicator > 0) {
        File20 << "   Polyethylene mulch cover. Transmissivity values are: "
               << endl;
        File20 << " For short waves:  ";
        File20.width(8);
        File20.precision(3);
        File20 << MulchTranSW;
        File20 << " For long waves:  ";
        File20.width(8);
        File20.precision(3);
        File20 << MulchTranLW << endl;
        File20 << " From Day of Year  ";
        File20.width(4);
        File20 << DayStartMulch;
        File20 << " to Day of Year  ";
        File20.width(4);
        File20 << DayEndMulch << endl;
        if (MulchIndicator == 1)
            File20 << " All soil surface covered by mulch." << endl;
        else {
            File20.width(6);
            File20.precision(2);
            if (MulchIndicator == 2)
                File20 << RowSpace / maxk;
            else if (MulchIndicator == 3)
                File20 << 2 * RowSpace / maxk;
            File20 << " cm on each side of planr rows not covered by mulch."
                   << endl;
        }
    }
    File20 << "    Cultural Input File:           " << AgrInputFileName << endl;
    File20 << "    Initial Soil Data Input File:  " << SoilInitFileName << endl;
    File20 << "    Soil Hydrology Input File:     " << SoilHydFileName << endl;
    if (PlantmapFileName.length() > 0)
        File20 << "    Plant Map Adjustment File:     " << PlantmapFileName
               << endl
               << endl;
    //   Write names of the site and the variety
    File20 << "    Site...     " << endl;
    File20 << "    Variety...  " << endl << endl;
}