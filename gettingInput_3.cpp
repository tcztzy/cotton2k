//  File GettingInput_3.cpp
//    Functions in this file:
// OpenClimateFile()
// ReadClimateData()
// tdewest()
// ReadAgriculturalInput()
// SlabLoc()
// ReadPlantMapInput()
//
#include <filesystem>
#include "global.h"
#include "exceptions.h"
#include "GeneralFunctions.h"
#include "Simulation.hpp"

namespace fs = std::filesystem;

static int ReadClimateData(ifstream &, const int &, ClimateStruct[]);

extern "C"
{
    double tdewest(double, double, double);
    int SlabLoc(int, double);
}



///////////////////////////////////////////////////////////////////////////////
static int OpenClimateFile(const string &ActWthFileName, const string &PrdWthFileName, const int &DayStart, ClimateStruct Clim[400])
//     This function gets the climate data file. It is called by ReadInput(),
// and it calls ReadClimateData().
//  Global variables set:  LastDayWeatherData.
//  Return value: LastDayOfActualWeather.
//
{
    int LastDayOfPredictedWeather = 0;    // last DOY of predicted weather data
    int LastDayOfActualWeather = 0;       // last DOY of actual weather data
//     Open file of predicted weather data.
    if (PrdWthFileName.length() > 0) {
        fs::path strFileName = fs::path("climate") / PrdWthFileName;
//     If file does not exist, display message.
        if (!fs::exists(strFileName))
            throw FileNotExists(strFileName);
        else {
            ifstream DataFile(strFileName, ios::in);
            if (DataFile.fail())
                throw FileNotOpened(strFileName);
            else {
                LastDayOfPredictedWeather = ReadClimateData(DataFile, DayStart, Clim);
                DataFile.close();
            }
        }
    }
//     Open file of actual weather data.
    if (ActWthFileName.length() > 0) {
        fs::path strFileName = fs::path("climate") / ActWthFileName;
//     If file does not exist, display message.
        if (!fs::exists(strFileName))
            throw FileNotExists(strFileName);
        else {
            ifstream DataFile1(strFileName, ios::in);
            if (DataFile1.fail())
                throw FileNotOpened(strFileName);
            else {
                LastDayOfActualWeather = ReadClimateData(DataFile1, DayStart, Clim);
                DataFile1.close();
            }
        }
    }
//     Define the last day of weather data
    LastDayWeatherData = LastDayOfActualWeather;
    if (LastDayWeatherData < LastDayOfPredictedWeather)
        LastDayWeatherData = LastDayOfPredictedWeather;
    return LastDayOfActualWeather;
}

/////////////////////////////////////////////////////
int ReadClimateData(ifstream &DataFile, const int &DayStart, ClimateStruct Clim[400])
//     This function reads the climate data file. It is called by OpenClimateFile(), and it
//  calls GetLineData().
//     It stores the data read from the file in the structure Clim.
//
{
//     Read 1st line
    string Dummy = GetLineData(DataFile);
    int nLength = Dummy.length();
//     The following variables inform if the data in the climate file are in metric or "English" units.
    int iswRad = 0, iswTmp = 0, iswRain = 0, iswDewt = 0, iswWind = 0;
    double AverageWind = 0;
    if (nLength >= 31)
        iswRad = atoi(Dummy.substr(31, 3).c_str());
    if (nLength >= 34)
        iswTmp = atoi(Dummy.substr(34, 3).c_str());
    if (nLength >= 37)
        iswRain = atoi(Dummy.substr(37, 3).c_str());
    if (nLength >= 40)
        iswWind = atoi(Dummy.substr(40, 3).c_str());
    if (nLength >= 43)
        iswDewt = atoi(Dummy.substr(43, 3).c_str());
    if (nLength >= 61)
        AverageWind = atof(Dummy.substr(61, 10).c_str());
//     Read all other lines
    int jdd = 0;    // day of year
    while (DataFile.eof() == 0) {
//     Test for end of file, and get the day of year (jdd)
        Dummy = GetLineData(DataFile);
        if (DataFile.eof())
            break;
        jdd = atoi(Dummy.substr(0, 4).c_str());
        int j = jdd - DayStart;   // days from start of simulation
//     Daily weather data are stored in the structure Clim, derived from
//  struct ClimateStruct {int nDay; float Rad, Tmax, Tmin, Rain, Wind, Tdew; }
        if (j >= 0) {
            Clim[j].nDay = jdd;
//     Get float values
            Clim[j].Rad = atof(Dummy.substr(21, 7).c_str());//  Radiation data
            Clim[j].Tmax = atof(Dummy.substr(28, 7).c_str());//  maximum temperature data
            Clim[j].Tmin = atof(Dummy.substr(35, 7).c_str());//  minimum temperature data
            Clim[j].Rain = atof(Dummy.substr(42, 7).c_str());//  rain data
//     If iswRad is 1, convert solar radiation from MJ/m2 to langleys.
            if (iswRad == 1)
                Clim[j].Rad = Clim[j].Rad * 23.884;
//     If iswTmp is 0, convert maximum and minimum temperatures from F to C.
            if (iswTmp == 0) {
                Clim[j].Tmax = (Clim[j].Tmax - 32) / 1.8;
                Clim[j].Tmin = (Clim[j].Tmin - 32) / 1.8;
            }
//     Check if wind data are available
            double f6 = atof(Dummy.substr(49, 7).c_str());
            if (iswWind < 0 || f6 <= 0)
                Clim[j].Wind = AverageWind;
            else {
                Clim[j].Wind = f6;
//     If iswWind is 0, convert daily windrun from miles to km.
                if (iswWind == 0)
                    Clim[j].Wind = Clim[j].Wind * 1.609;
            }
//     Check if dew point temperature data are available
            if (iswDewt < 0)
                Clim[j].Tdew = -100;
            else {
                Clim[j].Tdew = atof(Dummy.substr(56, 7).c_str());
//     If iswDewt is 0, convert dewpoint temperatures from F to C.
                if (iswTmp == 0)
                    Clim[j].Tdew = (Clim[j].Tdew - 32) / 1.8;
            }
//     If iswRain is 0, convert rainfall from inches to mm.
            if (iswRain == 0)
                Clim[j].Rain = Clim[j].Rain * 25.4;
        }
//     Estimate dewpoint temperature when it is not available:
        if (Clim[j].Tdew <= -100)
            Clim[j].Tdew = tdewest(Clim[j].Tmax, SitePar[5], SitePar[6]);
    }
    DataFile.close();
    return jdd;
}

///////////////////////////////////////////////////////////////////////////////////////////
static void ReadAgriculturalInput(Simulation &sim, const string &AgrInputFileName)
//     This function opens the agricultural inputs data file and reads it.
//  It is called by ReadInput() once at the beginning of the simulation.
//
//     The following common variables are referenced:
//  AgrInputFileName, iyear
//     The following common variables are set:
//  CultivationDate, CultivationDepth, DayFirstDef, DayStartPredIrrig, DayStopPredIrrig,
//  DayWaterTableInput, DefoliationDate, DefoliationMethod, DefoliantAppRate,
//  ElCondSatSoil, Irrig (structure), IrrigationDepth, IrrigMethod, LevelsOfWaterTable,
//  LocationColumnDrip, LocationLayerDrip, MaxIrrigation, MinDaysBetweenIrrig,
//  NFertilizer (structure), NumIrrigations, NumNitApps, NumWaterTableData.
//
{
//     Open the input file
    fs::path m_FilePath = fs::path("aginput") / AgrInputFileName;
    ifstream DataFile(m_FilePath, ios::in);
    if (DataFile.fail())
        throw FileNotOpened(m_FilePath);
//     Line #1: Read file description.
    string Dummy = GetLineData(DataFile);
    string m_AgrInptDesc; // Description of the Profile file
    if (Dummy.length() > 20) {
        m_AgrInptDesc = Dummy.substr(20, 55);
        m_AgrInptDesc.erase(m_AgrInptDesc.find_last_not_of(" \r\n\t\f\v") + 1);
    } else
        m_AgrInptDesc = "";
//
    NumNitApps = 0;
    NumIrrigations = 0;
    NumWaterTableData = 0;
    int icult = 0;     // number of cultivation.
    int idef = 0;      // count of this defoliant application.
//     The next lines are read
    while (DataFile.eof() == 0) {
        int isddph; // vertical placement of DRIP, cm from soil surface.
        int isdhrz; // horizontal placement of DRIP, cm from left edge of soil slab.
        string cdate; // date of this application.
//
        string StrTemp = GetLineData(DataFile);
//     The type of input is defined from the first 5 characters of the line.
        if (StrTemp.substr(0, 5) == "IRRIG") {
//     Structure Irrigation {int day, method, LocationColumnDrip, LocationLayerDrip;
//                           double amount; }   Irrig[150];
            cdate = StrTemp.substr(19, 11);
            sim.irrigation[NumIrrigations].day = DateToDoy(cdate.c_str(), sim.year);     // day of year of this irrigation
            sim.irrigation[NumIrrigations].amount = atof(StrTemp.substr(30, 15).c_str()); // net amount of water applied, mm
            sim.irrigation[NumIrrigations].method = atoi(StrTemp.substr(45, 5).c_str());  // method of irrigation: 1=  2=drip
            isdhrz = atoi(StrTemp.substr(50, 5).c_str());              // horizontal placement cm
            isddph = atoi(StrTemp.substr(55, 5).c_str());              // vertical placement cm
//     If this is a drip irrigation, convert distances to soil
//  layer and column numbers by calling SlabLoc.
            if (sim.irrigation[NumIrrigations].method == 2) {
                sim.irrigation[NumIrrigations].LocationColumnDrip = SlabLoc(isdhrz, sim.row_space);
                sim.irrigation[NumIrrigations].LocationLayerDrip = SlabLoc(isddph, 0);
            }
            NumIrrigations++;
        }
//
        else if (StrTemp.substr(0, 5) == "FERTI") {
//     Structure NFertilizer {int day; double amtamm, amtnit, amtura,
//		                      int mthfrt, ksdr, lsdr; }    NFertilizer[150];
            cdate = StrTemp.substr(19, 11);
            NFertilizer[NumNitApps].day = DateToDoy(cdate.c_str(), sim.year);             // day of year
            NFertilizer[NumNitApps].amtamm = (float) atof(StrTemp.substr(30, 10).c_str()); //ammonium N (kg/ha)
            NFertilizer[NumNitApps].amtnit = (float) atof(StrTemp.substr(40, 10).c_str()); //nitrate N (kg/ha)
            NFertilizer[NumNitApps].amtura = (float) atof(StrTemp.substr(50, 10).c_str()); //urea N (kg/ha)
            NFertilizer[NumNitApps].mthfrt = (int) atoi(StrTemp.substr(60, 5).c_str());    // method of application
            isdhrz = atoi(StrTemp.substr(65, 5).c_str());
            // NOTE: string.substr will check boundary
            if (StrTemp.length() > 70)
                isddph = atoi(StrTemp.substr(70, 5).c_str());
//      If this is a side dressing or drip, convert distances to soil
//  layer and column numbers by calling SlabLoc.
            if (NFertilizer[NumNitApps].mthfrt == 1 || NFertilizer[NumNitApps].mthfrt == 3) {
                NFertilizer[NumNitApps].ksdr = SlabLoc(isdhrz, sim.row_space);
                NFertilizer[NumNitApps].lsdr = SlabLoc(isddph, 0);
            } else {
                NFertilizer[NumNitApps].ksdr = 0; // column of application
                NFertilizer[NumNitApps].lsdr = 0; // layer of application
            }
            NumNitApps++;
        }
//
        else if (StrTemp.substr(0, 5) == "CULTI") {
            cdate = StrTemp.substr(19, 11);
            CultivationDate[icult] = DateToDoy(cdate.c_str(), sim.year);
            CultivationDepth[icult] = (float) atof(StrTemp.substr(30, 5).c_str());
            icult++;
        }
//
        else if (StrTemp.substr(0, 3) == "DEF") {
            cdate = StrTemp.substr(19, 11);
            int pgrmth = atoi(StrTemp.substr(30, 5).c_str()); // code number for method of application.
            DefoliationDate[idef] = DateToDoy(cdate.c_str(), sim.year);
//     If this is input for defoliation prediction, define the rate as -99.9,
//  and pgrmth in this case is the percentage of boll opening for which defoliation
//  will be activated, and cdate is the latest date for defoliation application.
            if (StrTemp.substr(4, 5) == "PREDI") {
                DefoliantAppRate[idef] = -99.9;
            }
//     If this is input for actual defoliation, convert application
//  date to Julian date and store it as defdate(idef). Convert
//  application rate to pints per acre and store it as DEFPPA(IDEF).
//  Store the method of application code number as DEFMTH(IDEF).
            else {
                double rtepgr = (float) atof(StrTemp.substr(40, 10).c_str());  // rate of application
                int pgunit = atoi(StrTemp.substr(50, 5).c_str()) - 1; // code number for rate units used.
                if (rtepgr > 0.01) {
                    if (pgunit == 1)
                        rtepgr = rtepgr * 8;
                    else if (pgunit == 2)
                        rtepgr = rtepgr / 16;
                    else if (pgunit == 4)
                        rtepgr = 1 / rtepgr;
                    else if (pgunit == 5)
                        rtepgr = 8 / rtepgr;
                    DefoliantAppRate[idef] = rtepgr;
                }
            }
            DayFirstDef = DefoliationDate[0];
            DefoliationMethod[idef] = pgrmth;
            idef++;
        }
//
        else if (StrTemp.substr(0, 5) == "WATER") {
            cdate = StrTemp.substr(19, 11);
            LevelsOfWaterTable[NumWaterTableData] = (float) atof(StrTemp.substr(30, 5).c_str());
            ElCondSatSoil[NumWaterTableData] = (float) atof(StrTemp.substr(35, 10).c_str());
            DayWaterTableInput[NumWaterTableData] = DateToDoy(cdate.c_str(), sim.year);
            NumWaterTableData++;
        }
//
        else if (StrTemp.substr(0, 5) == "PREDI") {
            MaxIrrigation = atof(StrTemp.substr(20, 7).c_str());
            IrrigMethod = atoi(StrTemp.substr(28, 2).c_str());
            string datstrir = StrTemp.substr(30, 15);
            datstrir.erase(remove(datstrir.begin(), datstrir.end(), ' '), datstrir.end());
            DayStartPredIrrig = DateToDoy(datstrir.c_str(), sim.year); // date to start the predicted irrigation.
            string datstpir = StrTemp.substr(45, 15);
            datstpir.erase(remove(datstpir.begin(), datstpir.end(), ' '), datstpir.end());
            DayStopPredIrrig = DateToDoy(datstpir.c_str(), sim.year); // date to stop the predicted irrigation.
            MinDaysBetweenIrrig = atoi(StrTemp.substr(60, 5).c_str());
            isdhrz = atoi(StrTemp.substr(65, 5).c_str());
            isddph = atoi(StrTemp.substr(70, 5).c_str());
            IrrigationDepth = (float) atof(StrTemp.substr(75, 5).c_str());
//     If this is a drip irrigation, convert distances (input in cm)
//  to soil layer and column numbers by calling SlabLoc.
            if (IrrigMethod == 2) {
                LocationColumnDrip = SlabLoc(isdhrz, sim.row_space);
                LocationLayerDrip = SlabLoc(isddph, 0);
            }
        }
    }  //  end while DataFile
}
