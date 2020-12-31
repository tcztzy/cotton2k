//  File GettingInput_3.cpp
//    Functions in this file:
// OpenClimateFile()
// ReadClimateData()
// tdewest()
// ReadAgriculturalInput()
// SlabLoc()
// ReadPlantMapInput()
//
#include "global.h"
#include "GeneralFunctions.h"

int ReadClimateData(ifstream &, const int &);

double tdewest(double);

int SlabLoc(int, int);

///////////////////////////////////////////////////////////////////////////////
int OpenClimateFile(const string &ActWthFileName, const string &PrdWthFileName, const int &DayStart)
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
                LastDayOfPredictedWeather = ReadClimateData(DataFile, DayStart);
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
                LastDayOfActualWeather = ReadClimateData(DataFile1, DayStart);
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
int ReadClimateData(ifstream &DataFile, const int &DayStart)
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
//  struct Climstruct {int nDay; float Rad, Tmax, Tmin, Rain, Wind, Tdew; }
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
            Clim[j].Tdew = tdewest(Clim[j].Tmax);
    }
    DataFile.close();
    return jdd;
}

//////////////////////////////////////////////////////////
double tdewest(double maxt)
//     This function estimates the approximate daily average dewpoint temperature when 
//  it is not available. It is called by ReadClimateData().
//     Global variables referenced: SitePar[5] and SitePar[6]
//     Argument used:  maxt = maximum temperature of this day.
//
{
    double esttdew; // the value to estimate.
    if (maxt <= 20)
        esttdew = SitePar[5];
    else if (maxt >= 40)
        esttdew = SitePar[6];
    else
        esttdew = ((40 - maxt) * SitePar[5] + (maxt - 20) * SitePar[6]) / 20;
//
    return esttdew;
}

///////////////////////////////////////////////////////////////////////////////////////////
void ReadAgriculturalInput(const string &ProfileName, const string &AgrInputFileName)
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
//  NFertilizer (structure), NumIrrigations, NumNitApps, NumWaterTableData, 
//  pixday, pixmth, pixppa.
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
    int ipx = 0;       // count of this PIX application.
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
            Irrig[NumIrrigations].day = DateToDoy(cdate, iyear);     // day of year of this irrigation
            Irrig[NumIrrigations].amount = atof(StrTemp.substr(30, 15).c_str()); // net amount of water applied, mm
            Irrig[NumIrrigations].method = atoi(StrTemp.substr(45, 5).c_str());  // method of irrigation: 1=  2=drip
            isdhrz = atoi(StrTemp.substr(50, 5).c_str());              // horizontal placement cm
            isddph = atoi(StrTemp.substr(55, 5).c_str());              // vertical placement cm
//     If this is a drip irrigation, convert distances to soil
//  layer and column numbers by calling SlabLoc.
            if (Irrig[NumIrrigations].method == 2) {
                Irrig[NumIrrigations].LocationColumnDrip = SlabLoc(isdhrz, 1);
                Irrig[NumIrrigations].LocationLayerDrip = SlabLoc(isddph, 2);
            }
            NumIrrigations++;
        }
//
        else if (StrTemp.substr(0, 5) == "FERTI") {
//     Structure NFertilizer {int day; double amtamm, amtnit, amtura,
//		                      int mthfrt, ksdr, lsdr; }    NFertilizer[150];
            cdate = StrTemp.substr(19, 11);
            NFertilizer[NumNitApps].day = DateToDoy(cdate, iyear);             // day of year
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
                NFertilizer[NumNitApps].ksdr = SlabLoc(isdhrz, 1);
                NFertilizer[NumNitApps].lsdr = SlabLoc(isddph, 2);
            } else {
                NFertilizer[NumNitApps].ksdr = 0; // column of application
                NFertilizer[NumNitApps].lsdr = 0; // layer of application
            }
            NumNitApps++;
        }
//
        else if (StrTemp.substr(0, 5) == "CULTI") {
            cdate = StrTemp.substr(19, 11);
            CultivationDate[icult] = DateToDoy(cdate, iyear);
            CultivationDepth[icult] = (float) atof(StrTemp.substr(30, 5).c_str());
            icult++;
        }
//
        else if (StrTemp.substr(0, 3) == "DEF") {
            cdate = StrTemp.substr(19, 11);
            int pgrmth = atoi(StrTemp.substr(30, 5).c_str()); // code number for method of application.
            DefoliationDate[idef] = DateToDoy(cdate, iyear);
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
        else if (StrTemp.substr(0, 5) == "PIX  ") {
            cdate = StrTemp.substr(19, 11);
            pixmth[ipx] = atoi(StrTemp.substr(30, 5).c_str());
            double rtepgr = (float) atof(StrTemp.substr(40, 10).c_str()); // rate of application.
            int pgunit = atoi(StrTemp.substr(50, 5).c_str()) - 1;         // code number for rate units used.
            if (rtepgr > 0.01) {
                pixday[ipx] = DateToDoy(cdate, iyear);
                if (pgunit == 1)
                    rtepgr = rtepgr * 8;
                if (pgunit == 2)
                    rtepgr = rtepgr / 16;
                if (pgunit == 4)
                    rtepgr = 1 / rtepgr;
                if (pgunit == 5)
                    rtepgr = 8 / rtepgr;
                pixppa[ipx] = rtepgr;
            }
            ipx++;
        }
//
        else if (StrTemp.substr(0, 5) == "WATER") {
            cdate = StrTemp.substr(19, 11);
            LevelsOfWaterTable[NumWaterTableData] = (float) atof(StrTemp.substr(30, 5).c_str());
            ElCondSatSoil[NumWaterTableData] = (float) atof(StrTemp.substr(35, 10).c_str());
            DayWaterTableInput[NumWaterTableData] = DateToDoy(cdate, iyear);
            NumWaterTableData++;
        }
//
        else if (StrTemp.substr(0, 5) == "PREDI") {
            MaxIrrigation = atof(StrTemp.substr(20, 7).c_str());
            IrrigMethod = atoi(StrTemp.substr(28, 2).c_str());
            string datstrir = StrTemp.substr(30, 15);
            datstrir.erase(remove(datstrir.begin(), datstrir.end(), ' '), datstrir.end());
            DayStartPredIrrig = DateToDoy(datstrir, iyear); // date to start the predicted irrigation.
            string datstpir = StrTemp.substr(45, 15);
            datstpir.erase(remove(datstpir.begin(), datstpir.end(), ' '), datstpir.end());
            DayStopPredIrrig = DateToDoy(datstpir, iyear); // date to stop the predicted irrigation.
            MinDaysBetweenIrrig = atoi(StrTemp.substr(60, 5).c_str());
            isdhrz = atoi(StrTemp.substr(65, 5).c_str());
            isddph = atoi(StrTemp.substr(70, 5).c_str());
            IrrigationDepth = (float) atof(StrTemp.substr(75, 5).c_str());
//     If this is a drip irrigation, convert distances (input in cm)
//  to soil layer and column numbers by calling SlabLoc.
            if (IrrigMethod == 2) {
                LocationColumnDrip = SlabLoc(isdhrz, 1);
                LocationLayerDrip = SlabLoc(isddph, 2);
            }
        }
    }  //  end while DataFile
//
//     Writing of Agricultural Input Data to file B01:
    ofstream File20(fs::path("output") / (ProfileName + ".B01"), ios::app);
    if (NumIrrigations > 0)
//     Write irrigations set by input:
    {
        File20 << endl << "    Irrigations set by input file: " << endl;
        File20 << "       Date           Amount        Method          Drip Location " << endl;
        if (OutIndex[1] == 0)
            File20 << "                        mm                        Column       Layer " << endl;
        else
            File20 << "                      inches                      Column       Layer " << endl;
        for (int ii = 0; ii < NumIrrigations; ii++) {
            File20.width(14);
            File20 << DoyToDate(Irrig[ii].day, iyear);     // date of this irrigation
            File20.setf(ios::fixed);
            File20.width(14);
            File20.precision(2);
            if (OutIndex[1] == 0)
                File20 << Irrig[ii].amount; // net amount of water applied, mm
            else
                File20 << Irrig[ii].amount / 25.4; // net amount of water applied, inches
//     Methods of irrigation: (0, "SPRINKLER") (1, "FURROW") (2, "DRIP")
            string IrrMethod = "";
            if (Irrig[ii].method == 0)
                IrrMethod = "Sprinkler";
            else if (Irrig[ii].method == 1)
                IrrMethod = "Furrow";
            else if (Irrig[ii].method == 2)
                IrrMethod = "Drip";
            File20.width(14);
            File20 << IrrMethod; // net amount of water applied, mm
            if (Irrig[ii].method == 2) {
                File20.width(12);
                File20 << Irrig[ii].LocationColumnDrip;              // horizontal placement cm
                File20.width(12);
                File20 << Irrig[ii].LocationLayerDrip;              // vertical placement cm
            }
            File20 << endl;
        }
    }
    if (NumNitApps > 0)
//     Write Nitrogen applications set by input:
    {
        File20 << endl << "    Nitrogen applications set by input file: " << endl;
        File20 << "       Date           NH4 N        NO3 N       Urea N        Method       Drip Location " << endl;
        if (OutIndex[1] == 0)
            File20 << "                                  kg / ha                               Column       Layer "
                   << endl;
        else
            File20 << "                                  lbs / acre                            Column       Layer "
                   << endl;
        for (int ii = 0; ii < NumIrrigations; ii++) {
            File20.width(14);
            File20 << DoyToDate(NFertilizer[ii].day, iyear);     // date of this irrigation
            File20.setf(ios::fixed);
            File20.width(13);
            File20.precision(2);
            if (OutIndex[1] == 0) {
                File20 << NFertilizer[ii].amtamm; // net amount of NH4 N, kg/ha
                File20.width(13);
                File20 << NFertilizer[ii].amtnit; // net amount of NO3 N, kg/ha
                File20.width(13);
                File20 << NFertilizer[ii].amtura; // net amount of Urea N, kg/ha
            } else {
                File20 << NFertilizer[ii].amtamm / 1.121F; // net amount of NH4 N, lbs/ac
                File20.width(13);
                File20 << NFertilizer[ii].amtnit / 1.121F; // net amount of NO3 N, lbs/ac
                File20.width(13);
                File20 << NFertilizer[ii].amtura / 1.121F; // net amount of Urea N, lbs/ac
            }
//     Methods of application: (0, "Broadcast") (1, "Side-dress") (2, "Foliar") (3, "DRIP")
            string FMethod = "";
            if (NFertilizer[ii].mthfrt == 0)
                FMethod = "Broadcast";
            else if (NFertilizer[ii].mthfrt == 1)
                FMethod = "Sidedress";
            else if (NFertilizer[ii].mthfrt == 2)
                FMethod = "Foliar";
            else if (NFertilizer[ii].mthfrt == 3)
                FMethod = "Drip";
            File20.width(13);
            File20 << FMethod;
            if (NFertilizer[ii].mthfrt == 1 || NFertilizer[ii].mthfrt == 3) {
                File20.width(12);
                File20 << NFertilizer[ii].ksdr;    // horizontal placement cm
                File20.width(12);
                File20 << NFertilizer[ii].lsdr;     // vertical placement cm
            }
            File20 << endl;
        }
    }
    if (icult > 0)
//     Write cultivations set by input:
    {
        File20 << endl << "    Cultivations set by input file: " << endl;
        if (OutIndex[1] == 0)
            File20 << "       Date         Depth, cm " << endl;
        else
            File20 << "       Date         Depth, inches " << endl;
        for (int ii = 0; ii < icult; ii++) {
            File20.width(14);
            File20 << DoyToDate(CultivationDate[ii], iyear);     // date of this irrigation
            File20.setf(ios::fixed);
            File20.width(14);
            File20.precision(2);
            if (OutIndex[1] == 0)
                File20 << CultivationDepth[ii]; // cm
            else
                File20 << CultivationDepth[ii] / 25.4; //  inches
            File20 << endl;
        }
    }
    if (ipx > 0)
//     Write pix applications set by input:
    {
        File20 << endl << "    Pix applications set by input file: " << endl;
        File20 << "       Date           Method          Rate" << endl;
        for (int ii = 0; ii < ipx; ii++) {
            File20.width(14);
            File20 << DoyToDate(pixday[ii], iyear);     // date of this application
            File20.setf(ios::fixed);
            File20.width(14);
            File20 << pixmth[ii]; //  method code
            File20.width(14);
            File20.precision(2);
            File20 << pixppa[ii]; //  rate of app
            File20 << endl;
        }
    }
    if (idef > 0)
//     Write defoliation applications set by input:
    {
        File20 << endl << "    Defoliation applications set by input file: " << endl;
        File20 << "       Date           Method          Rate" << endl;
        for (int ii = 0; ii < idef; ii++) {
            if (DefoliantAppRate[ii] == -99.9)
                continue;
            File20.width(14);
            File20 << DoyToDate(DefoliationDate[ii], iyear);  // date of this application
            File20.setf(ios::fixed);
            File20.width(14);
            File20 << DefoliationMethod[ii]; //  method code
            File20.width(14);
            File20.precision(2);
            File20 << DefoliantAppRate[ii]; //  rate of app
            File20 << endl;
        }
        for (int ii = 0; ii < idef; ii++) {
            if (DefoliantAppRate[ii] == -99.9) {
//     Write setting defoliation prediction data
                File20 << endl << " Setting prediction of defoliation at percentage boll opening of ";
                File20.width(10);
                File20 << DefoliationMethod[ii] << endl;  // boll opening percentage
                File20 << "         Last date of defoliation prediction:  ";
                File20 << DoyToDate(DefoliationDate[ii], iyear) << endl;  // last date of prediction
            }
        }
    }
    if (NumWaterTableData > 0)
//     Write defoliation applications set by input:
    {
        File20 << endl << "    Water Table and Salinity Data set by input file: " << endl;
        File20 << "       Date       Water Table, cm   Electr-Conductivity" << endl;
        for (int ii = 0; ii < NumWaterTableData; ii++) {
            File20.width(14);
            File20 << DoyToDate(DayWaterTableInput[ii], iyear);  // date of this application
            File20.setf(ios::fixed);
            File20.width(14);
            File20 << LevelsOfWaterTable[ii]; //  Water table level
            File20.width(14);
            File20.precision(2);
            File20 << ElCondSatSoil[ii]; //  salinity
            File20 << endl;
        }
    }
    if (MaxIrrigation > 0)
//     Write data for setting predicted irrigations:
    {
        string IrrMethod = "";
        if (IrrigMethod == 0)
            IrrMethod = "Sprinkler";
        else if (IrrigMethod == 1)
            IrrMethod = "Furrow";
        else if (IrrigMethod == 2)
            IrrMethod = "Drip";
        File20 << endl << "    Setting Prediction of Irrigation. Method: ";
        File20.width(14);
        File20 << IrrMethod << endl;
        File20 << " Starting Date " << DoyToDate(DayStartPredIrrig, iyear);
        File20 << " Stopping Date " << DoyToDate(DayStopPredIrrig, iyear) << endl;
        File20 << " Minimum days between irrigations = ";
        File20.width(14);
        File20 << MinDaysBetweenIrrig << endl;
        if (IrrigMethod == 2) {
            File20 << " Location of Drip Column = ";
            File20.width(8);
            File20 << LocationColumnDrip;
            File20 << " Location of Drip Layer = ";
            File20.width(8);
            File20 << LocationLayerDrip << endl;
        } else {
            File20 << " Depth of Irrigation, cm = ";
            File20.width(12);
            File20.precision(2);
            File20 << IrrigationDepth << endl;
        }
        File20 << endl;
    }
}

/////////////////////////////////////////////////////////////////
int SlabLoc(int isd, int index)
//     This function computes the layer (lsdr) or column (ksdr) where the emitter 
//  of drip irrigation, or the fertilizer side - dressing is located. It is called
//  from ReadAgriculturalInput().
//     The following input arguments are used:
//               isd = horizontal or vertical distance
//               index = 1 if horizontal or = 2 if vertical 
//     The following global variables are referenced here:
//                dl, nk, nl, wk.
//
{
    int result = 0; // the resulting column or layer number.
    if (index == 1)  // horizontal
    {
//     Define the column of this location
        double sumwk = 0; // sum of soil column widths.
        for (int k = 0; k < nk; k++) {
            sumwk += wk[k];
            if (sumwk >= isd) {
                result = k;
                break;
            }
        }
    } else if (index == 2)  // vertical
    {
//	    Define the layer of this location
        double sumdl = 0; // sum of soil layer depths.
        for (int l = 0; l < nl; l++) {
            sumdl += dl[l];
            if (sumdl >= isd) {
                result = l;
                break;
            }
        }
    }
    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadPlantMapInput(const string &PlantmapFileName)
//     This sunbroutine opens and reads an ascii file with input of
//  observed plant map adjustment data. It is used to adjust the simulation.
//     It is called by ReadInput().
//     The following global variables are referenced:    iyear, 
//     The following global variables are set:
//       MapDataGreenBollNum[], MapDataDate[], MapDataMainStemNodes[], 
//       MapDataPlantHeight[], MapDataSquareNum[], MapDataAllSiteNum[];
//
{
    if (PlantmapFileName.length() <= 0)
        return;
    string Exten = ".MAP";
//     Check file extension
    string plantMapFileName = PlantmapFileName;
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
        // NOTE: Even PlantmapFileName had been changed here, it never be used then
        plantMapFileName = PlantmapFileName.substr(0, newlen) + Exten;
    }
    fs::path m_FilePath = fs::path("plantmap") / plantMapFileName;
//
    ifstream DataFile(m_FilePath, ios::in);
    if (DataFile.fail())
        throw FileNotOpened(m_FilePath);
//     Line #1: Read file description.
    string Dummy = GetLineData(DataFile);
    string m_PmapDesc; // Description of the Profile file
    if (Dummy.length() > 20) {
        m_PmapDesc = Dummy.substr(20, 55);
        m_PmapDesc.erase(m_PmapDesc.find_last_not_of(" \r\n\t\f\v") + 1);
    } else
        m_PmapDesc = "";
//     Read other lines -
    int i = 0;
    while (DataFile.eof() == 0) {
        string StrTemp = GetLineData(DataFile);
        if (StrTemp.length() <= 0)
            break;
//
        MapDataDate[i] = DateToDoy(StrTemp.substr(0, 11), iyear);   // day of year
        MapDataPlantHeight[i] = (double) atof(StrTemp.substr(11, 9).c_str());  // Plant height, cm
        MapDataMainStemNodes[i] = (double) atof(StrTemp.substr(20, 10).c_str()); // Number of mainstem nodes
        MapDataSquareNum[i] = (double) atof(StrTemp.substr(30, 10).c_str()); // Number of squares per plant
        MapDataGreenBollNum[i] = (double) atof(StrTemp.substr(40, 10).c_str()); // Number of green bolls per plant
        if (StrTemp.length() >= 61)  // old style *.MAP file from older versions
            MapDataAllSiteNum[i] = (double) atof(StrTemp.substr(60, 10).c_str()); // Number of total sites per plant
        else
            MapDataAllSiteNum[i] = (double) atof(StrTemp.substr(50, 10).c_str()); // Number of total sites per plant
//
        i++;
        if (i >= 20)
            return;
    }
}
