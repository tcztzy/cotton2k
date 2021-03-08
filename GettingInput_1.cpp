//  GettingInput_1.cpp
//
//   functions in this file:
// ReadInput()
// ReadProfileFile()
// ReadCalibrationData()
// InitializeGrid()
// WriteInitialInputData()
//
#include <filesystem>
#include <vector>
#include "global.h"
#include "exceptions.h"
#include "GeneralFunctions.h"
#include "Output.h"

using namespace std;
namespace fs = std::filesystem;

static Simulation ReadProfileFile(const char *, vector<string> &);

static void ReadCalibrationData();

static void InitializeGrid(Simulation &);

extern "C"
{
    void WriteInitialInputData(Simulation &, bool, double, double, double, const char *, int, const char *, const char *, const char *, const char *, const char *, const char *);
}
// GettingInput_2
static void InitSoil(const string &);

static void ReadSoilImpedance(Simulation &);

static void InitializeSoilData(Simulation &, const string &);

static void InitializeSoilTemperature();

static void InitializeRootData(Simulation &);

// GettingInput_3
static int OpenClimateFile(const string &, const string &, const int &, ClimateStruct[400]);

static void ReadAgriculturalInput(Simulation &, const string &);

//
// Definitions of File scope variables:
static int nSiteNum,               // index number for site.
    nVarNum;                // index number for cultivar.
static double SkipRowWidth,        // the smaller distance between skip rows, cm
    PlantsPerM;             // average number of plants pre meter of row.
static string
    VarName,                // name of the cultivar
    SiteName;               // name of the site

/////////////////////////////////////////////////////////////////////////////
Simulation ReadProfileFile(const char *ProfileName, vector<string> &filenames)
//     This function opens and reads the profile file. It is called from ReadInput().
//  It calls GetLineData(), DateToDoy() and OpenOutputFiles().
//     The following global or file-scope variables are set here:
//  bLat, bLong, CO2EnrichmentFactor,
//  DayEndCO2, DayStartCO2, DayStartPlantMaps, DayStartSoilMaps,
//  DayStopPlantMaps, DayStopSoilMaps, Elevation, isw, iyear, Latitude, Longitude, m_mulchdata,
//  MulchIndicator, nSiteNum, nVarNum, OutIndex, PlantMapFreq, PlantsPerM,
//  RowSpace, SkipRowWidth, SoilHydFileName, SoilInitFileName, SoilMapFreq.
//
{
    fs::path strFileName = fs::path("profiles") / (string(ProfileName) + ".pro"); // file name with path
                                                                                  //     If file does not exist, or can not be opened, display message
    if (!fs::exists(strFileName))
        throw FileNotExists(strFileName);
    ifstream DataFile(strFileName, ios::in);
    if (DataFile.fail())
        throw FileNotOpened(strFileName);
    //     Line #1: Read file description.
    string Dummy = GetLineData(DataFile);
    //     Line #2: Read dates of emergence, start and end of simulation, and planting date.
    string DateEmerge, DateSimStart, DateSimEnd, DatePlant;
    Dummy = GetLineData(DataFile);
    int nLength = Dummy.length();
    if (nLength >= 14)
    {
        DateEmerge = Dummy.substr(0, 14);
        DateEmerge.erase(remove(DateEmerge.begin(), DateEmerge.end(), ' '), DateEmerge.end());
    }
    if (nLength >= 23)
    {
        DateSimStart = Dummy.substr(15, 14);
        DateSimStart.erase(remove(DateSimStart.begin(), DateSimStart.end(), ' '), DateSimStart.end());
    }
    if (nLength >= 38)
    {
        DateSimEnd = Dummy.substr(30, 14);
        DateSimEnd.erase(remove(DateSimEnd.begin(), DateSimEnd.end(), ' '), DateSimEnd.end());
    }
    if (nLength >= 53)
    {
        DatePlant = Dummy.substr(45, 14);
        DatePlant.erase(remove(DatePlant.begin(), DatePlant.end(), ' '), DatePlant.end());
    }
    else
        DatePlant = "";
    int year = atoi(DateSimStart.substr(7, 4).c_str());
    //     For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates
    //	for start and stop of enrichment (these are left blank if there is no CO2 enrichment).
    double CO2EnrichmentFactor = 0;
    uint32_t DayStartCO2 = 0;
    uint32_t DayEndCO2 = 0;
    if (nLength > 76)
    {
        string ttt = Dummy.substr(60, 10);
        ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
        CO2EnrichmentFactor = atof(ttt.c_str());
        ttt = Dummy.substr(70, 5);
        ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
        DayStartCO2 = atoi(ttt.c_str());
        ttt = Dummy.substr(75);
        ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
        DayEndCO2 = atoi(ttt.c_str());
    }
    else
        CO2EnrichmentFactor = 0;
    //     Line #3: Names of weather files: actual and predicted.
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    if (nLength > 1)
    {
        string ActWthFileName = Dummy.substr(0, 20);
        ActWthFileName.erase(remove(ActWthFileName.begin(), ActWthFileName.end(), ' '), ActWthFileName.end());
        filenames[0] = ActWthFileName;
    }
    if (nLength > 21)
    {
        string PrdWthFileName = Dummy.substr(20, 20);
        PrdWthFileName.erase(remove(PrdWthFileName.begin(), PrdWthFileName.end(), ' '), PrdWthFileName.end());
        filenames[1] = PrdWthFileName;
    }
    //     For advanced users only: If soil mulch is used, read relevant parameters.
    uint32_t MulchIndicator = 0;
    uint32_t DayStartMulch = 0;
    uint32_t DayEndMulch = 0;
    double MulchTranSW = 0;
    double MulchTranLW = 0;
    string m_mulchdata;
    if (nLength > 41)
    {
        m_mulchdata = Dummy.substr(40);
        MulchIndicator = atoi(m_mulchdata.substr(0, 10).c_str());
        if (MulchIndicator > 0)
        {
            MulchTranSW = atof(m_mulchdata.substr(10, 10).c_str());
            MulchTranLW = atof(m_mulchdata.substr(20, 10).c_str());
            DayStartMulch = atoi(m_mulchdata.substr(30, 5).c_str());
            DayEndMulch = atoi(m_mulchdata.substr(35, 5).c_str());
            if (DayEndMulch <= 0)
                DayEndMulch = DateToDoy(DateSimEnd.c_str(), year);
        }
    }
    else
    {
        m_mulchdata = "";
        MulchIndicator = 0;
    }
    DataFile.close();
    //     Calendar dates of emergence, planting, start and stop of simulation, start and stop of
    // output of soil slab and plant maps are converted to DOY dates by calling function DateToDoy.
    uint32_t DayStart = DateToDoy(DateSimStart.c_str(), year);
    uint32_t DayEmerge = DateToDoy(DateEmerge.c_str(), year);
    uint32_t DayFinish = DateToDoy(DateSimEnd.c_str(), year);
    uint32_t DayPlant = DateToDoy(DatePlant.c_str(), year);
    Simulation sim = { ProfileName };
    sim.profile_name_length = strlen(ProfileName);
    sim.year = year;
    sim.day_emerge = DayEmerge;
    sim.day_start = DayStart;
    sim.day_finish = DayFinish;
    sim.day_plant = DayPlant;
    sim.day_start_co2 = DayStartCO2;
    sim.day_end_co2 = DayEndCO2;
    sim.co2_enrichment_factor = CO2EnrichmentFactor;
    sim.day_start_mulch = DayStartMulch;
    sim.day_end_mulch = DayEndMulch;
    sim.mulch_indicator = MulchIndicator;
    sim.mulch_transmissivity_short_wave = MulchTranSW;
    sim.mulch_transmissivity_long_wave = MulchTranLW;
    return sim;
}

//////////////////////////////////////////////////////////
void ReadCalibrationData()
//     This function reads the values of the calibration parameters
//  from input files. It is called from ReadInput(). It calls GetLineData().
//
//     The following global or file-scope variables are set here:
//  SiteName, SitePar, VarName, VarPar
{
    //     Open file of variety file list.
    fs::path strFileName = fs::path("data") / "vars" / "varlist.dat";
    //     If file does not exist, display message and and open a new file
    if (!fs::exists(strFileName))
        throw FileNotExists(strFileName);
    ifstream DataFile(strFileName, ios::in);
    if (DataFile.fail())
        throw FileNotOpened(strFileName);
    //
    string Dummy, VarFile;
    // FIXME: if it go through to the last blank line, it wouldn't break because it doesn't reach the eof
    for (int m_idx = 0; m_idx < 1000; m_idx++)
    {
        if (DataFile.eof() == 1)
            break;
        Dummy = GetLineData(DataFile);
        int nLength = Dummy.length();
        int num;
        string Name, FileName;
        if (nLength >= 4)
        {
            num = atoi(Dummy.substr(0, 4).c_str());
        }
        if (nLength >= 25)
        {
            Name = Dummy.substr(5, 20);
            Name.erase(remove(Name.begin(), Name.end(), ' '), Name.end());
        }
        if (nLength >= 45)
        {
            FileName = Dummy.substr(40, 20);
            FileName.erase(remove(FileName.begin(), FileName.end(), ' '), FileName.end());
        }
        if (num == nVarNum)
        {
            VarFile = FileName;
            VarName = Name;
            break;
        }
    }
    DataFile.close();
    //
    strFileName = fs::path("data") / "vars" / VarFile;
    //  If file does not exist, or can not be opened, display message
    if (!fs::exists(strFileName))
        throw FileNotExists(strFileName);
    ifstream DataFile1(strFileName, ios::in);
    if (DataFile1.fail())
        throw FileNotOpened(strFileName);
    //     Read values of variety related parameters
    GetLineData(DataFile1); // skip 1st line
    for (int i = 1; i <= 60; i++)
    {
        Dummy = GetLineData(DataFile1);
        VarPar[i] = atof(Dummy.substr(0, 20).c_str());
    }
    DataFile1.close();
    //     Open file of site file list.
    strFileName = fs::path("data") / "site" / "sitelist.dat";
    //     If file does not exist, or can not be opened, display message
    if (!fs::exists(strFileName))
        throw FileNotExists(strFileName);
    string SiteFile;
    ifstream DataFile2(strFileName, ios::in);
    if (DataFile2.fail())
        throw FileNotOpened(strFileName);
    //
    for (int m_idx = 0; m_idx < 1000; m_idx++)
    {
        if (DataFile2.eof() == 1)
            break;

        Dummy = GetLineData(DataFile2);
        int nLength = Dummy.length();
        int num;
        string Name, FileName;
        if (nLength >= 4)
        {
            num = atoi(Dummy.substr(0, 4).c_str());
        }
        if (nLength >= 25)
        {
            Name = Dummy.substr(5, 20);
            Name.erase(remove(Name.begin(), Name.end(), ' '), Name.end());
        }
        if (nLength >= 45)
        {
            FileName = Dummy.substr(40, 20);
            FileName.erase(remove(FileName.begin(), FileName.end(), ' '), FileName.end());
        }
        if (num == nSiteNum)
        {
            SiteFile = FileName;
            SiteName = Name;
            break;
        }
    }
    DataFile2.close();
    //
    strFileName = fs::path("data") / "site" / SiteFile;
    //     If file does not exist, or can not be opened, display message
    if (!fs::exists(strFileName))
        throw FileNotExists(strFileName);
    ifstream DataFile3(strFileName, ios::in);
    if (DataFile3.fail())
        throw FileNotOpened(strFileName);
    //     Read values of site related parameters
    GetLineData(DataFile3); // skip 1st line
    for (int i = 1; i <= 20; i++)
    {
        Dummy = GetLineData(DataFile3);
        SitePar[i] = atof(Dummy.substr(0, 20).c_str());
    }
    DataFile3.close();
}

//////////////////////////////////////////////////////////
void InitializeGrid(Simulation &sim)
//     This function initializes the soil grid variables. It is executed once
//  at the beginning of the simulation. It is called from ReadInput().
//
//     The following global or file-scope variables are set here:
//  DensityFactor, dl, nk, nl, PerPlantArea, PlantPopulation,
//  PlantRowColumn, PlantRowLocation, RowSpace, wk.
//     The following global variables are referenced here:
//  PlantsPerM, SkipRowWidth, VarPar, maxk, maxl.
{
    //     PlantRowLocation is the distance from edge of slab, cm, of the plant row.
    PlantRowLocation = 0.5 * sim.row_space;
    if (SkipRowWidth > 1)
    {
        //     If there is a skiprow arrangement, RowSpace and PlantRowLocation are redefined.
        sim.row_space = 0.5 * (sim.row_space + SkipRowWidth); // actual width of the soil slab (cm)
        PlantRowLocation = 0.5 * SkipRowWidth;
    }
    //     Compute PlantPopulation - number of plants per hectar, and PerPlantArea - the average
    //  surface area per plant, in dm2, and the empirical plant density factor (DensityFactor).
    //  This factor will be used to express the effect of plant density on some plant
    //  growth rate functions.  Note that DensityFactor = 1 for 5 plants per sq m (or 50000 per ha).
    PlantPopulation = PlantsPerM / sim.row_space * 1000000;
    PerPlantArea = 1000000 / PlantPopulation;
    DensityFactor = exp(VarPar[1] * (5 - PlantPopulation / 10000));
    //     Define the numbers of rows and columns in the soil slab (nl, nk).
    //  Define the depth, in cm, of consecutive nl layers.
    //     Note that maxl and maxk are defined as constants in file "global.h".
    nl = maxl;
    nk = maxk;
    //      The width of the slab columns is computed by dividing the row
    //  spacing by the number of columns. It is assumed that slab width is
    //  equal to the average row spacing, and column widths are uniform.
    //      Note: wk is an array - to enable the option of non-uniform
    //  column widths in the future.
    //      PlantRowColumn (the column including the plant row) is now computed from
    //  PlantRowLocation (the distance of the plant row from the edge of the slab).
    double sumwk = 0; // sum of column widths
    sim.plant_row_column = 0;
    for (int k = 0; k < nk; k++)
    {
        sumwk = sumwk + wk(k, sim.row_space);
        if (sim.plant_row_column == 0 && sumwk > PlantRowLocation)
        {
            if ((sumwk - PlantRowLocation) > (0.5 * wk(k, sim.row_space)))
                sim.plant_row_column = k - 1;
            else
                sim.plant_row_column = k;
        }
    }
}
