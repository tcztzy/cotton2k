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
