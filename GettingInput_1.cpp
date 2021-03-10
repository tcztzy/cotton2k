//  GettingInput_1.cpp
//
//   functions in this file:
// ReadInput()
// ReadProfileFile()
// InitializeGrid()
//
#include <cmath>
#include <filesystem>
#include <vector>
#include "global.h"
#include "exceptions.h"
#include "GeneralFunctions.h"
#include "Output.h"

using namespace std;
namespace fs = std::filesystem;

static void InitializeGrid(Simulation &);

// GettingInput_2
static void InitSoil(const string &);

static void ReadSoilImpedance(Simulation &);

static void InitializeSoilData(Simulation &, const string &);

static void InitializeSoilTemperature();

static void InitializeRootData(Simulation &);

// GettingInput_3
static int OpenClimateFile(const string &, const string &, const int &, ClimateStruct[400]);

static void ReadAgriculturalInput(Simulation &, const string &);

static double SkipRowWidth,        // the smaller distance between skip rows, cm
    PlantsPerM;             // average number of plants pre meter of row.

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
