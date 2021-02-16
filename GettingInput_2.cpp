//  File GettingInput_2.cpp
//     Consists of the following functions:
// ReadSoilImpedance()
// InitSoil()
// InitializeSoilData()
// ReadSoilHydraulicData()
// InitializeRootData()
// InitializeSoilTemperature()
// form()
//
#include "global.h"
#include "GeneralFunctions.h"
#include "Input.h"

int ReadSoilHydraulicData(const string &);

//
// Definitions of variables with file scope:
//
double condfc[9],        // hydraulic conductivity at field capacity of horizon layers, cm per day.
h2oint[14],          // initial soil water content, percent of field capacity,
// defined by input for consecutive 15 cm soil layers.
ldepth[9],           // depth from soil surface to the end of horizon layers, cm.
oma[14],             // organic matter at the beginning of the season, percent of soil weight,
// defined by input for consecutive 15 cm soil layers.
pclay[9],            // percentage of clay in soil horizon of horizon layers.
psand[9],            // percentage of sand in soil horizon of horizon layers.
psidra,              // soil matric water potential, bars, for which immediate drainage
// will be simulated (suggested value -0.25 to -0.1).
psisfc,              // soil matric water potential at field capacity,
// bars (suggested value -0.33 to -0.1).
rnnh4[14],           // residual nitrogen as ammonium in soil at beginning of season, kg per ha.
// defined by input for consecutive 15 cm soil layers.
rnno3[14];           // residual nitrogen as nitrate in soil at beginning of season, kg per ha.
// defined by input for consecutive 15 cm soil layers.
////////////////////////////////////////////////////////////////////////////////////
void ReadSoilImpedance(Simulation & sim)
//     This function opens the soil root impedance data file and reads it. 
//  It is called from ReadInput(), and executed once at the beginning of the simulation.
//  The variables read here are later used to compute soil impedance to root growth.
//
//     Global or file scope variables set:
//        gh2oc, impede, inrim, ncurve, tstbd.
//
{
    fs::path strFileName = fs::path("data") / "soilimpd.dat";
//  If file does not exist, or can not be opened, display message.
    if (!fs::exists(strFileName)) throw FileNotExists(strFileName);
    ifstream DataFile(strFileName, ios::in);
    if (DataFile.fail()) throw FileNotOpened(strFileName);
//     Read data from file
    string Dummy = GetLineData(DataFile); // 1st line
    string SoilName = Dummy;
    Dummy = GetLineData(DataFile); // 2nd line
    int ncurve = atoi(Dummy.c_str());
    if (ncurve > 10)
        ncurve = 10;
    sim.num_curve = ncurve;
    for (int i = 0; i < ncurve; i++) {
        DataFile >> inrim >> gh2oc[i];
        if (inrim > 10)
            inrim = 10;
        for (int j = 0; j < inrim; j++)
            DataFile >> tstbd[j][i] >> impede[j][i];
    }
    DataFile.close();
}

//////////////////////////////////////////////////////////////////
void InitSoil(const string &SoilInitFileName)
//     This function opens the initial soil data file and reads it. It is executed
//  once at the beginning of the simulation. It is called by ReadInput().
//
//     Global or file scope variables set:
//        rnnh4, rnno3, oma, h2oint.
{
//     The the file containing the initial soil data is opened for input.
    fs::path strFileName = fs::path("soils") / SoilInitFileName;
    ifstream fstr(strFileName, ios::in);
    if (fstr.fail())
        throw FileNotOpened(strFileName);
//     The description of this file is read from the first line of the file.
    char dscrip[100];
    fstr.getline(dscrip, 99);
//     The next 12 lines contain data for 12 consecutive 15 cm soil
//  layers (the last layer is 35 cm).
//      (1) residual ammonium N, kgs per ha.
//      (2) residual nitrate N, kgs per ha.
//      (3) organic matter, percent of soil dry weight.
//      (4) soil water content, percent of field capacity.
    for (int i = 0; i < 12; i++)
        fstr >> rnnh4[i] >> rnno3[i] >> oma[i] >> h2oint[i];
    fstr.close();
//     Data for soil layers below 180 cm are assumed to be the same
//  as those for the last defined layer (layer 11 is 165-200 cm).
    for (int i = 12; i < 14; i++) {
        h2oint[i] = h2oint[11];
        rnnh4[i] = rnnh4[11];
        rnno3[i] = rnno3[11];
        oma[i] = oma[11];
    }
}

//////////////////////////////////////////////////////////
void InitializeSoilData(Simulation &sim, const string &SoilHydFileName)
//     This function computes and sets the initial soil data. It is
//  executed once at the beginning of the simulation, after the soil
//  hydraulic data file has been read. It is called by ReadInput().
//     It calls ReadSoilHydraulicData().
//
//     Global and file scope variables referenced:
//        airdr, alpha, vanGenuchtenBeta, BulkDensity, condfc, Date, dl, h2oint, ldepth, 
//        nl, oma, OutIndex, psisfc, wk.
//     Global and file scope variables set:
//        FieldCapacity, FreshOrganicMatter, HumusOrganicMatter, InitialTotalSoilWater, 
//        MaxWaterCapacity, NO3FlowFraction, PoreSpace, rnnh4, rnno3, SaturatedHydCond, 
//        SoilHorizonNum, thad, thetar, thetas, thts, TotalSoilNh4N, TotalSoilNo3N, TotalSoilUreaN, 
//        VolNh4NContent, VolNo3NContent, VolUreaNContent, VolWaterContent
{
    int lyrsol;  // the number of soil horizons in the slab (down to 2 m).
    lyrsol = ReadSoilHydraulicData(SoilHydFileName);
//
    int j = 0; // horizon number
    double sumdl = 0; // depth to the bottom this layer (cm);
    double rm = 2.65; // density of the solid fraction of the soil (g / cm3)
    double bdl[40];  // array of bulk density of soil layers
    for (int l = 0; l < nl; l++) {
//     Using the depth of each horizon layer (ldepth), the horizon
//  number (SoilHorizonNum) is computed for each soil layer.
        sumdl += dl[l];
        while ((sumdl > ldepth[j]) && (j < lyrsol))
            j++;
        SoilHorizonNum[l] = j;
//     bdl, thad, thts are defined for each soil layer, using the
//  respective input variables BulkDensity, airdr, thetas.
//      FieldCapacity, MaxWaterCapacity and thetar are computed for each layer, as water
//  content (cm3 cm-3) of each layer corresponding to matric potentials
//  of psisfc (for field capacity), psidra (for free drainage) and -15
//  bars (for permanent wilting point), respectively, using function
//  qpsi.
//      pore space volume (PoreSpace) is also computed for each layer.
//      make sure that saturated water content is not more than pore space.
        bdl[l] = BulkDensity[j];
        PoreSpace[l] = 1 - BulkDensity[j] / rm;
        if (thetas[j] > PoreSpace[l])
            thetas[j] = PoreSpace[l];
        thad[l] = airdr[j];
        thts[l] = thetas[j];
        FieldCapacity[l] = qpsi(psisfc, thad[l], thts[l], alpha[j], vanGenuchtenBeta[j]);
        MaxWaterCapacity[l] = qpsi(psidra, thad[l], thts[l], alpha[j], vanGenuchtenBeta[j]);
        thetar[l] = qpsi(-15., thad[l], thts[l], alpha[j], vanGenuchtenBeta[j]);
//     When the saturated hydraulic conductivity (SaturatedHydCond) is not
//  given, it is computed from the hydraulic conductivity at field
//  capacity (condfc), using the wcond function.
        if (SaturatedHydCond[j] <= 0)
            SaturatedHydCond[j] = condfc[j] / wcond(FieldCapacity[l], thad[l], thts[l], vanGenuchtenBeta[j], 1, 1);
    }
//     Loop for all soil layers. Compute depth from soil surface to
//  the end of each layer (sumdl).
    sumdl = 0;
    for (int l = 0; l < nl; l++) {
        sumdl += dl[l];
//     At start of simulation compute estimated movable fraction of
//  nitrates in each soil layer, following the work of:
//     Bowen, W.T., Jones, J.W., Carsky, R.J., and Quintana, J.O. 1993.
//  Evaluation of the nitrogen submodel of CERES-maize following legume
//  green manure incorporation. Agron. J. 85:153-159.
//     The fraction of total nitrate in a layer that is in solution and
//  can move from one layer to the next with the downward flow of
//  water, FLOWNO3[l], is a function of the adsorption coefficient,
//  soil bulk density, and the volumetric soil water content at the
//  drained upper limit.
//     Adsorption coefficients are assumed to be 0.0 up to 30 cm depth,
//  and deeper than 30 cm - 0.2, 0.4, 0.8, 1.0, 1.2, and 1.6 for each
//  successive 15 cm layer.
        double coeff; // Adsorption coefficient
        if (sumdl <= 30)
            coeff = 0;
        else if (sumdl <= 45)
            coeff = 0.2;
        else if (sumdl <= 60)
            coeff = 0.4;
        else if (sumdl <= 75)
            coeff = 0.6;
        else if (sumdl <= 90)
            coeff = 0.8;
        else if (sumdl <= 105)
            coeff = 1.0;
        else if (sumdl <= 120)
            coeff = 1.2;
        else
            coeff = 1.6;
        NO3FlowFraction[l] = 1 / (1 + coeff * bdl[l] / MaxWaterCapacity[l]);
//     Determine the corresponding 15 cm layer of the input file.
//     Compute the initial volumetric water content (VolWaterContent) of each
//  layer, and check that it will not be less than the air-dry value or
//  more than pore space volume.
        j = (int) ((sumdl - 1) / 15);
        if (j > 13)
            j = 13;
        int n = SoilHorizonNum[l];
        VolWaterContent[l][0] = FieldCapacity[l] * h2oint[j] / 100;
        if (VolWaterContent[l][0] < airdr[n])
            VolWaterContent[l][0] = airdr[n];
        if (VolWaterContent[l][0] > PoreSpace[l])
            VolWaterContent[l][0] = PoreSpace[l];
//     Initial values of ammonium N (rnnh4, VolNh4NContent) and nitrate N
//  (rnno3, VolNo3NContent) are converted from kgs per ha to mg / cm3 for each
//  soil layer, after checking for minimal amounts.
        if (rnno3[j] < 2.0)
            rnno3[j] = 2.0;
        if (rnnh4[j] < 0.2)
            rnnh4[j] = 0.2;
        VolNo3NContent[l][0] = rnno3[j] / 15 * .01;
        VolNh4NContent[l][0] = rnnh4[j] / 15 * .01;
        double om; // organic matter in mg / cm3 units.
        om = (oma[j] / 100) * bdl[l] * 1000;
//     potom is the proportion of readily mineralizable om. it is a
//  function of soil depth (sumdl, in cm), modified from GOSSYM (where
//  it probably includes the 0.4 factor for organic C in om).
        double potom = max(0., 0.15125 - 0.02878 * log(sumdl));
//     FreshOrganicMatter is the readily mineralizable organic matter (= "fresh organic
//  matter" in CERES models). HumusOrganicMatter is the remaining organic matter, which
//  is mineralized very slowly.
        FreshOrganicMatter[l][0] = om * potom;
        HumusOrganicMatter[l][0] = om * (1 - potom);
    }
//     Since the initial value has been set for the first column only
//  in each layer, these values are now assigned to all the other columns.
    for (int l = 0; l < nl; l++)
        for (int k = 1; k < nk; k++) {
            VolWaterContent[l][k] = VolWaterContent[l][0];
            VolNo3NContent[l][k] = VolNo3NContent[l][0];
            VolNh4NContent[l][k] = VolNh4NContent[l][0];
            FreshOrganicMatter[l][k] = FreshOrganicMatter[l][0];
            HumusOrganicMatter[l][k] = HumusOrganicMatter[l][0];
        }
//     Total amounts of water (InitialTotalSoilWater), nitrate N (TotalSoilNo3N), ammonium
//  N (TotalSoilNh4N), and urea N (TotalSoilUreaN) are computed for the whole slab.
    InitialTotalSoilWater = 0;
    TotalSoilNo3N = 0;
    TotalSoilNh4N = 0;
    TotalSoilUreaN = 0;
//
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++) {
            InitialTotalSoilWater += VolWaterContent[l][k] * dl[l] * wk[k];
            TotalSoilNo3N += VolNo3NContent[l][k] * dl[l] * wk[k];
            TotalSoilNh4N += VolNh4NContent[l][k] * dl[l] * wk[k];
            VolUreaNContent[l][k] = 0;
        }
//     InitialTotalSoilWater is converted from cm3 per slab to mm.
    InitialTotalSoilWater = 10 * InitialTotalSoilWater / sim.row_space;
}

/////////////////////////////////////////////////////////////////////
int ReadSoilHydraulicData(const string &SoilHydFileName)
//     This functions reads the soil hydrology file. It is called by InitializeSoilData().
//     Global or file scope Variables set here:
//  airdr, alpha, vanGenuchtenBeta, BulkDensity, condfc, conmax, ldepth, pclay, psand, 
//  psidra, psisfc, RatioImplicit, SaturatedHydCond, thetas
//     Return value = number of soilhorizons (lyrsol)
{
//     The the file containing the hydraulic soil data is opened for input.
    fs::path m_FileName = fs::path("soils") / SoilHydFileName;
    ifstream DataFile(m_FileName, ios::in);
    if (DataFile.fail())
        throw FileNotOpened(m_FileName);
//     Zeroise arrays of data values.
    for (int i = 0; i < 9; i++) {
        ldepth[i] = 0;
        airdr[i] = 0;
        thetas[i] = 0;
        alpha[i] = 0;
        vanGenuchtenBeta[i] = 0;
        SaturatedHydCond[i] = 0;
        condfc[i] = 0;
        BulkDensity[i] = 0;
        pclay[i] = 0;
        psand[i] = 0;
    }
//     The following code reads data from the file to the variables: 
//     -Reading line 1--------------------------------------------------------------//
    string Dummy = GetLineData(DataFile);
    string m_fileDesc = ""; // Description of the Hydraulic file
    if (Dummy.length() > 1)
        m_fileDesc = Dummy.erase(Dummy.find_last_not_of(" \r\n\t\f\v") + 1);
//     -Reading line 2--------------------------------------------------------------//
    Dummy = GetLineData(DataFile);
    int lyrsol;   // the number of soil horizons in the slab (down to 2 m).
    lyrsol = atoi(Dummy.substr(0, 10).c_str());
    RatioImplicit = atof(Dummy.substr(10, 10).c_str());
    conmax = atof(Dummy.substr(20, 10).c_str());
    psisfc = atof(Dummy.substr(30, 10).c_str());
    psidra = atof(Dummy.substr(40, 10).c_str());
    if (lyrsol > 9)
        lyrsol = 9; // not more than 9 layers
    if (psisfc > 0)
        psisfc = -psisfc; // make sure it is negative
    if (psidra > 0)
        psidra = -psidra; // make sure it is negative
//     -Reading next lines  -------------------------------------------
    for (int il = 0; il < lyrsol; il++) {
//     First line for each layer
        Dummy = GetLineData(DataFile);
        ldepth[il] = atof(Dummy.substr(0, 10).c_str());
        airdr[il] = atof(Dummy.substr(10, 10).c_str());
        thetas[il] = atof(Dummy.substr(20, 10).c_str());
        alpha[il] = atof(Dummy.substr(30, 10).c_str());
        vanGenuchtenBeta[il] = atof(Dummy.substr(40, 10).c_str());
        SaturatedHydCond[il] = atof(Dummy.substr(50, 10).c_str());
        condfc[il] = atof(Dummy.substr(60, 10).c_str());
//     Second line for each layer
        Dummy = GetLineData(DataFile);
        BulkDensity[il] = atof(Dummy.substr(0, 10).c_str());
        pclay[il] = atof(Dummy.substr(10, 10).c_str());
        psand[il] = atof(Dummy.substr(20, 10).c_str());
    }
    DataFile.close();
    return lyrsol;
}

void init_root_data(SoilCell soil_cells[40][20], uint32_t plant_row_column, double mul) {
    for (int l = 0; l < 40; l++) {
        for (int k = 0; k < 20; k++) {
            soil_cells[l][k] = {{0, 1, 0, 0, 0, {0, 0, 0}}};
        }
    }
    // FIXME: I consider the value is incorrect
    soil_cells[0][plant_row_column - 1].root.weight[0] = 0.0020 * mul;
    soil_cells[0][plant_row_column].root.weight[0] = 0.0070 * mul;
    soil_cells[0][plant_row_column + 1].root.weight[0] = 0.0070 * mul;
    soil_cells[0][plant_row_column + 2].root.weight[0] = 0.0020 * mul;
    soil_cells[1][plant_row_column - 1].root.weight[0] = 0.0040 * mul;
    soil_cells[1][plant_row_column].root.weight[0] = 0.0140 * mul;
    soil_cells[1][plant_row_column + 1].root.weight[0] = 0.0140 * mul;
    soil_cells[1][plant_row_column + 2].root.weight[0] = 0.0040 * mul;
    soil_cells[2][plant_row_column - 1].root.weight[0] = 0.0060 * mul;
    soil_cells[2][plant_row_column].root.weight[0] = 0.0210 * mul;
    soil_cells[2][plant_row_column + 1].root.weight[0] = 0.0210 * mul;
    soil_cells[2][plant_row_column + 2].root.weight[0] = 0.0060 * mul;
    soil_cells[3][plant_row_column].root.weight[0] = 0.0200 * mul;
    soil_cells[3][plant_row_column + 1].root.weight[0] = 0.0200 * mul;
    soil_cells[4][plant_row_column].root.weight[0] = 0.0150 * mul;
    soil_cells[4][plant_row_column + 1].root.weight[0] = 0.0150 * mul;
    soil_cells[5][plant_row_column].root.weight[0] = 0.0100 * mul;
    soil_cells[5][plant_row_column + 1].root.weight[0] = 0.0100 * mul;
    soil_cells[6][plant_row_column].root.weight[0] = 0.0050 * mul;
    soil_cells[6][plant_row_column + 1].root.weight[0] = 0.0050 * mul;
    for (int l = 0; l < 3; l++) {
        for (int k = plant_row_column - 1; k < plant_row_column + 3; k++) {
            soil_cells[l][k].root.age = 0.01;
        }
    }
    for (int l = 3; l < 7; l++) {
        soil_cells[l][plant_row_column].root.age = 0.01;
        soil_cells[l][plant_row_column + 1].root.age = 0.01;
    }
}

//////////////////////////////////////////////////////////
void InitializeRootData(Simulation & sim)
//     This function initializes the root submodel parameters and variables. It is called
//  by ReadInput(). it is executed once at the beginning of the simulation. 
//
//     Global or file scope variables referenced:
//        dl, PlantRowColumn, nk, nl, PerPlantArea.
//     Global or file scope variables set:
//        ActualRootGrowth[maxl][maxk], DepthLastRootLayer,
//	      LastTaprootLayer, LateralRootFlag[maxl], NumLayersWithRoots, NumRootAgeGroups, 
//        PotGroRoots[maxl][maxk], RootAge[maxl][maxk], RootColNumLeft[maxl], 
//        RootColNumRight[maxl], RootGroFactor[maxl][maxk], RootWeight[maxl][maxk][3], 
//        TapRootLength, TotalRootWeight.
//
{
//     The parameters of the root model are defined for each root class:
//       grind(i), cuind(i), thtrn(i), trn(i), thdth(i), dth(i).
    double rlint = 10; // Vertical interval, in cm, along the taproot, for
    // initiating lateral roots.
    int ll = 1; // Counter for layers with lateral roots.
    double sumdl = 0; // Distance from soil surface to the middle of a soil layer.
    for (int l = 0; l < nl; l++) {
//     Using the value of rlint (interval between lateral roots), the
//  layers from which lateral roots may be initiated are now computed.
//  LateralRootFlag[l] is assigned a value of 1 for these layers.
        LateralRootFlag[l] = 0;
        if (l > 0)
            sumdl += 0.5 * dl[l - 1];
        sumdl += 0.5 * dl[l];
        if (sumdl >= ll * rlint) {
            LateralRootFlag[l] = 1;
            ll++;
        }
    }
//     All the state variables of the root system are initialized to zero.
    for (int l = 0; l < nl; l++) {
        if (l < 3) {
            RootColNumLeft[l] = sim.plant_row_column - 1;
            RootColNumRight[l] = sim.plant_row_column + 2;
        } else if (l < 7) {
            RootColNumLeft[l] = sim.plant_row_column;
            RootColNumRight[l] = sim.plant_row_column + 1;
        } else {
            RootColNumLeft[l] = 0;
            RootColNumRight[l] = 0;
        }
    }
    init_root_data(sim.states[0].soil_cells, sim.plant_row_column, 0.01 * sim.row_space / PerPlantArea);
//     Start loop for all soil layers containing roots.
    DepthLastRootLayer = 0;
    TotalRootWeight = 0;
    for (int l = 0; l < 7; l++) {
        DepthLastRootLayer += dl[l]; //compute total depth to the last layer with roots (DepthLastRootLayer).
//     For each soil soil cell with roots, compute total root weight 
//  per plant (TotalRootWeight), and convert RootWeight from g per plant to g per cell.
        for (int k = 0; k < nk; k++) {
            for (int i = 0; i < 3; i++) {
                TotalRootWeight += sim.states[0].soil_cells[l][k].root.weight[i] * 100 / sim.row_space * PerPlantArea;
            }
        }
    }
//     Initial value of taproot length, TapRootLength, is computed to the
// middle of the last layer with roots. The last soil layer with
// taproot, LastTaprootLayer, is defined.
    int NumLayersWithRoots = 7;
    TapRootLength = (DepthLastRootLayer - 0.5 * dl[NumLayersWithRoots - 1]);
    LastTaprootLayer = 6;
}

//////////////////////////////////////////////////////////
void InitializeSoilTemperature()
//     This function initializes the variables needed for the simulation of 
// soil temperature, and variables used by functions ThermalCondSoil() and SoilHeatFlux().
//     It is executed once at the beginning of the simulation. It is called by 
// ReadInput() and it calls form().
//
//     Global and file scope variables set:
//	      MarginalWaterContent[maxl], dclay, dsand, HeatCapacitySoilSolid[maxl], 
//        HeatCondDrySoil[maxl], ClayVolumeFraction[maxl], SandVolumeFraction[maxl];
//
//     Global and file scope variables referenced:
//        dl[maxl], nl, oma[14], pclay[9], psand[9], PoreSpace[maxl], SoilHorizonNum[maxl].
//
{
    double bsand = 20;   // heat conductivity of sand and silt (mcal cm-1 s-1 C-1).
    double bclay = 7;    // heat conductivity of clay (mcal cm-1 s-1 C-1).
    double cka = 0.0615; // heat conductivity of air (mcal cm-1 s-1 C-1).
    double ckw = 1.45;   // heat conductivity of water (mcal cm-1 s-1 C-1).
    double cmin = 0.46;  // heat capacity of the mineral fraction of the soil.
    double corg = 0.6;   // heat capacity of the organic fraction of the soil.
    double ga = 0.144;   // shape factor for air in pore spaces.
    double rm = 2.65;    // specific weight of mineral fraction of soil.
    double ro = 1.3;     // specific weight of organic fraction of soil.
//     Compute aggregation factors:
    dsand = form(bsand, ckw, ga);    // aggregation factor for sand in water
    dclay = form(bclay, ckw, ga);    // aggregation factor for clay in water
    double dsandair = form(bsand, cka, ga); // aggregation factor for sand in air
    double dclayair = form(bclay, cka, ga); // aggregation factor for clay in air
//     Loop over all soil layers, and define indices for some soil arrays.
    double sumdl = 0; // sum of depth of consecutive soil layers.
    for (int l = 0; l < nl; l++)  // loop by soil layers
    {
        sumdl += dl[l];
        int j = (int) ((sumdl + 14) / 15) - 1;   //  layer definition for oma
        if (j > 13)
            j = 13;
//     Using the values of the clay and organic matter percentages in the soil, compute
//   mineral and organic fractions of the soil, by weight and by volume.
        double mmo = oma[j] / 100; // organic matter fraction of dry soil (by weight).
        double mm = 1 - mmo;       // mineral fraction of dry soil (by weight).
//     MarginalWaterContent is set as a function of the sand fraction of the soil.
        int i1 = SoilHorizonNum[l];   //  layer definition as in soil hydrology input file.
        MarginalWaterContent[l] = 0.1 - 0.07 * psand[i1] / 100;
//     The volume fractions of clay (ClayVolumeFraction) and of sand plus silt (SandVolumeFraction), are
//  calculated.
        double ra = (mmo / ro) / (mm / rm);        // volume ratio of organic to mineral soil fractions.
        double xo = (1 - PoreSpace[l]) * ra / (1 + ra); // organic fraction of soil (by volume).
        double xm = (1 - PoreSpace[l]) - xo;             // mineral fraction of soil (by volume).
        ClayVolumeFraction[l] = pclay[i1] * xm / mm / 100;
        SandVolumeFraction[l] = 1 - PoreSpace[l] - ClayVolumeFraction[l];
//     Heat capacity of the solid soil fractions (mineral + organic, by volume )
        HeatCapacitySoilSolid[l] = xm * cmin + xo * corg;
//     The heat conductivity of dry soil (HeatCondDrySoil) is computed using the
//  procedure suggested by De Vries.
        HeatCondDrySoil[l] = 1.25 * (PoreSpace[l] * cka + dsandair * bsand * SandVolumeFraction[l] +
                                     dclayair * bclay * ClayVolumeFraction[l])
                             / (PoreSpace[l] + dsandair * SandVolumeFraction[l] + dclayair * ClayVolumeFraction[l]);
    }
}
