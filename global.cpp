// global.cpp : Defines common variables.
// This file also serves as a "dictionary" for these variables.
//
#include "global.h"

//
// Structures:
//
struct scratch Scratch21[400]; // structure used to store daily values of many state variables,
//  see details in file global.h
struct NitrogenFertilizer NFertilizer[150]; // nitrogen fertilizer application information for each day.
// int day = date of application (DOY)
// int mthfrt = method of application ( 0 = broadcast; 1 = sidedress; 2 = foliar; 3 = drip fertigation);
// int ksdr = horizontal placement of side-dressed fertilizer, cm.
// int lsdr = vertical placement of side-dressed fertilizer, cm.
// double amtamm = ammonium N applied, kg N per ha;
// double amtnit = nitrate N applied, kg N per ha;
// double amtura = urea N applied, kg N per ha;
//
// Integer variables:
//
int CultivationDate[5],     // Dates (DOY) of cultivatrion.
    DayFirstDef,            // Date (DOY) of first defoliation.
    DayStartPredIrrig,      // Date (DOY) for starting predicted irrigation.
    DayStopPredIrrig,       // Date (DOY) for stopping predicted irrigation.
    DayWaterTableInput[20], // Dates (DOY) of water table input data.
    DefoliationDate[5],     // Dates (DOY) of defoliant applications.
    DefoliationMethod[5];   // code number of method of application of defoliants:
// 0 = 'banded'; 1 = 'sprinkler'; 2 = 'broaddcast'.

int inrim,       // number of input bulk-density data points for the impedance curve
    IrrigMethod, // method of predicted irrigation.
    isw,         // switch affecting the method of computing soil temperature.
    // 0 = one dimensional (no horizontal flux) - used to predict emergence when emergence date is not known;
    // 1 = one dimensional - used before emergence when emergence date is given;
    // 2 = two dimensional - used after emergence.
    Kday,                  // number of days since emergence.
    LastDayWeatherData,    // last date (DOY) with weather data.
    LastIrrigation,        // date (Doy) of last irrigation (for prediction).
    LastTaprootLayer,      // last soil layer with taproot.
    LateralRootFlag[maxl], // flags indicating presence of lateral roots in soil layers:
    //  0 = no lateral roots are possible.   1 = lateral roots may be initiated.
    //  2 = lateral roots have been initiated.
    LocationColumnDrip, // number of column in which the drip emitter is located
    LocationLayerDrip;  // number of layer in which the drip emitter is located.

int MainStemNodes,       // number of main stem nodes.
    MinDaysBetweenIrrig, // minimum number of days between consecutive irrigations (used for computing predicted irrigation).
    nk,                  // number of vertical columns of soil cells in the slab.
    nl,                  // number of horizontal layers of soil cells in the slab.
    noitr,               // number of iterations per day, for calling some soil water related functions.
    NumAbscisedLeaves,   // number of leaves, per plant, lost by abscission.
    NumFruitBranches[3], // number of fruiting branches at each vegetative branch.
    NumFruitSites,       // total number of fruiting sites per plant.
    NumIrrigations,      // number of irrigations.
    NumNitApps,          // number of applications of nitrogen fertilizer.
    NumNodes[3][30],     // number of nodes on each fruiting branch.
    NumPreFruNodes,      // number of prefruiting nodes, per plant.
    NumSheddingTags,     // number of 'box-car' units used for moving values in arrays defining fruit shedding
    // (AbscissionLag, ShedByCarbonStress, ShedByNitrogenStress and ShedByWaterStress).
    NumVegBranches,    // number of vegetative branches (including the main branch), per plant.
    NumWaterTableData; // number of water table level input data.

int OutIndex[24],          // output flags.
    pixday[10],            // Date (DOY) of PIX application.
    pixmth[10],            // code number for method of application of PIX: 0 = 'BANDED'; 1 = 'SPKLER'; 2 = 'BDCAST'.
    RootColNumLeft[maxl],  // first column with roots in a soil layer.
    RootColNumRight[maxl], // last column with roots in a soil layer.
    SoilHorizonNum[maxl],  // the soil horizon number associated with each soil layer in the slab.
    SoilMapFreq,           // frequency of output of soil slab maps
    WaterTableLayer;       // number of uppermost soil layer below water table.
//
// boolean variables:
//
bool bPollinSwitch; // pollination switch: false = no pollination, true = yes.

// double variables:
//
// A
double
    AbscissionLag[20],            // the time (in physiological days) from tagging fruiting sites for shedding.
    ActualBollGrowth,             // total actual growth of seedcotton in bolls, g plant-1 day-1.
    ActualBurrGrowth,             // total actual growth of burrs in bolls, g plant-1 day-1.
    ActualSquareGrowth,           // total actual growth of squares, g plant-1 day-1.
    ActualSoilEvaporation,        // actual evaporation from soil surface, mm day-1.
    ActualStemGrowth,             // actual growth rate of stems, g plant-1 day-1.
    ActualTranspiration,          // actual transpiration from plants, mm day-1.
    addwtbl,                      // water added to the slab (in mm) due to high water table.
    AgeOfPreFruNode[9],           // age of each prefruiting node, physiological days.
    airdr[9],                     // volumetric water content of soil at "air-dry" for each soil horizon, cm3 cm-3.
    AirTemp[24],                  // hourly air temperatures, C.
    albedo[24],                   // hourly albedo of a reference crop.
    alpha[9],                     // parameter of the Van Genuchten equation.
    AppliedWater,                 // the amount of water to apply (mm), computed for a predicted irrigation.
    AverageLwp,                   // running average of LwpMin + LwpMax for the last 3 days.
    AverageLwpMin,                // running average of LwpMin for the last 3 days.
    AverageSoilPsi,               // average soil matric water potential, bars, computed as the weighted average of the root zone.
    AvrgDailyTemp,                // average daily temperature, C, for 24 hours.
    AvrgNodeTemper[3][30][5];     // running average temperature of each node.
// B
double
    vanGenuchtenBeta[9],  // parameter of the Van Genuchten equation.
    BloomWeightLoss,      // cumulative weight lost due to petals shed after blooming, g per plant.
    BulkDensity[9],       // bulk density of soil in a horizon, g cm-3.
    BurrNConc,            // average nitrogen concentration in burrs.
    BurrNitrogen,         // nitrogen in burrs, g per plant.
    BurrWeightGreenBolls, // total weight of burrs in green bolls, g plant-1.
    BurrWeightOpenBolls;  // total weight of burrs in open bolls, g per plant.
// C
double
    CarbonAllocatedForRootGrowth, // available carbon allocated for root growth, g per plant.
    ClayVolumeFraction[maxl],     // fraction by volume of clay in the soil.
    CloudCoverRatio[24],          // cloud cover ratio (0 to 1).
    CloudTypeCorr[24],            // hourly cloud type correction.
    conmax,                       // the maximum value for non-dimensional hydraulic conductivity
    CottonWeightGreenBolls,       // total weight of seedcotton in green bolls, g plant-1.
    CottonWeightOpenBolls,        // total weight of seedcotton in open bolls, g per plant.
    CultivationDepth[5],          // depth of cultivation, in cm.
    CumEvaporation,               // cumulative evaporation from soil surface, mm.
    CumFertilizerN,               // cumulative amount of inorganic nitrogen applied in fertilizer, mg N per slab.
    CumNetPhotosynth,             // cumulative sum of photosynthate produced to date, g per plant.
    CumNitrogenUptake,            // cumulative total uptake of nitrogen by plants, mg N per slab.
    CumPlantNLoss,                // total cumulative nitrogen lost in sloughed roots, and abscised leaves, squares and bolls, g per plant.
    CumTranspiration,             // cumulative transpiration, mm.
    CumWaterAdded,                // cumulative water added to the slab by rainfall or irrigation, mm.
    CumWaterDrained;              // cumulative water drained out from the slab, mm.
// D
double
    DayTimeTemp,          // average day-time temperature, C.
    dclay,                // aggregation factor for clay in water.
    DeepSoilTemperature,  // boundary soil temperature of deepest layer (K)
    DefoliantAppRate[5],  // rate of defoliant application in pints per acre.
    DelayNewFruBranch[3], // cumulative effect of stresses on delaying the formation of a new fruiting branch.
    DelayNewNode[3][30],  // cumulative effect of stresses on delaying the formation of a new node on a fruiting branch.
    DensityFactor,        // empirical plant density factor.
    DepthLastRootLayer,   // the depth to the end of the last layer with roots (cm).
    DewPointTemp[24],     // hourly dew point temperatures, C.
    dl[maxl],             // the vertical width of a soil layer (cm).
    dsand;                // aggregation factor for sand in water.
// E F G H
double ElCondSatSoil[20],             // electrical conductivity of saturated soil extract (mmho/cm)
    ElCondSatSoilToday,               // electrical conductivity of saturated extract (mmho/cm) on this day.
    es1hour[24],                      // part of hourly Penman evapotranspiration affected by net radiation, in mm per hour.
    es2hour[24],                      // part of hourly Penman evapotranspiration affected by wind and vapor pressure deficit, in mm per hour.
    ExtraCarbon,                      // Extra carbon, not used for plant potential growth requirements, assumed to accumulate in taproot.
    FieldCapacity[maxl],              // volumetric water content of soil at field capacity for each soil layer, cm3 cm-3.
    FoliageTemp[maxk],                // average foliage temperature (oK).
    FreshOrganicMatter[maxl][maxk],   // fresh organic matter in the soil, mg / cm3.
    FreshOrganicNitrogen[maxl][maxk], // N in fresh organic matter in a soil cell, mg cm-3.
    FruitFraction[3][30][5],          // fraction of fruit remaining at each fruiting site (0 to 1).
    FruitGrowthRatio,                 // ratio between actual and potential square and boll growth.
    gh2oc[10],                        // input gravimetric soil water content, g g-1, in the soil mechanical impedance table. values have been read from the soil impedance file.
    ginp,                           // ginning percentage of an individual boll.
    Gintot,                         // weighted average ginning percentage of all open bolls.
    GreenBollsLost,                 // cumulative loss of green bolls, due to abscission, g per plant.
    HeatCapacitySoilSolid[maxl],    // heat capacity of the solid phase of the soil.
    HeatCondDrySoil[maxl],          // the heat conductivity of dry soil.
    HumusNitrogen[maxl][maxk],      // N in stable humic fraction material in a soil cells, mg/cm3.
    HumusOrganicMatter[maxl][maxk]; // humus fraction of soil organic matter, mg/cm3.
// I J K L
double
    impede[10][10],            // input table of soil impedance to root growth
    InitialTotalSoilWater,     // initial total soil water in the profile, mm.
    IrrigationDepth,           // depth of predicted irrigation, cm.
    LeafAreaIndex,             // leaf area index.
    LeafAreaMainStem[3][30],   // mainstem leaf area at each node, dm2.
    LeafAreaNodes[3][30][5],   // leaf area at each fruiting site, dm2.
    LeafAreaPreFru[9],         // area of prefruiting node leaves, dm2.
    LeafNConc,                 // average nitrogen concentration in leaves.
    LeafNitrogen,              // total leaf nitrogen, g per plant.
    LeafWeightAreaRatio,       // temperature dependent factor for converting leaf area to leaf weight during the day, g dm-1.
    LeafWeightMainStem[3][30], // mainstem leaf weight at each node, g.
    LeafWeightNodes[3][30][5], // leaf weight at each fruiting site, g.
    LeafWeightPreFru[9],       // weight of prefruiting node leaves, g.
    LevelsOfWaterTable[20],    // water table level input data (cm below soil surface).
    LightIntercept,            // ratio of light interception by plant canopy.
    LintYield,                 // yield of lint, kgs per hectare.
    LwpMax,                    // maximum (dawn) leaf water potential, MPa.
    LwpMin,                    // minimum (noon) leaf water potential, MPa.
    LwpMinX[3],                // array of values of LwpMin for the last 3 days.
    LwpX[3];                   // array of values of LwpMin + LwpMax for the last 3 days.
// M
double
    MarginalWaterContent[maxl], // marginal soil water content (as a function of soil texture) for computing soil heat conductivity.
    MaxIrrigation,              // maximum amount of applied water in a predicted irrigation
    MaxWaterCapacity[maxl],     // volumetric water content of a soil layer at maximum capacity, before drainage, cm3 cm-3.
    MineralizedOrganicN,        // cumulative amount of mineralized organic N, mgs per slab.
    MulchTemp[maxk];            // polyethylene mulch temperature at a column(oK).
// N
double
    NetPhotosynthesis,     // net photosynthetic rate, g per plant per day.
    NightTimeTemp,         // average night-time temperature, C.
    NitrogenStress,        // the average nitrogen stress coefficient for vegetative and reproductive organs
    NO3FlowFraction[maxl], // fraction of nitrate that can move to the next layer.
    NStressFruiting,       // nitrogen stress limiting fruit development.
    NStressRoots,          // nitrogen stress limiting root development.
    NStressVeg,            // nitrogen stress limiting vegetative development.
    NumGreenBolls,         // average number of retained green bolls, per plant.
    NumOpenBolls,          // number of open bolls, per plant.
    NumSquares;            // number of squares per plant.
// P
double
    PercentDefoliation,                 // percentage of leaves abscised as a result of defoliant application.
    PerPlantArea,                       // average soil surface area per plant, dm2.
    PetioleNConc,                       // average nitrogen concentration in petioles.
    PetioleNitrogen,                    // total petiole nitrogen, g per plant.
    PetioleNO3NConc,                    // average nitrate nitrogen concentration in petioles.
    PetioleWeightMainStem[3][30],       // weight of mainstem leaf petiole at each node, g.
    PetioleWeightPreFru[9],             // weight of prefruiting node petioles, g.
    pixcon,                             // PIX concentration.
    pixda,                              // effect of PIX application on the growth rate of leaves.
    pixdn,                              // effect of PIX application on delaying the appearance of new fruiting nodes.
    pixdz,                              // effect of PIX application on the growth rate of stem weight and of height.
    PixInPlants,                        // total amount of PIX in the plant, g per plant.
    pixppa[10],                         // rate of PIX application, pints per acre.
    PlantPopulation,                    // plant population, plants per hectar.
    PlantRowLocation,                   // distance of plant row from slab edge, cm.
    PlantWeight,                        // total plant weight, g.
    PlantWeightAtStart,                 // total plant dry matter at germination, g per plant.
    PoreSpace[maxl],                    // pore space of soil, volume fraction.
    PotGroAllBolls,                     // sum of potential growth rates of seedcotton in all bolls, g plant-1 day-1.
    PotGroAllBurrs,                     // sum of potential growth rates of burrs in all bolls, g plant-1 day-1.
    PotGroAllLeaves,                    // sum of potential growth rates of all leaves, g plant-1 day-1.
    PotGroAllPetioles,                  // sum of potential growth rates of all petioles, g plant-1 day-1.
    PotGroAllRoots,                     // potential growth rate of roots, g plant-1 day-1
    PotGroAllSquares,                   // sum of potential growth rates of all squares, g plant-1 day-1.
    PotGroBurrs[3][30][5],              // potential growth rate of burrs in an individual boll, g day-1.
    PotGroLeafAreaMainStem[3][30],      // potential growth in area of an individual main stem node leaf, dm2 day-1.
    PotGroLeafAreaNodes[3][30][5],      // potential growth in area of an individual fruiting node leaf, dm2 day-1.
    PotGroLeafAreaPreFru[9],            // potentially added area of a prefruiting node leaf, dm2 day-1.
    PotGroLeafWeightMainStem[3][30],    // potential growth in weight of an individual main stem node leaf, g day-1.
    PotGroLeafWeightNodes[3][30][5],    // potential growth in weight of an individual fruiting node leaf, g day-1.
    PotGroLeafWeightPreFru[9],          // potentially added weight of a prefruiting node leaf, g day-1.
    PotGroPetioleWeightMainStem[3][30], // potential growth in weight of an individual main stem node petiole, g day-1.
    PotGroPetioleWeightNodes[3][30][5], // potential growth in weight of an individual fruiting node petiole, g day-1.
    PotGroPetioleWeightPreFru[9],       // potentially added weight of a prefruiting node petiole, g day-1.
    PotGroRoots[maxl][maxk],            // potential root growth in a soil cell (g per day).
    PotGroSquares[3][30][5],            // potential growth rate of an individual square, g day-1.
    PotGroStem;                         // potential growth rate of stems, g plant-1 day-1.
// R
double
    Radiation[24],                 // hourly global radiation, W / m2.
    RatioImplicit,                 // the ratio for the implicit numerical solution of the water transport equation (used in FLUXI and in SFLUX.
    ReferenceETP[24],              // reference evapotranspiration, mm per hour.
    ReferenceTransp,               // daily sum of hourly reference evapotranspiration, mm per day.
    RelativeHumidity[24],          // hourly values of relative humidity (%).
    ReserveC,                      // reserve carbohydrates in leaves, g per plant.
    rlat1[maxl],                   // lateral root length (cm) to the left of the tap root
    rlat2[maxl],                   // lateral root length (cm) to the right of the tap root
    Rn,                            // daily total net radiation, W m-2.
    RootGroFactor[maxl][maxk],     // root growth correction factor in a soil cell (0 to 1).
    RootImpede[maxl][maxk],        // root mechanical impedance for a soil cell, kg cm-2.
    RootNConc,                     // average nitrogen concentration in roots.
    RootNitrogen,                  // total root nitroge, g per plant.
    RootWeightLoss;                // total cumulative weight of sloughed roots, g per plant.
// S
double
    SandVolumeFraction[maxl],      // fraction by volume of sand plus silt in the soil.
    SaturatedHydCond[9],           // saturated hydraulic conductivity, cm per day.
    SeedNConc,                     // average nitrogen concentration in seeds.
    SeedNitrogen,                  // total seed nitrogen, g per plant.
    ShedByCarbonStress[20],        // the effect of carbohydrate stress on shedding
    ShedByNitrogenStress[20],      // the effect of nitrogen stress on shedding.
    ShedByWaterStress[20],         // the effect of moisture stress on shedding.
    SitePar[21],                   // array of site specific constant parameters.
    SoilNitrogenAtStart,           // total soil nitrogen at the start of simulation, mg per slab.
    SoilNitrogenLoss,              // cumulative loss of nitrogen by drainage out of the lowest soil layer, mg per slab.
    SoilPsi[maxl][maxk],           // matric water potential of a soil cell, bars.
    SoilTemp[maxl][maxk],          // hourly soil temperature oK.
    SoilTempDailyAvrg[maxl][maxk], // daily average soil temperature, oK.
    SquareNConc,                   // average concentration of nitrogen in the squares.
    SquareNitrogen,                // total nitrogen in the squares, g per plant
    SquareWeight[3][30][5],        // weight of each square,g
    StemNConc,                     // ratio of stem nitrogen to dry matter.
    StemNitrogen,                  // total stem nitrogen, g per plant
    StemWeight[365],               // stem weight added at each day, g per plant.
    SumNO3N90,                     // sum of soil nitrate n, 0-90 cm depth, in kg/ha.
    SupplyNH4N,                    // uptake of ammonia N by the plant from the soil, mg N per slab per day.
    SupplyNO3N;                    // uptake of nitrate by the plant from the soil, mg N per slab per day.
// T
double
    TapRootLength,            // the length of the taproot, in cm.
    thad[maxl],               // residual volumetric water content of soil layers (at air-dry condition), cm3 cm-3.
    thetar[maxl],             // volumetric water content of soil layers at permanent wilting point (-15 bars), cm3 cm-3.
    thetas[9],                // volumetric saturated water content of soil horizon, cm3 cm-3.
    thts[maxl],               // saturated volumetric water content of each soil layer, cm3 cm-3.
    TotalActualLeafGrowth,    // actual growth rate of all the leaves, g plant-1 day-1.
    TotalActualPetioleGrowth, // actual growth rate of all the petioles, g plant-1 day-1.
    TotalLeafArea,            // total leaf area, dm2 per plant.
    TotalLeafWeight,          // total leaf weight, g per plant.
    TotalPetioleWeight,       // total petiole weight, g per plant.
    TotalRequiredN,           // total nitrogen required for plant growth, g per plant.
    TotalRootWeight,          // total dry weight of the roots, g per plant.
    TotalSoilNh4N,            // total ammonium in profile, mg N per slab.
    TotalSoilNitrogen,        // total soil nitrogen, mg N per slab.
    TotalSoilNo3N,            // total nitrate in profile, mg N per slab.
    TotalSoilUreaN,           // total urea in profile, mg N per slab.
    TotalSoilWater,           // total water in the soil profile, mm.
    TotalSquareWeight,        // total weight of squares, g per plant
    TotalStemWeight,          // total stem weight, g per plant.
    tstbd[10][10];            // input bulk density in the impedance table, g cm-3.
// U V W
double
    VarPar[61],                  // array of cultivar specific constant parameters.
    VolNh4NContent[maxl][maxk],  // volumetric ammonium nitrogen content of a soil cell, mg N cm-3.
    VolNo3NContent[maxl][maxk],  // volumetric nitrate nitrogen content of a soil cell, mg N cm-3.
    VolUreaNContent[maxl][maxk], // volumetric urea nitrogen content of a soil cell, mg N cm-3.
    VolWaterContent[maxl][maxk], // volumetric water content of a soil cell, cm3 cm-3.
    WaterStressStem,             // water stress index for stem growth (0 to 1).
    WindSpeed[24],               // Hourly wind velocity, m per second.
    wk[maxk];                    // horizontal width of a soil column (cm).
