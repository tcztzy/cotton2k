#include "global.h"

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
unsigned int ncurve;

int CultivationDate[5],     // Dates (DOY) of cultivatrion.
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
    LastDayWeatherData,    // last date (DOY) with weather data.
    LastIrrigation,        // date (Doy) of last irrigation (for prediction).
    LastTaprootLayer,      // last soil layer with taproot.
    LateralRootFlag[maxl], // flags indicating presence of lateral roots in soil layers:
    //  0 = no lateral roots are possible.   1 = lateral roots may be initiated.
    //  2 = lateral roots have been initiated.
    LocationColumnDrip,  // number of column in which the drip emitter is located
    LocationLayerDrip,   // number of layer in which the drip emitter is located.
    MainStemNodes,       // number of main stem nodes.
    MinDaysBetweenIrrig, // minimum number of days between consecutive irrigations (used for computing predicted irrigation).
    nk,                  // number of vertical columns of soil cells in the slab.
    nl,                  // number of horizontal layers of soil cells in the slab.
    noitr,               // number of iterations per day, for calling some soil water related functions.
    NumAbscisedLeaves,   // number of leaves, per plant, lost by abscission.
    NumIrrigations,      // number of irrigations.
    NumNitApps,          // number of applications of nitrogen fertilizer.
    NumSheddingTags,     // number of 'box-car' units used for moving values in arrays defining fruit shedding
    // (AbscissionLag, ShedByCarbonStress, ShedByNitrogenStress and ShedByWaterStress).
    NumWaterTableData; // number of water table level input data.

int SoilHorizonNum[maxl],  // the soil horizon number associated with each soil layer in the slab.
    WaterTableLayer;       // number of uppermost soil layer below water table.

double
    AbscissionLag[20],                // the time (in physiological days) from tagging fruiting sites for shedding.
    ActualBollGrowth,                 // total actual growth of seedcotton in bolls, g plant-1 day-1.
    ActualBurrGrowth,                 // total actual growth of burrs in bolls, g plant-1 day-1.
    ActualSquareGrowth,               // total actual growth of squares, g plant-1 day-1.
    ActualStemGrowth,                 // actual growth rate of stems, g plant-1 day-1.
    airdr[9],                         // volumetric water content of soil at "air-dry" for each soil horizon, cm3 cm-3.
    alpha[9],                         // parameter of the Van Genuchten equation.
    AverageLwp,                       // running average of LwpMin + LwpMax for the last 3 days.
    AverageLwpMin,                    // running average of LwpMin for the last 3 days.
    AverageSoilPsi,                   // average soil matric water potential, bars, computed as the weighted average of the root zone.
    vanGenuchtenBeta[9],              // parameter of the Van Genuchten equation.
    BulkDensity[9],                   // bulk density of soil in a horizon, g cm-3.
    BurrNConc,                        // average nitrogen concentration in burrs.
    BurrNitrogen,                     // nitrogen in burrs, g per plant.
    BurrWeightGreenBolls,             // total weight of burrs in green bolls, g plant-1.
    BurrWeightOpenBolls,              // total weight of burrs in open bolls, g per plant.
    CarbonAllocatedForRootGrowth,     // available carbon allocated for root growth, g per plant.
    ClayVolumeFraction[maxl],         // fraction by volume of clay in the soil.
    conmax,                           // the maximum value for non-dimensional hydraulic conductivity
    CottonWeightOpenBolls,            // total weight of seedcotton in open bolls, g per plant.
    CultivationDepth[5],              // depth of cultivation, in cm.
    CumFertilizerN,                   // cumulative amount of inorganic nitrogen applied in fertilizer, mg N per slab.
    CumNetPhotosynth,                 // cumulative sum of photosynthate produced to date, g per plant.
    CumNitrogenUptake,                // cumulative total uptake of nitrogen by plants, mg N per slab.
    CumWaterAdded,                    // cumulative water added to the slab by rainfall or irrigation, mm.
    CumWaterDrained,                  // cumulative water drained out from the slab, mm.
    DayTimeTemp,                      // average day-time temperature, C.
    dclay,                            // aggregation factor for clay in water.
    DeepSoilTemperature,              // boundary soil temperature of deepest layer (K)
    DefoliantAppRate[5],              // rate of defoliant application in pints per acre.
    DepthLastRootLayer,               // the depth to the end of the last layer with roots (cm).
    dsand,                            // aggregation factor for sand in water.
    ElCondSatSoil[20],                // electrical conductivity of saturated soil extract (mmho/cm)
    ElCondSatSoilToday,               // electrical conductivity of saturated extract (mmho/cm) on this day.
    FieldCapacity[maxl],              // volumetric water content of soil at field capacity for each soil layer, cm3 cm-3.
    FoliageTemp[maxk],                // average foliage temperature (oK).
    FreshOrganicNitrogen[maxl][maxk], // N in fresh organic matter in a soil cell, mg cm-3.
    gh2oc[10],                        // input gravimetric soil water content, g g-1, in the soil mechanical impedance table. values have been read from the soil impedance file.
    ginp,                             // ginning percentage of an individual boll.
    Gintot,                           // weighted average ginning percentage of all open bolls.
    GreenBollsLost,                   // cumulative loss of green bolls, due to abscission, g per plant.
    HeatCapacitySoilSolid[maxl],      // heat capacity of the solid phase of the soil.
    HeatCondDrySoil[maxl],            // the heat conductivity of dry soil.
    HumusNitrogen[maxl][maxk],        // N in stable humic fraction material in a soil cells, mg/cm3.
    HumusOrganicMatter[maxl][maxk],   // humus fraction of soil organic matter, mg/cm3.
    impede[10][10],                   // input table of soil impedance to root growth
    InitialTotalSoilWater,            // initial total soil water in the profile, mm.
    IrrigationDepth,                  // depth of predicted irrigation, cm.
    LevelsOfWaterTable[20],           // water table level input data (cm below soil surface).
    LightIntercept,                   // ratio of light interception by plant canopy.
    LwpMax,                           // maximum (dawn) leaf water potential, MPa.
    LwpMin,                           // minimum (noon) leaf water potential, MPa.
    LwpMinX[3],                       // array of values of LwpMin for the last 3 days.
    LwpX[3],                          // array of values of LwpMin + LwpMax for the last 3 days.
    MarginalWaterContent[maxl],       // marginal soil water content (as a function of soil texture) for computing soil heat conductivity.
    MaxIrrigation,                    // maximum amount of applied water in a predicted irrigation
    MaxWaterCapacity[maxl],           // volumetric water content of a soil layer at maximum capacity, before drainage, cm3 cm-3.
    MineralizedOrganicN,              // cumulative amount of mineralized organic N, mgs per slab.
    NetPhotosynthesis,                // net photosynthetic rate, g per plant per day.
    NightTimeTemp,                    // average night-time temperature, C.
    NO3FlowFraction[maxl],            // fraction of nitrate that can move to the next layer.
    PercentDefoliation,               // percentage of leaves abscised as a result of defoliant application.
    PetioleNitrogen,                  // total petiole nitrogen, g per plant.
    PetioleNO3NConc,                  // average nitrate nitrogen concentration in petioles.
    PetioleWeightPreFru[9],           // weight of prefruiting node petioles, g.
    PlantRowLocation,                 // distance of plant row from slab edge, cm.
    PoreSpace[maxl],                  // pore space of soil, volume fraction.
    PotGroAllBolls,                   // sum of potential growth rates of seedcotton in all bolls, g plant-1 day-1.
    PotGroAllBurrs,                   // sum of potential growth rates of burrs in all bolls, g plant-1 day-1.
    PotGroAllLeaves,                  // sum of potential growth rates of all leaves, g plant-1 day-1.
    PotGroAllPetioles,                // sum of potential growth rates of all petioles, g plant-1 day-1.
    PotGroAllRoots,                   // potential growth rate of roots, g plant-1 day-1
    PotGroAllSquares,                 // sum of potential growth rates of all squares, g plant-1 day-1.
    PotGroLeafAreaPreFru[9],          // potentially added area of a prefruiting node leaf, dm2 day-1.
    PotGroLeafWeightPreFru[9],        // potentially added weight of a prefruiting node leaf, g day-1.
    PotGroPetioleWeightPreFru[9],     // potentially added weight of a prefruiting node petiole, g day-1.
    PotGroStem,                       // potential growth rate of stems, g plant-1 day-1.
    RatioImplicit,                    // the ratio for the implicit numerical solution of the water transport equation (used in FLUXI and in SFLUX.
    ReserveC,                         // reserve carbohydrates in leaves, g per plant.
    rlat1[maxl],                      // lateral root length (cm) to the left of the tap root
    rlat2[maxl],                      // lateral root length (cm) to the right of the tap root
    RootImpede[maxl][maxk],           // root mechanical impedance for a soil cell, kg cm-2.
    RootNitrogen,                     // total root nitrogen, g per plant.
    RootWeightLoss,                   // total cumulative weight of sloughed roots, g per plant.
    SandVolumeFraction[maxl],         // fraction by volume of sand plus silt in the soil.
    SaturatedHydCond[9],              // saturated hydraulic conductivity, cm per day.
    SeedNitrogen,                     // total seed nitrogen, g per plant.
    ShedByCarbonStress[20],           // the effect of carbohydrate stress on shedding
    ShedByNitrogenStress[20],         // the effect of nitrogen stress on shedding.
    ShedByWaterStress[20],            // the effect of moisture stress on shedding.
    SitePar[21],                      // array of site specific constant parameters.
    SoilNitrogenLoss,                 // cumulative loss of nitrogen by drainage out of the lowest soil layer, mg per slab.
    SoilPsi[maxl][maxk],              // matric water potential of a soil cell, bars.
    SoilTemp[maxl][maxk],             // hourly soil temperature oK.
    SoilTempDailyAvrg[maxl][maxk],    // daily average soil temperature, oK.
    SquareNitrogen,                   // total nitrogen in the squares, g per plant
    StemNConc,                        // ratio of stem nitrogen to dry matter.
    StemWeight[365],                  // stem weight added at each day, g per plant.
    SumNO3N90,                        // sum of soil nitrate n, 0-90 cm depth, in kg/ha.
    SupplyNH4N,                       // uptake of ammonia N by the plant from the soil, mg N per slab per day.
    SupplyNO3N,                       // uptake of nitrate by the plant from the soil, mg N per slab per day.
    TapRootLength,                    // the length of the taproot, in cm.
    thad[maxl],                       // residual volumetric water content of soil layers (at air-dry condition), cm3 cm-3.
    thetar[maxl],                     // volumetric water content of soil layers at permanent wilting point (-15 bars), cm3 cm-3.
    thetas[9],                        // volumetric saturated water content of soil horizon, cm3 cm-3.
    thts[maxl],                       // saturated volumetric water content of each soil layer, cm3 cm-3.
    TotalActualLeafGrowth,            // actual growth rate of all the leaves, g plant-1 day-1.
    TotalActualPetioleGrowth,         // actual growth rate of all the petioles, g plant-1 day-1.
    TotalPetioleWeight,               // total petiole weight, g per plant.
    TotalRootWeight,                  // total dry weight of the roots, g per plant.
    TotalSoilNh4N,                    // total ammonium in profile, mg N per slab.
    TotalSoilNitrogen,                // total soil nitrogen, mg N per slab.
    TotalSoilNo3N,                    // total nitrate in profile, mg N per slab.
    TotalSoilUreaN,                   // total urea in profile, mg N per slab.
    TotalSoilWater,                   // total water in the soil profile, mm.
    TotalSquareWeight,                // total weight of squares, g per plant
    tstbd[10][10],                    // input bulk density in the impedance table, g cm-3.
    VolNh4NContent[maxl][maxk],       // volumetric ammonium nitrogen content of a soil cell, mg N cm-3.
    VolUreaNContent[maxl][maxk],      // volumetric urea nitrogen content of a soil cell, mg N cm-3.
    VolWaterContent[maxl][maxk];      // volumetric water content of a soil cell, cm3 cm-3.
