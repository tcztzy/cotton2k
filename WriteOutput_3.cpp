//  File WriteOutput_3.cpp
//
//   functions in this file:
// outputplt()
// output2()
// output3()
// output4()
// output5()
// output6()
// output7()
// OutputForSoilMaps()
//
#include <numeric>
#include <filesystem>
#include "global.h"
#include "GeneralFunctions.h"

namespace fs = std::filesystem;

void OutputForSoilMaps(State &, int, int, int, int, const string &);

/////////////////////////////////////////////////////////////
void outputplt(Simulation &sim)
//     This function writes the output file that will be used to plot
//  plant charts at the end of the simulation. It is called from DataOutput().
//     The values of the variables to plot are taken from structure Scratch21.
//
//     The following global variables are referenced here:
//  OutIndex, PlantPopulation, RowSpace.
{
//     indexPLT indicates what type of data to write in the .PLT file:
//  indexPLT = 0 : metric (OutIndex(1)=0),  weights and numbers per plant (OutIndex(2)=0)
//  indexPLT = 1 : English (OutIndex(1)=1), weights and numbers per plant (OutIndex(2)=0)
//  indexPLT = 2 : metric (OutIndex(1)=0),  weights and numbers per area  (OutIndex(2)=1)
//  indexPLT = 3 : English (OutIndex(1)=1), weights and numbers per area  (OutIndex(2)=1)
    int indexPLT = OutIndex[1] + 2 * OutIndex[2];
//
    ofstream File25(fs::path("output") / (string(sim.profile_name) + ".PLT"), ios::out);
//     Write heading (first line) in the '*.PLT' file (File25). This line includes:
//  Name of simulation (Profile name)in col 7-26; index in col 30, and
//  date of emergence (day of year) in cols 31-40
    File25 << "      ";
    File25.setf(ios::left);
    File25.width(20);
    File25 << sim.profile_name;
    File25.unsetf(ios::left);
    File25.width(4);
    File25 << indexPLT;
    File25.width(10);
    File25 << sim.day_emerge << endl << endl;
//     Multipliers for output on basis of per area:
    double multi0 = 1; // multiplier for site numbers
    double multi1 = 1; // multiplier for dry weight of plant parts
    if (OutIndex[2] == 1) // per area
    {
        if (OutIndex[1] == 0) // metric units
        {
            multi0 = PlantPopulation * 0.0001; // Convert to sites per m2
            multi1 = PlantPopulation * 0.0001; // Convert to kg per ha
        } else  // English units
        {
            multi0 = PlantPopulation * 0.000405; // convert to 1000s per acre
            multi1 = PlantPopulation * 0.000893; // convert to lbs per acre
        }
    }
//     The structure Scratch21 is extracted for each day of simulation, and
//  data to be used for plotting are written to file '.PLT' (unit 25).
    for (int irec = 0; irec < sim.day_finish - sim.day_start + 1; irec++) {
        State &state = sim.states[irec];
        File25.unsetf(ios::left);
        File25.setf(ios::fixed);
        File25.width(4);   // I4
        File25 << Scratch21[irec].kday;
        File25.width(7);   // F7.1
        File25.precision(1);
        if (OutIndex[1] == 0)   //  metric
        {
            File25.precision(1);
            File25 << state.plant_height; // cm
        } else {
            File25.precision(2);
            File25 << state.plant_height / 2.54; // inches
        }
        File25.width(4);
        File25.precision(0);  // I4
        File25 << Scratch21[irec].mainStemNodes;
        File25.width(6);     // F6.2
        if (OutIndex[1] == 1 && OutIndex[2] == 1)
            File25.precision(1);
        else
            File25.precision(2);
        File25 << Scratch21[irec].leafAreaIndex;
        File25.width(7);     // 4F7.2 or 4F7.1
        File25 << multi0 * state.number_of_squares;
        File25.width(7);
        File25 << multi0 * state.number_of_green_bolls;
        File25.width(7);
        File25 << multi0 * state.number_of_open_bolls;
        File25.width(7);
        File25 << multi0 * state.abscised_fruit_sites;
        File25.precision(1);
        File25.width(7);     // F7.1
        if (OutIndex[2] == 1)
            File25 << Scratch21[irec].lintYield;
        else
            File25 << Scratch21[irec].lintYield * 10000 / PlantPopulation;
        File25.precision(3);
        File25.width(7);     // 5F7.3
        File25 << state.carbon_stress;
        File25.width(7);
        File25 << Scratch21[irec].nStressFruiting;
        File25.width(7);
        File25 << Scratch21[irec].nStressVeg;
        File25.width(7);
        File25 << state.leaf_nitrogen_concentration;
        File25.width(7);
        File25 << state.water_stress;
        File25.width(8);     // F8.3
        File25 << Scratch21[irec].lwpMin;
        File25.precision(2);
        File25.width(8);     // F8.2
        File25 << Scratch21[irec].averageSoilPsi;
        if (OutIndex[2] == 1)
            File25.precision(0); // per area
        else
            File25.precision(1); // per plant
        File25.width(7);     // F7.0 or F7.1
        File25 << multi1 * Scratch21[irec].cumNetPhotosynth;
        if (OutIndex[2] == 1)
            File25.precision(0); // per area
        else
            File25.precision(2); // per plant
        File25.width(7);     // 6F7.0 or 6F7.2
        File25 << multi1 * (Scratch21[irec].totalLeafWeight + Scratch21[irec].reserveC);
        File25.width(7);
        File25 << multi1 * Scratch21[irec].totalStemWeight;
        File25.width(7);
        File25 << multi1 * Scratch21[irec].totalRootWeight;
        File25.width(7);
        File25 << multi1 * Scratch21[irec].totalSquareWeight;
        File25.width(7);
        File25 << multi1 * Scratch21[irec].gbw;
        File25.width(7);
        double cxx = multi1 * (Scratch21[irec].cottonWeightOpenBolls
                               + Scratch21[irec].burrWeightOpenBolls);
        File25 << cxx;  // weight of open bolls including burrs
        File25.width(7);     // 3F7.1 or 3F7.3
        if (OutIndex[1] == 0)   //  metric
        {
            File25.precision(1); // 3F7.1
            File25 << Scratch21[irec].cumEvaporation;  // mm
            File25.width(7);
            File25 << Scratch21[irec].cumTranspiration;
            File25.width(7);
            File25 << Scratch21[irec].cumWaterAdded;
        } else // English
        {
            File25.precision(3); // 3F7.3
            File25 << Scratch21[irec].cumEvaporation / 25.4;  // inches
            File25.width(7);
            File25 << Scratch21[irec].cumTranspiration / 25.4;
            File25.width(7);
            File25 << Scratch21[irec].cumWaterAdded / 25.4;
        }
        if (OutIndex[2] == 1)
            File25.precision(0); // per area
        else
            File25.precision(2); // per plant
        File25.width(7);     // F7.0 or F7.2
        File25 << multi1 * Scratch21[irec].totalPetioleWeight;
//    Compute cumulative uptake and supply of soil N, in kg/ha.
        double supn; // cumulative supply of soil N (from fertilizer and
        //  mineralization, minus loss) in kg/ha or lb/ac.
        double CumUptN; // cumulative uptake of n from soil, kg/ha or lb/ac
        if (OutIndex[1] == 0)   //  metric
        {
            CumUptN = Scratch21[irec].cumNitrogenUptake * 100 / sim.row_space;
            supn = (Scratch21[irec].cumFertilizerN + Scratch21[irec].mineralizedOrganicN -
                    Scratch21[irec].soilNitrogenLoss) * 100 / sim.row_space;
        } else   // English
        {
            CumUptN = Scratch21[irec].cumNitrogenUptake * 89.3 / sim.row_space;
            supn = (Scratch21[irec].cumFertilizerN + Scratch21[irec].mineralizedOrganicN -
                    Scratch21[irec].soilNitrogenLoss) * 89.3 / sim.row_space;
        }
        File25.precision(2);
        File25.width(7);     // 2F7.2
        File25 << CumUptN;
        File25.width(7);
        File25 << supn;
        File25.precision(3);
        File25.width(7);     // F7.3
        File25 << Scratch21[irec].nitrogenStress;
        File25.precision(2);
        File25.width(9);     // 3F9.2
        if (OutIndex[1] == 0)   //  metric
            File25 << Scratch21[irec].sumNO3N90;  // kg / ha
        else
            File25 << Scratch21[irec].sumNO3N90 * 0.893;  // lbs / acre
        File25.width(9);
        File25 << multi1 * Scratch21[irec].cottonWeightGreenBolls;
        File25.width(9);
        File25 << multi1 * Scratch21[irec].cottonWeightOpenBolls;
        File25.precision(3);
        File25.width(7);     // 6F7.3
        File25 << state.petiole_nitrogen_concentration;
        File25.width(7);
        File25 << Scratch21[irec].petioleNO3NConc;
        File25.width(7);
        File25 << Scratch21[irec].rootNConc;
        File25.width(7);
        File25 << Scratch21[irec].stemNConc;
        File25.width(7);
        File25 << Scratch21[irec].burrNConc;
        File25.width(7);
        File25 << Scratch21[irec].seedNConc;
        File25 << endl;
    }
}

///////////////////////
void output2(Simulation &sim)
//     This function is optionally called from DataOutput().  It writes output of stress
//  factors to file F01, when the output flag OutIndex(6) is non-zero.
//     The values of the variables to plot are taken from structure Scratch21.
{
//     Write header lines.
    ofstream File46(fs::path("output") / (string(sim.profile_name) + ".F01"), ios::app);
    File46.unsetf(ios::left);
    File46 << endl << "                            Stress factors output" << endl;
    File46 << "                %                            Leaf           Veg.    Leaf    Soil" << endl;
    File46 << "     Date     Light  C Stres --N-Stress--    N      Water   Water  Water   Water" << endl;
    File46 << "              Interc.         Fruit   Veg.   conc.  stress  stress  -- Pot. MPa --" << endl;
//     Extract data from Scratch21 and write line for each day.
    for (int irec = 0; irec < sim.day_finish - sim.day_start + 1; irec++) {
        if (Scratch21[irec].kday <= 0)
            continue;
//
        File46.unsetf(ios::left);
        File46.width(12);
        File46 << sim.states[irec].date;
        File46.setf(ios::fixed);
        File46.precision(1);
        File46.width(6);
        File46 << Scratch21[irec].lightIntercept * 100;
        File46.precision(3);
        File46.width(8);
        File46 << sim.states[irec].carbon_stress;
        File46.width(8);
        File46 << Scratch21[irec].nStressFruiting;
        File46.width(8);
        File46 << Scratch21[irec].nStressVeg;
        File46.width(8);
        File46 << sim.states[irec].leaf_nitrogen_concentration;
        File46.width(8);
        File46 << sim.states[irec].water_stress;
        File46.width(8);
        File46 << sim.states[irec].water_stress_stem;
        File46.width(8);
        File46 << Scratch21[irec].lwpMin;
        File46.width(8);
        File46 << Scratch21[irec].averageSoilPsi * 0.1;  // converted to MPa
        File46 << endl;
    }
}

///////////////////////
void output3(Simulation &sim)
//     This procedure is optionally called from DataOutput(). It writes output of weights of 
//  dry matter of plant parts to file F01. Data are extracted from structure Scratch21. 
//     Global variables referenced:    OutIndex, PlantPopulation
{
//  Write header lines.
    double multi;  // multiplier for some data
    ofstream File46(fs::path("output") / (string(sim.profile_name) + ".F01"), ios::app);
    File46.unsetf(ios::left);
    if (OutIndex[2] == 0) {
        multi = 1;
        File46 << endl << "                          Weights in grams dry matter per plant " << endl;
    } else if (OutIndex[1] == 0) {
        multi = PlantPopulation * 0.001;
        File46 << endl << "                          Weights in kgs dry matter per hectare " << endl;
    } else {
        multi = PlantPopulation * 0.000893;
        File46 << endl << "                          Weights in lbs dry matter per acre " << endl;
    }
    File46 << "     Date      Photosynth.  Plant  Leaf  Peti. Stem  Root Squares -- Bolls --  Deadwt" << endl;
    File46 << "               net    cumul.                                      Green  Open" << endl << endl;
//   Write data for each day.
    for (int irec = 0; irec < sim.day_finish - sim.day_start + 1; irec++) {
        if (Scratch21[irec].kday <= 0)
            continue;
        File46.unsetf(ios::left);
        File46.width(12);
        File46 << sim.states[irec].date;
        File46.setf(ios::fixed);
        if (OutIndex[2] == 0)
            File46.precision(2);
        else
            File46.precision(0);
        File46.width(7);
        File46 << multi * Scratch21[irec].netPhotosynthesis;
        File46.width(7);
        File46 << multi * Scratch21[irec].cumNetPhotosynth;
        File46.width(7);
        File46 << multi * Scratch21[irec].plantWeight;
        File46.width(6);
        File46 << multi * Scratch21[irec].totalLeafWeight + Scratch21[irec].reserveC;
        File46.width(6);
        File46 << multi * Scratch21[irec].totalPetioleWeight;
        File46.width(6);
        File46 << multi * Scratch21[irec].totalStemWeight;
        File46.width(6);
        File46 << multi * Scratch21[irec].totalRootWeight;
        File46.width(6);
        File46 << multi * Scratch21[irec].totalSquareWeight;
        File46.width(7);
        File46 << multi * Scratch21[irec].gbw;
        File46.width(7);
        double cxx = multi * (Scratch21[irec].cottonWeightOpenBolls
                              + Scratch21[irec].burrWeightOpenBolls);
        File46 << cxx; // weight of open bolls including burrs
        File46.width(7);
        File46 << multi * Scratch21[irec].deadwt << endl;
    }
}

///////////////////////
void output4(Simulation &sim)
//     This procedure is optionally called from DataOutput(). It writes output of water
//  and evapotranspiration data to file F01, when the output flag OutIndex(4) is non-zero.
//     Data are extracted from structure Scratch21. 
//     Global variables referenced:         OutIndex
{
//     Write header lines.
    double multi;  // multiplier for some data
    ofstream File46(fs::path("output") / (string(sim.profile_name) + ".F01"), ios::app);
    File46.unsetf(ios::left);
    File46 << endl << "                          Water and ET variables" << endl;
    File46 << "    Date      Potential   Net    ---- Actual Cumulative --  Total      " << endl;
    File46 << "               ES    EP   Rad.    ES     EP   Drain  Irrig. Water  H2OBAL" << endl;
    if (OutIndex[1] == 0) {
        multi = 1;
        File46 << "              --- mm ---  MJm-2   ---------------  mm  ------------------ " << endl << endl;
    } else {
        multi = 1 / 25.4; // convert from mm to inches
        File46 << "              - inches -   Ly.    -------------- inches ----------------- " << endl << endl;
    }
//
    for (int irec = 0; irec < sim.day_finish - sim.day_start + 1; irec++) {
        File46.unsetf(ios::left);
        File46.width(12);
        File46 << sim.states[irec].date;
        File46.setf(ios::fixed);
        File46.precision(2);
        File46.width(6);
        File46 << Scratch21[irec].es * multi;
        File46.width(6);
        File46 << Scratch21[irec].ep * multi;
        File46.precision(1);
        File46.width(7);
        File46.precision(1);
        if (OutIndex[1] == 0) {
            File46 << sim.states[irec].net_radiation * 0.0036;
            File46.precision(1);
        } else {
            File46 << sim.states[irec].net_radiation * 0.0036 * 23.884;
            File46.precision(3);
        }
        File46.width(6);
        File46 << Scratch21[irec].cumEvaporation * multi;
        File46.width(7);
        File46 << Scratch21[irec].cumTranspiration * multi;
        File46.width(7);
        File46 << Scratch21[irec].cumWaterDrained * multi;
        File46.width(7);
        File46 << Scratch21[irec].cumWaterAdded * multi;
        File46.width(7);
        File46 << Scratch21[irec].totalSoilWater * multi;
        File46.width(9);
        File46 << Scratch21[irec].h2obal * multi;
        File46 << endl;
    }
}

///////////////////////
void output5(Simulation &sim)
//     This procedure is optionally called from DataOutput() when the output flag
// OutIndex[5] is non-zero. It writes output of weather variables to file F01.
//     OutIndex[1] determines if the units are in  metric or English.
//
{
//     Write header lines.
    ofstream File46(fs::path("output") / (string(sim.profile_name) + ".F01"), ios::app);
    File46.unsetf(ios::left);
    if (OutIndex[1] == 0) {
        File46 << endl << "                    Weather variables in metric units" << endl;
        File46 << "               -- Temperature, degrees C --" << endl;
        File46 << "     Date      Max   Min   Avg   Avg   Avg   Solar    Rain  +Irrig  Runoff  Wind" << endl;
        File46 << "                          Daily  Day  Night (MJm-2)   (mm)   (mm)    (mm)   (km)" << endl << endl;
    } else {
        File46 << endl << "                    Weather variables in English units" << endl;
        File46 << "               -- Temperature, degrees F --" << endl;
        File46 << "     Date      Max   Min   Avg   Avg   Avg   Solar    Rain  +Irrig  Runoff  Wind" << endl;
        File46 << "                          Daily  Day  Night  (Ly)     (in)   (in)    (in)  (miles)" << endl << endl;
    }
//
    for (int irec = 0; irec < sim.day_finish - sim.day_start + 1; irec++) {
        File46.unsetf(ios::left);
        File46.width(12);
        File46 << sim.states[irec].date;
        File46.setf(ios::fixed);
        File46.precision(1);
        if (OutIndex[1] == 0) {
            File46.width(6);
            File46 << sim.climate[irec].Tmax;
            File46.width(6);
            File46 << sim.climate[irec].Tmin;
            File46.width(6);
            File46 << Scratch21[irec].avrgDailyTemp;
            File46.width(6);
            File46 << Scratch21[irec].dayTimeTemp;
            File46.width(6);
            File46 << Scratch21[irec].nightTimeTemp;
            File46.precision(2);
            File46.width(8);
            File46 << sim.climate[irec].Rad / 23.884;
            File46.width(8);
            File46 << sim.climate[irec].Rain;
            File46.width(8);
            File46 << Scratch21[irec].amitri;
            File46.width(7);
            File46 << sim.states[irec].runoff;
            File46.width(8);
            File46 << sim.climate[irec].Wind << endl;
        } else {
            File46.width(6);
            File46 << sim.climate[irec].Tmax * 1.8 + 32;
            File46.width(6);
            File46 << sim.climate[irec].Tmin * 1.8 + 32;
            File46.width(6);
            File46 << Scratch21[irec].avrgDailyTemp * 1.8 + 32;
            File46.width(6);
            File46 << Scratch21[irec].dayTimeTemp * 1.8 + 32;
            File46.width(6);
            File46 << Scratch21[irec].nightTimeTemp * 1.8 + 32;
            File46.precision(2);
            File46.width(8);
            File46 << sim.climate[irec].Rad;
            File46.width(8);
            File46 << sim.climate[irec].Rain / 25.4;
            File46.width(8);
            File46 << Scratch21[irec].amitri / 25.4;
            File46.width(7);
            File46 << sim.states[irec].runoff / 25.4;
            File46.width(8);
            File46 << sim.climate[irec].Wind / 1.609 << endl;
        }
    }
}

///////////////////////
void output6(const string &ProfileName)
//     This procedure is always called from DataOutput().
//     It writes output of final yields to file F01 and file S01.
{
    ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::app);
    File46 << endl << "      Lint Yield:  kgs / ha  lbs / acre   bales / acre" << endl;
    File46.unsetf(ios::left);
    File46.setf(ios::fixed);
    File46.precision(1);
    File46.width(26);
    File46 << LintYield;
    File46.width(12);
    File46 << LintYield * 0.893;
    File46.precision(2);
    File46.width(13);
    File46 << LintYield * 0.893 / 500 << endl;
//
    ofstream File22(fs::path("output") / (ProfileName + ".S01"), ios::app);
    File22 << endl << "      Lint Yield:  kgs / ha  lbs / acre   bales / acre" << endl;
    File22.unsetf(ios::left);
    File22.setf(ios::fixed);
    File22.precision(1);
    File22.width(26);
    File22 << LintYield;
    File22.width(12);
    File22 << LintYield * 0.893;
    File22.precision(2);
    File22.width(13);
    File22 << LintYield * 0.893 / 500 << endl;
}

///////////////////////
void output7(Simulation &sim)
//     This function is called by DataOutput(). It writes soil map output by calling
//  OutputForSoilMaps(). If the day of year is between DayStartSoilMaps and DayStopSoilMaps,
//  the following will be executed at SoilMapFreq day intervals. It will also be executed 
//  if the simulation is ended.
//
//     The following global variables are referenced:
//       DayStartSoilMaps, DayStopSoilMaps, SoilMapFreq,
{
    for (int nDaynum = sim.day_start_soil_maps; nDaynum <= sim.day_stop_soil_maps; nDaynum++) {
        int idum = nDaynum - sim.day_start_soil_maps;
        if ((idum % SoilMapFreq) == 0 || nDaynum >= sim.day_finish) {
            int irec = nDaynum - sim.day_start;
            State &state = sim.states[irec];
//     If the output flag (OutIndex(8)) indicates that this output is
//  requested - produce soil slab maps for VolWaterContent.
            if (OutIndex[8] > 0)
                OutputForSoilMaps(state, irec, 2, nDaynum, sim.year, sim.profile_name);
//     If the output flag (OutIndex(9)) indicates that this output is
//  requested - produce soil slab maps for root weight per cell and for RootWtCapblUptake.
            if (OutIndex[9] > 0) {
                OutputForSoilMaps(state, irec, 3, nDaynum, sim.year, sim.profile_name);
                OutputForSoilMaps(state, irec, 6, nDaynum, sim.year, sim.profile_name);
            }
//     If the output flag (OutIndex(10)) indicates that this output is
//  requested - produce soil slab maps for VolNo3NContent and for VolNh4NContent.
            if (OutIndex[10] > 0) {
                OutputForSoilMaps(state, irec, 1, nDaynum, sim.year, sim.profile_name);
                OutputForSoilMaps(state, irec, 7, nDaynum, sim.year, sim.profile_name);
            }
//     If the output flag (OutIndex(11)) indicates that this output is
//  requested - produce soil slab maps for soil water potential.
            if (OutIndex[11] > 0)
                OutputForSoilMaps(state, irec, 4, nDaynum, sim.year, sim.profile_name);
//     If the output flag (OutIndex(12)) indicates that this output is
//  requested - produce soil slab maps for soil temperature.
            if (OutIndex[12] > 0)
                OutputForSoilMaps(state, irec, 5, nDaynum, sim.year, sim.profile_name);
        }
    }
}

//////////////////////////////////////////////////////////////
void OutputForSoilMaps(State &state, int irec, int igo, int nday, int year, const string &ProfileName)
//     This function plots the soil slab and the array variables in each cell. It is called from
//  function output7().
//
//     The following global variables are referenced here:
//       Date, dl, FieldCapacity, nl, PoreSpace, thad, thetar, wk.
//     The following variables are used as arguments:
//       irec - record number for extracting the data from structure Scratch21.
//       igo - indicates which variable is to be plotted: (1) - nitrate content;
//             (2) - volumetric water content;   (3) - total root weight;
//             (4) - soil water potential;       (5) - soil temperature;
//             (6) - roots capable of uptake;    (7) - ammonium content.
//       nDay - day of year for these data.
{
//     The following constant variables are used:
    string ka[12] = {"  ", " 0", " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8",
                     " 9", " *"}; // character symbols used in plots.
    double psisca[11] = {-50, -15, -8, -4, -2, -1, -0.6, -0.4,
                         -0.2, -0.1, 0}; // range of values for plotting soil water potential.
    double roocup[11] = {0, 0.0005, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
                         0.70, 0.80, 1}; // range of values for plotting roots capable of uptake.
    double roosca[11] = {0, 0.001, 0.2, 0.4, 0.6, 0.8, 1,
                         1.2, 1.4, 1.6, 2}; // range of values for plotting total root weight.
    double tempsca[11] = {0, 18, 20, 22, 24, 26, 28, 30, 32,
                          34, 36}; // range of values for plotting soil temperature.
    double vnhsca[11] = {0, 0.00005, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030,
                         0.040, 0.050, 0.10}; // range of values for plotting soil ammonium.
    double vnosca[11] = {0, 0.0001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10,
                         0.20}; // range of values for plotting soil nitrate.
    double Array[maxl][maxk];
    string khar[maxl][maxk]; // array of symbols plotted.
    double range[11]; // array of range values used for plotting.
    string tl1; // keyword in first line of plot title.
//     Define the array of values to plot, titles, units, and range of
//  values used for plotting soil nitrate.
    if (igo == 1) {
        for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
                Array[l][k] = Scratch21[irec].volNo3NContent[l][k];
        tl1 = "NITRATE N";
        for (int i = 0; i < 11; i++)
            range[i] = vnosca[i];
    }
//     Define titles, units, and range of values used for plotting
//  volumetric soil water content.
//     Compute average water content at air-dry (wad), wilting point
//  (wwp), field capacity (wfc), and saturated (wsa).
    else if (igo == 2) {
        for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
                Array[l][k] = Scratch21[irec].volWaterContent[l][k];
        double wad = 0; // average value of air-dry water content.
        double wfc = 0; // average value of field-capacity water content.
        double wsa = 0; // average value of saturated water content.
        double wwp = 0; // average value of wilting point water content.
        for (int l = 0; l < nl; l++) {
            wad += thad[l];
            wwp += thetar[l];
            wfc += FieldCapacity[l];
            wsa += PoreSpace[l];
        }
        wad = wad / nl;
        wwp = wwp / nl;
        wfc = wfc / nl;
        wsa = wsa / nl - 0.005;
//    Set scaling
        double capsca[11]; // range of values for plotting volumetric water content.
        capsca[0] = 0;
        capsca[1] = 0.001;
        capsca[2] = wad;
        capsca[3] = wwp;
        for (int i = 4; i < 9; i++)
            capsca[i] = capsca[i - 1] + (wfc - wwp) / 5;
        capsca[9] = wsa;
        capsca[10] = 0.99;
        tl1 = "WATER    ";
        for (int i = 0; i < 11; i++)
            range[i] = capsca[i];
    }
//     Define titles, units, and range of values used for plotting
//  total root weight.
    else if (igo == 3) {
        for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
                Array[l][k] = accumulate(state.root[l][k].weight, state.root[l][k].weight + 3, double(0));
        tl1 = "ROOT TOT ";
        for (int i = 0; i < 11; i++)
            range[i] = roosca[i];
    }
//     Define titles, units, and range of values used for plotting
//  soil water potential.
    else if (igo == 4) {
        for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
                Array[l][k] = Scratch21[irec].soilPsi[l][k];
        tl1 = "PSIS     ";
        for (int i = 0; i < 11; i++)
            range[i] = psisca[i];
    }
//     Define titles, units, and range of values used for plotting
//  average soil temperature.
    else if (igo == 5) {
        for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
                Array[l][k] = Scratch21[irec].soilTempDailyAvrg[l][k];
        tl1 = "SOILT    ";
        for (int i = 0; i < 11; i++)
            range[i] = tempsca[i];
    }
//     Define titles, units, and range of values used for plotting
//  root weight capable of uptake.
    else if (igo == 6) {
        for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
                Array[l][k] = state.root[l][k].weight_capable_uptake;
        tl1 = "ROOT UPT ";
        for (int i = 0; i < 11; i++)
            range[i] = roocup[i];
    }
//     Define titles, units, and range of values used for plotting
//  soil ammonia.
    else if (igo == 7) {
        for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
                Array[l][k] = Scratch21[irec].volNh4NContent[l][k];
        tl1 = "AMMONIA N";
        for (int i = 0; i < 11; i++)
            range[i] = vnhsca[i];
    }
//
    for (int k = 0; k < maxk; k++)
        for (int l = 0; l < maxl; l++)
            khar[l][k] = "  ";
    for (int k = 0; k < maxk; k++)
        for (int l = 0; l < maxl; l++) {
            double araylk; // array element for cell (l,k).
//     Assign the value for this cell to araylk. 
            if (igo == 5) // For soil temperature - convert from K to C.
                araylk = Array[l][k] - 273.161;
            else if (igo == 3 || igo == 6)// For root weight capable of uptake or total
                // roots - convert from g to mg mass per cm3 soil.
                araylk = Array[l][k] * 1000 / (dl[l] * wk[k]);
            else
                araylk = Array[l][k];
//     Determine to which range of values belongs araylk, and assign accordingly
//  the appropriate symbol (ka) to the element of khar for this cell.
            int i;
            for (i = 0; i < 11; i++) {
                if (araylk <= range[i])
                    break;
            }
            khar[l][k] = ka[i];
        } //
//     Write the data of the soil slab to file SMP, first line includes title
//  (keyword) and date.
    ofstream File23(fs::path("output") / (ProfileName + ".SMP"), ios::app);
    File23.width(9);
    File23 << endl << tl1 << "           " << DoyToDate(nday, year) << endl << endl;
//     Next line includes range data.
    File23.setf(ios::fixed);
    for (int i = 0; i < 11; i++) {
        File23.width(7);
        File23.precision(3);
        File23 << range[i];
    }
    File23 << endl;
//     Loop of all layers to print the slab diagram.
    for (int l = 0; l < maxl; l++) {
        File23.width(3);
        File23 << l + 1;
        File23 << "  ";
        for (int k = 0; k < maxk; k++)
            File23 << khar[l][k];
        File23 << endl;
    }
}
