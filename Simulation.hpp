#ifndef SIMULATION_TYPE
#define SIMULATION_TYPE
#include "State.h"
#include "Climate.h"
#include "Irrigation.h"
typedef struct Simulation
{
    const char *profile_name;               // name of input file with profile data (without the extension ".PRO")
    int year;                               // Simulation start year
    unsigned int day_emerge;                // Date of emergence (DOY).
    unsigned int day_start;                 // Date (DOY) to start simulation.
    unsigned int day_finish;                // Date (DOY) to finish simulation.
    unsigned int day_plant;                 // Date (DOY) of planting.
    unsigned int day_start_soil_maps;       // Date (DOY) to start soil slab maps output.
    unsigned int day_stop_soil_maps;        // Date (DOY) to stop soil slab maps output.
    double latitude;                        // degree
    double longitude;                       // degree
    double elevation;                       // meter
    double row_space;                       // average row spacing, cm.
    unsigned int num_curve;                 // number of input soil-moisture curves in the impedance table.
    unsigned int first_bloom;               // Date (DOY) of first bloom.
    unsigned int first_square;              // Date of first square (DOY), if no squares have been formed, FirstSquare = 0.
    unsigned int plant_row_column;          // column number to the left of plant row location.
    ClimateStruct climate[400];             // structure containing the following daily weather data:
                                            // int nDay =    day of year.
                                            // double Rad =  daily global radiation, in langleys.
                                            // double Tmax = maximum daily temperature, C.
                                            // double Tmin = minimum daily temperature, C.
                                            // double Tdew = dew point temperature, C.
                                            // double Rain = daily rainfall, mm.
                                            // double Wind = daily wind run, km.
    Irrigation irrigation[150];
    State *states;
} Simulation;
#endif
