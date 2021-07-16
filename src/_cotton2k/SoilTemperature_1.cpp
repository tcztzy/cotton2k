// File SoilTemperature_1.cpp
//   List of functions in this file:
//       ColumnShading()
//       SoilTemperature()
//       SoilTemperatureInit()
//
#include <cmath>
#include "global.h"
#include "exceptions.h"
#include "Simulation.hpp"

using namespace std;


////////////////////////////////////////////////////////////////////////
void SoilTemperatureInit(Simulation &sim)
//     This function is called from SoilTemperature() at the start of the simulation. It sets
//  initial values to soil and canopy temperatures.
//     The following arguments are set in this function:
//  jt1, jt2 - input of start and stop of output of soil temperatures.
//
//     The following global variables are referenced here:
//  Clim (structure), nl, SitePar.
//     The following global variables are set here:
//  SoilTemp.
{
    //     Compute initial values of soil temperature: It is assumed that at the start of simulation
    //  the temperature of the first soil layer (upper boundary) is equal to the average air temperature
    //  of the previous five days (if climate data not available - start from first climate data).
    //     NOTE: For a good simulation of soil temperature, it is recommended to start simulation at
    //  least 10 days before planting date. This means that climate data should be available for
    //  this period. This is especially important if emergence date has to be simulated.
    int idd = 0;     // number of days minus 4 from start of simulation.
    double tsi1 = 0; // Upper boundary (surface layer) initial soil temperature, C.
    for (int i = 0; i < 5; i++)
        tsi1 += sim.climate[i].Tmax + sim.climate[i].Tmin;
    tsi1 = tsi1 / 10;
    //     The temperature of the last soil layer (lower boundary) is computed as a sinusoidal function
    //  of day of year, with site-specific parameters.
    sim.states[0].deep_soil_temperature = SitePar[9] + SitePar[10] * sin(2 * pi * (sim.day_start - SitePar[11]) / 365) + 273.161;
    //     SoilTemp is assigned to all columns, converted to degrees K.
    tsi1 += 273.161;
    for (int l = 0; l < nl; l++)
    {
        //     The temperatures of the other soil layers are linearly interpolated.
        //  tsi = computed initial soil temperature, C, for each layer
        double tsi = ((nl - l - 1) * tsi1 + l * sim.states[0].deep_soil_temperature) / (nl - 1);
        for (int k = 0; k < nk; k++)
            SoilTemp[l][k] = tsi;
    }
}
