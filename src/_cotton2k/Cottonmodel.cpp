// CottonModel.cpp  "CottonModel" is the simulation module of the Cotton2K model.
// Version 4.0 was written by A. Marani June 2004.
// Compiled by Microsoft Visual C++.Net 2003.
// This file defines the following functions in this file:
//       Class C2KApp:
//    Message map and constructor
//    InitInstance()
//    ExitInstance()
//    GetJobFile()
//    GetProfilesList()
//    RunTheModel()
//    DailySimulation()
//    SimulateThisDay()
//    OnAppAbout()
//       Class CAoutDlg
//
#include <bits/stdint-uintn.h>
#include <cstring>
#include "State.hpp"
#include "global.h"
#include "Cottonmodel.h"
#include "GeneralFunctions.h"
#include "CottonPhenology.h"
#include "DailyClimate.h"
#include "Input.h"
#include "PlantGrowth.h"
#include "PlantNitrogen.h"
#include "SoilNitrogen.h"
#include "SoilProcedures.h"
#include "SoilTemperature.h"


//////////////////////////////////////////////////
void SimulateThisDay(Simulation &sim, uint32_t u)
//     This function executes all the simulation computations in a day. It is called from
//  DailySimulation().   It calls the following functions:
//     ColumnShading(), DayClim(), SoilTemperature(), SoilProcedures(),
//     SoilNitrogen(), SoilSum(), PhysiologicalAge(), Defoliate(), Stress(),
//     GetNetPhotosynthesis(), PlantGrowth(), CottonPhenology(), PlantNitrogen(),
//     CheckDryMatterBal();
//
//     The following global variables are referenced here:
//  iyear, Kday, LastDayWeatherData, LeafAreaIndex.
//
//     The following global variables are set here:
//  isw, Kday.
//
{
    State &state = sim.states[u];
    double rracol[20]; // the relative radiation received by a soil column, as affected by shading by plant canopy.
    // days from emergence
    if (sim.day_emerge <= 0)
        state.kday = 0;
    else
        state.kday = sim.day_start - sim.day_emerge + u + 1;
    if (state.kday < 0)
        state.kday = 0;
    //     The following functions are executed each day (also before emergence).
    ColumnShading(state, rracol, sim.day_emerge, sim.row_space, sim.plant_row_column); // computes light interception and soil shading.
    DayClim(sim, u);                                                                   // computes climate variables for today.
    SoilTemperature(sim, u, rracol);                                                   // executes all modules of soil and canopy temperature.
    SoilProcedures(sim, u);                                                            // executes all other soil processes.
    SoilNitrogen(sim, u);                                                              // computes nitrogen transformations in the soil.
    SoilSum(state, sim.row_space);                                                     // computes totals of water and N in the soil.
                                                                                       //     The following is executed each day after plant emergence:
    if (state.daynum >= sim.day_emerge && isw > 0)
    {
        //     If this day is after emergence, assign to isw the value of 2.
        isw = 2;
        state.day_inc = PhysiologicalAge(state.hours); // physiological days increment for this day. computes physiological age
        Defoliate(sim, u);                                     // effects of defoliants applied.
        Stress(sim, u);                                        // computes water stress factors.
        GetNetPhotosynthesis(sim, u, state.day_length);        // computes net photosynthesis.
        PlantGrowth(sim, u, 3, state.day_length);              // executes all modules of plant growth.
        CottonPhenology(sim, u);                               // executes all modules of plant phenology.
        PlantNitrogen(sim, u);                                 // computes plant nitrogen allocation.
        CheckDryMatterBal(state);            // checks plant dry matter balance.
                                                               //     If the relevant output flag is not zero, compute soil nitrogen balance and soil
                                                               //  nitrogen averages by layer, and write this information to files.
    }
    //     Check if the date to stop simulation has been reached, or if this is the last day
    //  with available weather data. Simulation will also stop when no leaves remain on the plant.
    //
    if (state.daynum >= LastDayWeatherData || (state.kday > 10 && state.leaf_area_index < 0.0002))
        throw SimulationEnd();
}
