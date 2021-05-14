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


void DailySimulation(Simulation &sim)
//     This function controls the dynamic phase of the simulation.
//     It calls the functions:
//        SimulateThisDay()
//
//     The following global variable are referenced:   Kday.
//
{
    try
    {
        initialize_state0(sim.states[0], sim.day_start);
        SimulateThisDay(sim, 0);
        for (int i = 1; i < sim.day_finish - sim.day_start + 1; i++)
        {
            memcpy(&sim.states[i], &sim.states[i - 1], sizeof(State));
            sim.states[i].daynum++;
            SimulateThisDay(sim, i);
        }
    }
    catch (SimulationEnd)
    {
    }
}

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


void initialize_state0(State &state0, uint32_t day_start) {
    state0.daynum = day_start;
    state0.lint_yield = 0;
    state0.soil.number_of_layers_with_root = 7;
    state0.plant_height = 4.0;
    state0.plant_weight = 0;
    state0.stem_weight = 0.2;
    state0.green_bolls_weight = 0;
    state0.bloom_weight_loss = 0;
    state0.abscised_fruit_sites = 0;
    state0.abscised_leaf_weight = 0;
    state0.cumulative_nitrogen_loss = 0;
    state0.cumulative_transpiration = 0;
    state0.cumulative_evaporation = 0;
    state0.applied_water = 0;
    state0.water_stress = 1;
    state0.water_stress_stem = 1;
    state0.carbon_stress = 1;
    state0.extra_carbon = 0;
    state0.leaf_area_index = 0.001;
    state0.leaf_area = 0;
    state0.leaf_weight = 0.20;
    state0.leaf_nitrogen = 0.0112;
    state0.number_of_vegetative_branches = 1;
    state0.number_of_squares = 0;
    state0.number_of_green_bolls = 0;
    state0.number_of_open_bolls = 0;
    state0.nitrogen_stress = 1;
    state0.nitrogen_stress_vegetative = 1;
    state0.nitrogen_stress_fruiting = 1;
    state0.nitrogen_stress_root = 1;
    state0.total_required_nitrogen = 0;
    state0.leaf_nitrogen_concentration = .056;
    state0.petiole_nitrogen_concentration = 0;
    state0.seed_nitrogen_concentration = 0;
    state0.root_nitrogen_concentration = .026;
    state0.square_nitrogen_concentration = 0;
    state0.square_nitrogen = 0;
    state0.stem_nitrogen = 0.0072;
    state0.fruit_growth_ratio = 1;
    state0.ginning_percent = 0.35;
    state0.number_of_pre_fruiting_nodes = 1;
    for (int i = 0; i < 9; i++) {
        state0.age_of_pre_fruiting_nodes[i] = 0;
        state0.leaf_area_pre_fruiting[i] = 0;
        state0.leaf_weight_pre_fruiting[i] = 0;
    }
    for (int k = 0; k < 3; k++)
    {
        state0.vegetative_branches[k].number_of_fruiting_branches = 0;
        for (int l = 0; l < 30; l++)
        {
            state0.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes = 0;
            state0.vegetative_branches[k].fruiting_branches[l].delay_for_new_node = 0;
            state0.vegetative_branches[k].fruiting_branches[l].main_stem_leaf = {0, 0, 0, 0, 0, 0};
            for (int m = 0; m < 5; m++)
                state0.vegetative_branches[k].fruiting_branches[l].nodes[m] = {0, 0, 0, Stage::NotYetFormed, {0, 0, 0, 0}, {0, 0}, {0, 0, 0}, {0, 0}, {0, 0}};
        }
    }
}
