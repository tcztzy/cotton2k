#pragma once

#include <string>
#include "Simulation.hpp"

using namespace std;

extern "C"
{
    double PhysiologicalAge(Hour[24]);
}

void Stress(State &, double);

void GetNetPhotosynthesis(Simulation &, uint32_t, const double &);

void PlantGrowth(State &state, double density_factor, double per_plant_area, double row_space, double cultivar_parameters[61], int NumRootAgeGroups, unsigned int day_emerge, unsigned int day_topping, unsigned int first_square, unsigned int plant_row_column);

void CheckDryMatterBal(State &);

void Defoliate(Simulation &, uint32_t);
