#pragma once

#include <string>
#include "Simulation.hpp"

using namespace std;

extern "C"
{
    double PhysiologicalAge(Hour[24]);
}

void PlantGrowth(State &state, double density_factor, double per_plant_area, double row_space, double cultivar_parameters[61], int NumRootAgeGroups, unsigned int day_emerge, unsigned int day_topping, unsigned int first_square, unsigned int plant_row_column);

void CheckDryMatterBal(State &);

void Defoliate(Simulation &, uint32_t);

void LeafWaterPotential(State &, double);

void PotentialLeafGrowth(State &, double, double[61]);

void PotentialFruitGrowth(State &, double[61]);

void DryMatterBalance(State &, double &, double &, double &, double &, double);

void ActualFruitGrowth(State &);

void ActualLeafGrowth(State &);
