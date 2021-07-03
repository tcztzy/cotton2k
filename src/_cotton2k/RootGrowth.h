#pragma once
#include <cinttypes>
#include "Simulation.hpp"

double PotentialRootGrowth(SoilCell[40][20], int, int, double);

void ComputeActualRootGrowth(State &state, double sumpdr, double row_space, double per_plant_area, int NumRootAgeGroups, unsigned int day_emerge, unsigned int plant_row_column);

void RootImpedance(SoilCell[40][20]);
