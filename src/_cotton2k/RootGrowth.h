#pragma once
#include <cinttypes>
#include "Simulation.hpp"

double PotentialRootGrowth(SoilCell[40][20], int, int);

void ComputeActualRootGrowth(Simulation &, uint32_t, double, int);
