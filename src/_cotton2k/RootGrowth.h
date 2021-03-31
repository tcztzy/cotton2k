#pragma once
#include <cinttypes>
#include "Simulation.hpp"

double PotentialRootGrowth(SoilCell[40][20], const int &, const int &);

void ComputeActualRootGrowth(Simulation &, const uint32_t &, double, const int &);
