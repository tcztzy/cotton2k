#pragma once

#include "Simulation.h"

using namespace std;

double PotentialRootGrowth(SoilCell[40][20], const int &, const int &, const int &);

void ComputeActualRootGrowth(Simulation &, const uint32_t &, double, const int &);
