#pragma once

#include <string>
#include "Simulation.hpp"

using namespace std;

extern "C"
{
    double PhysiologicalAge(Hour[24]);
}

void Stress(Simulation &, unsigned int);

void GetNetPhotosynthesis(Simulation &, uint32_t, const double &);

void PlantGrowth(Simulation &, const uint32_t &, const int &, const double &);

void CheckDryMatterBal(State &);

void Defoliate(Simulation &, uint32_t);