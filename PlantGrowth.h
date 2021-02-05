#pragma once

#include <string>

using namespace std;

extern "C"
{
    double PhysiologicalAge(Hour[24]);
    void Pix();
}

void Stress(Simulation &, unsigned int);

void GetNetPhotosynthesis(Simulation &, uint32_t, const double &);

void PlantGrowth(Simulation &, const uint32_t &, const int &, const double &);

void CheckDryMatterBal(const string &, const string &, const double &);

void Defoliate(Simulation &, uint32_t);
