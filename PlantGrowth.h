#pragma once

#include <string>

using namespace std;

extern "C"
{
    double PhysiologicalAge(const double[24]);
}

tuple<double> Stress(const string &, const double &, const int &);

void GetNetPhotosynthesis(Simulation &, uint32_t, const double &);

tuple<int, double>
PlantGrowth(Simulation &, const uint32_t &, const int &, int, double, const double &, const double &, const double &);

void CheckDryMatterBal(const string &, const string &, const double &);

extern "C"
{
    void Pix();
}

void Defoliate(const string &, const string &, const int &, const int &);
