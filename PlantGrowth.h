#pragma once

#include <string>

using namespace std;

extern "C" {
    double PhysiologicalAge(const double[24]);
}

tuple<double> Stress(const string &, const double &, const int &);

void GetNetPhotosynthesis(const int &, const int &, const int &, const int &, const double &, const double &, const Climstruct[400]);

tuple<int, double>
PlantGrowth(const string &, const string &, const int &, const int &, const int &, const int &, const int &, int,
            const int &, double, const double &, const double &, const double &, double[40][20][3], double[40][20]);

void CheckDryMatterBal(const string &, const string &, const double &);

extern "C" {
    void Pix();
}

void Defoliate(const string &, const string &, const int &, const int &);
