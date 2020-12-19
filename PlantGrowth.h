#pragma once
#include <string>

using namespace std;

double PhysiologicalAge();
tuple<double> Stress(const string&, const double&, const int&);
void GetNetPhotosynthesis(const int&, const int&, const double&);
tuple<int, double> PlantGrowth(const string&, const string&, const int&, const int&, const int&, const int&, int, double, const double&, const double&, const double&, double[40][20][3], double[40][20]);
void CheckDryMatterBal(const string&, const string&, const double&);
void Pix();
void Defoliate(const string&, const string&, const int&, const int&);
