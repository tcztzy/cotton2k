#pragma once
#include <string>

using namespace std;

double PhysiologicalAge();
void Stress(const string&, const double&, const int&);
void GetNetPhotosynthesis(const int&, const int&);
tuple<int, double> PlantGrowth(const string&, const string&, const int&, const int&, int, double, const double&);
void CheckDryMatterBal(const string&, const string&);
void Pix();
void Defoliate(const string&, const string&, const int&, const int&);
