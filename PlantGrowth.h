#pragma once
#include <string>

using namespace std;

double PhysiologicalAge();
void Stress(const string&, const double&);
void GetNetPhotosynthesis(const int&, const int&);
tuple<double> PlantGrowth(const string&, const string&, const int&, const int&, double);
void CheckDryMatterBal(const string&, const string&);
void Pix();
void Defoliate(const string&, const string&, const int&, const int&);
