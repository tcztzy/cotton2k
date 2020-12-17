#pragma once
#include <string>

using namespace std;

double PotentialRootGrowth(const int&, const double[40][20][3]);
tuple<int> ComputeActualRootGrowth(double, const string&, const int&, const int&, int, double[40][20][3]);
