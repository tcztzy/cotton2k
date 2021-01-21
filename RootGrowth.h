#pragma once

#include <string>

using namespace std;

double
PotentialRootGrowth(const int &, const int &, const int &, Root[40][20], const double RootAge[40][20]);

tuple<int> ComputeActualRootGrowth(double, Simulation &, const int &, const int &, const int &, int, const int &, double[40][20]);
