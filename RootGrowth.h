#pragma once

#include <string>

using namespace std;

double
PotentialRootGrowth(const int &, const int &, const int &, const double[40][20][3], const double RootAge[40][20]);

tuple<int> ComputeActualRootGrowth(double, const string &, const int &, const int &, const int &, int, const int &,
                                   double[40][20][3], double[40][20]);

typedef struct RootStruct {
    double potential_growth;
    double growth_factor;
    double actual_growth;
    double age;
    double weight[3];
} Root;
