#pragma once
#include <string>
#include <tuple>

void WriteStateVariables(bool, const std::string&, const int&, const int&, const double&, const double&, const double[40][20][3], const double[40][20]);
std::tuple<string, int, int, double, double> PlantAdjustments(int, int, const std::string&, const std::string&, const int&, int, double, double, double[40][20][3], double[40][20]);
