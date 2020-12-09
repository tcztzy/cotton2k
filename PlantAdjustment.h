#pragma once
#include <string>
#include <tuple>

void WriteStateVariables(bool, const std::string&, const int&, const int&, const double&);
std::tuple<string, int, int, double> PlantAdjustments(int, int, const std::string&, const std::string&, const int&, int, double);
