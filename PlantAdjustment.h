#pragma once
#include <string>
#include <tuple>

void WriteStateVariables(bool, const std::string&, const int&);
std::tuple<string, int> PlantAdjustments(int, int, const std::string&, const std::string&, const int&);
