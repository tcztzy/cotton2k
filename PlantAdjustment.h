#pragma once

#include <string>
#include <tuple>

void WriteStateVariables(bool, const std::string &, const int &, const int &, const int &, const int &, const int &,
                         const double &, const double &, const double &, const double &, const double[40][20][3],
                         const double[40][20], const ClimateStruct[400]);
