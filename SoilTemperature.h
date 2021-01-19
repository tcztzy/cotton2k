#pragma once

#include <string>

void ColumnShading(Simulation &, const int &, const int &, const double &, double[20]);

tuple<int>
SoilTemperature(Simulation &, const std::string &, const int &, const int &, int, const int &, const int &, const int &, const int &,
                const int &, const int &, const double &, const double &, const double &, double[20], const ClimateStruct[400]);
