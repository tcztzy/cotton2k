#pragma once

#include <string>

void ColumnShading(const int &, const int &, const double &, double[20]);

tuple<int>
SoilTemperature(const std::string &, const int &, const int &, int, const int &, const int &, const int &, const int &,
                const int &, const int &, const double &, const double &, const double &, double[20], const ClimateStruct[400]);
