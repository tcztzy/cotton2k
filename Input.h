#pragma once

#include <string>
#include <tuple>

std::tuple<int, int, int, int, int, int, int, int, int, int, int, int, double, double, double, double, double>
ReadInput(const std::string &, double[40][20][3], double[40][20], ClimateStruct[400]);

// GettingInput_2
extern "C"
{
    double form(double c0, double d0, double g0);
}
