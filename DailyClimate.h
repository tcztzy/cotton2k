#pragma once

#include <string>
#include <tuple>
#include "cotton2k.h"

std::tuple<double>
DayClim(const std::string &, const std::string &, const int &, const int &, const int &, const int &, const double &,
        const double &, Climstruct[400]);

extern "C"
{
    double VaporPressure(double);
    double clearskyemiss(double, double);
}
