#pragma once

#include <string>
#include <tuple>
#include "cotton2k.h"

std::tuple<double> DayClim(Simulation &, uint32_t u);

extern "C"
{
    double VaporPressure(double);
    double clearskyemiss(double, double);
}
