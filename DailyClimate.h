#pragma once

#include <string>
#include <tuple>
#include "cotton2k.h"

void DayClim(Simulation &, uint32_t u);

extern "C"
{
    double VaporPressure(double);
    double clearskyemiss(double, double);
}
