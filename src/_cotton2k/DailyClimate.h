#pragma once

#include "exceptions.h"
#include "Simulation.hpp"

void DayClim(Simulation &, uint32_t u);

extern "C"
{
    double VaporPressure(double);
    double clearskyemiss(double, double);
}
