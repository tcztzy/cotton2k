#pragma once

#include "exceptions.h"
#include "Simulation.h"

void DayClim(Simulation &, uint32_t u);

extern "C"
{
    double VaporPressure(double);
    double clearskyemiss(double, double);
}
