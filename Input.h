#pragma once

#include "Simulation.h"

Simulation ReadInput(const char *);

// GettingInput_2
extern "C"
{
    double form(double c0, double d0, double g0);
}
