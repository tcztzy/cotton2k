#pragma once
#include <cinttypes>
#include "Simulation.hpp"

void ColumnShading(State &, double[20], double, double, unsigned int);

void SoilTemperature(Simulation &, uint32_t, double[20]);