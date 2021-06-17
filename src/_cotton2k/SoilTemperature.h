#pragma once
#include <cinttypes>
#include "Simulation.hpp"

void SoilTemperatureInit(Simulation &);

void PredictEmergence(Simulation &, unsigned int, int);

void SoilHeatFlux(State &, double, int, int, int, int, double);

void EnergyBalance(Simulation &, uint32_t, int, int, double, double, const double &, double[20]);
