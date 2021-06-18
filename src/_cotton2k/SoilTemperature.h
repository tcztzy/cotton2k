#pragma once
#include <cinttypes>
#include "Simulation.hpp"

void SoilTemperatureInit(Simulation &);

void PredictEmergence(Simulation &, unsigned int, int);

void SoilHeatFlux(State &, double, int, int, int, int, double);

void CanopyBalance(int, int, double, double, double, double, double, double, double, double &, const int &);

void SoilSurfaceBalance(State &, int, int, double, double, double, double, double, double &, double &, double &, double, double);

double SensibleHeatTransfer(double, double, double, double);
