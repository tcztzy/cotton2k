#pragma once
#include <cinttypes>
#include "Simulation.hpp"

void SoilTemperatureInit(Simulation &);

void SoilHeatFlux(State &, double, int, int, int, int, double);

void CanopyBalance(int, int, double, double, double, double, double, double, double, double &, const int &);

double SensibleHeatTransfer(double, double, double, double);

double ThermalCondSoil(double, double, int);
