#pragma once
#include <cinttypes>
#include "Simulation.hpp"

void ApplyFertilizer(Simulation &, unsigned int);
void RootsCapableOfUptake(Soil &);
double AveragePsi(const State &, double);
void WaterUptake(Simulation &, unsigned int); // UPTAKE
void WaterTable(Simulation &, unsigned int); // WATERTBL
void GravityFlow(SoilCell[40][20], double, double);
// SoilProcedure_2
void CapillaryFlow(Simulation &, unsigned int);
void DripFlow(SoilCell[40][20], double, double);
void SoilSum(State &, double);
