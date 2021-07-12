#pragma once
#include <cinttypes>
#include "Simulation.hpp"

void RootSummation(State &state, int NumRootAgeGroups, double row_space, double per_plant_area);

void RootImpedance(SoilCell[40][20]);

double RootCultivation(SoilCell[40][20], int, double, double, double);

double RootDeath(SoilCell &, int, int, double);

void RootAging(SoilCell &, int, int);

void LateralRootGrowthLeft(State &, int, int, unsigned int, double);

void LateralRootGrowthRight(State &, int, int, unsigned int, double);

void InitiateLateralRoots();
