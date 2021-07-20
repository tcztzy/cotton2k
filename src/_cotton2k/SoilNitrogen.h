#pragma once
#include "Simulation.hpp"

void UreaHydrolysis(SoilCell &, int, int);
void MineralizeNitrogen(SoilCell &, int, int, const int &, const int &, double);
void Nitrification(SoilCell &, int, int, double);
void Denitrification(SoilCell &, int, int, double);
