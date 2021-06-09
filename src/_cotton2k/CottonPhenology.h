#pragma once
#include <cinttypes>
#include "Simulation.hpp"

void SimulateFruitingSite(Simulation &, uint32_t, int, int, int, int &, const double &);

void AddFruitingNode(State &, int, int, double, double, double, double[61], double);

void AddFruitingBranch(State &, int, double, double, double, double[61], double);

void AddVegetativeBranch(State &, double, double, double, double[61], double);

void CreateFirstSquare(State &, double, double[61]);

void PreFruitingNode(State &, double, double[61]);
