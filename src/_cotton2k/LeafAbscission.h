#pragma once

#include <cstdint>
#include "Simulation.hpp"

void PreFruitLeafAbscission(State &, double, unsigned int, unsigned int, unsigned int, double);

void MainStemLeafAbscission(State &, int, int, double, unsigned int, unsigned int);

void DefoliationLeafAbscission(State &, unsigned int);
