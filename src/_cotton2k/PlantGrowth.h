#pragma once

#include <string>
#include "Simulation.hpp"

using namespace std;

void DryMatterBalance(State &, double &, double &, double &, double &, double);

void ActualFruitGrowth(State &);

void ActualLeafGrowth(State &);
