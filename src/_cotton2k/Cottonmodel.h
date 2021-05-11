// CottonModel.h : main header file for the COTTONMODEL application
// "CottonModel" is the simulation module of the Cotton2K model.
// Version 4.0 was written by A. Marani June 2004.
// Compiled by Microsoft Visual C++.Net 2003.
//  This file contains declartions for class C2K.
#pragma once

#include <cstdint>
#include "Simulation.hpp"

void DailySimulation(Simulation &);

void SimulateThisDay(Simulation &, uint32_t);

void initialize_state0(State &, uint32_t);
