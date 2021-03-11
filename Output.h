#pragma once
#include "Simulation.hpp"

void OpenOutputFiles(const char *, const char *, int, int);

void DailyOutput(Simulation &, uint32_t);

void DataOutput(Simulation &);
