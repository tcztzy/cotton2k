#pragma once

void OpenOutputFiles(const std::string &, const std::string &, const int &);

void
DailyOutput(Simulation &, uint32_t, const int &, const double &, const double &, const double &, const double &);

void
DataOutput(Simulation &);
