#pragma once

void OpenOutputFiles(const std::string &, const std::string &, const int &);

void
DailyOutput(const std::string &, const std::string &, const int &, const int &, const int &, const int &, const int &,
            const int &, const int &, const double &, const double &, const double &, const double &, const ClimateStruct[400]);

void
DataOutput(Simulation &);
