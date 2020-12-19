#pragma once
#include <string>
#include <tuple>
std::tuple<double> DayClim(const std::string&, const std::string&, const int&, const int&, const int&, const int&, const double&);
double VaporPressure(double);
double clearskyemiss(double, double);
