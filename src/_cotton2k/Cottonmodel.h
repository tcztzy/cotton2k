// CottonModel.h : main header file for the COTTONMODEL application
// "CottonModel" is the simulation module of the Cotton2K model.
// Version 4.0 was written by A. Marani June 2004.
// Compiled by Microsoft Visual C++.Net 2003.
//  This file contains declartions for class C2K.
#pragma once

#include "Simulation.hpp"

/////////////////////////////////////////////////////////////////////////////
// C2KApp:
//
class C2KApp
{
public:
    C2KApp();

    void DailySimulation(Simulation &);

    static void SimulateThisDay(Simulation &, const int &);
};