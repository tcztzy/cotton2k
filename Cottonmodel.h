// CottonModel.h : main header file for the COTTONMODEL application
// "CottonModel" is the simulation module of the Cotton2K model.
// Version 4.0 was written by A. Marani June 2004.
// Compiled by Microsoft Visual C++.Net 2003.
//  This file contains declartions for class C2K.
#pragma once

#include <string>
#include <vector>

/////////////////////////////////////////////////////////////////////////////
// C2KApp:
//
class C2KApp {
   public:
    C2KApp();
    void RunTheModel(const char *);
    void DailySimulation();
    void SimulateThisDay();
    bool DoAdjustments();

    std::vector<std::string> ProfileArray;  // array of the profile names
};
