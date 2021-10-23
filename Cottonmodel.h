// CottonModel.h : main header file for the COTTONMODEL application
// "CottonModel" is the simulation module of the Cotton2K model.
// Version 4.0 was written by A. Marani June 2004.
// Compiled by Microsoft Visual C++.Net 2003.
//  This file contains declartions for class C2K.
#pragma once
#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include <string>
#include "stdafx.h"
/////////////////////////////////////////////////////////////////////////////
// C2KApp:
//
class C2KApp
{
public:
	C2KApp();
    void GetProfilesList(CString JobFileName);
    void RunTheModel();
	void DailySimulation();
	void SimulateThisDay();
    bool DoAdjustments();
	virtual bool InitInstance(std::string);

    CStringArray ProfileArray;  // array of the profile names
};
