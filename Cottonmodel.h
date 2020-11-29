// CottonModel.h : main header file for the COTTONMODEL application
// "CottonModel" is the simulation module of the Cotton2K model.
// Version 4.0 was written by A. Marani June 2004.
// Compiled by Microsoft Visual C++.Net 2003.
//  This file contains declartions for class C2K.
#pragma once
#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "stdafx.h"
#include "resource.h"       // main symbols
#include "dialogs.h"       // main symbols
#include <string>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;
/////////////////////////////////////////////////////////////////////////////
// C2KApp:
//
class C2KApp : public CWinApp
{
public:
	C2KApp();
    CString GetJobFile();
    void GetProfilesList(fs::path JobFileName);
    void RunTheModel();
	tuple<string> DailySimulation(string ProfileName, const string& Date);
	tuple<string> SimulateThisDay(string ProfileName);
    tuple<BOOL, string> DoAdjustments(string ProfileName, const string& Date);
	virtual BOOL InitInstance();
	virtual int ExitInstance();
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()

    CProgCtrlDlg* pdlg;         // pointer to the progress control dialog
	CStringArray ProfileArray;  // array of the profile names
};
