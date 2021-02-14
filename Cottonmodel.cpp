// CottonModel.cpp  "CottonModel" is the simulation module of the Cotton2K model.
// Version 4.0 was written by A. Marani June 2004.
// Compiled by Microsoft Visual C++.Net 2003.
// This file defines the following functions in this file:
//       Class C2KApp:
//    Message map and constructor
//    InitInstance()
//    ExitInstance()
//    GetJobFile()
//    GetProfilesList()
//    RunTheModel()
//    DailySimulation()
//    SimulateThisDay()
//    OnAppAbout()
//       Class CAoutDlg
//
#include "stdafx.h"
#include "Cottonmodel.h"
#include "MainFrm.h"
#include "GeneralFunctions.h"
#include "CottonPhenology.h"
#include "DailyClimate.h"
#include "Input.h"
#include "Output.h"
#include "PlantGrowth.h"
#include "PlantNitrogen.h"
#include "SoilNitrogen.h"
#include "SoilProcedures.h"
#include "SoilTemperature.h"

///////////////////////////////////////////////////////////////////////////////
//    Class C2KApp
///////////////////////////////////////////////////////////////////////////////
BEGIN_MESSAGE_MAP(C2KApp, CWinApp)
ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// C2KApp construction
C2KApp::C2KApp() = default;

C2KApp theApp;

/////////////////////////////////////////////////////////////////////////////
BOOL C2KApp::InitInstance()
//     InitInstance() starts the application. it calls the functions:
//  GetJobFile(), GetProfilesList() and RunTheModel().
//     Global variables set:  FrameTitle
//
{
    AfxEnableControlContainer();
    SetRegistryKey(_T("CottonModel"));
    //
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
    //     Parse command line for standard shell commands.
    CCommandLineInfo cmdInfo;
    ParseCommandLine(cmdInfo);
    CString JobFileName; // Name of the Job file (with full directory).
                         //     The name of the "Job" file is usually taken from the command line.
    JobFileName = cmdInfo.m_strFileName;
    //     If the "Job" file is not defined in the command line, call GetJobFile().
    if (JobFileName.GetLength() == 0)
        JobFileName = GetJobFile();
    //     Build and set title of the main frame.
    string FrameTitle = "";
    CString Slash = "\\";
    int nLen = JobFileName.GetLength();
    //
    for (int i = nLen - 1; i >= 0; i--)
    {
        if (JobFileName[i] == Slash)
            break;
        FrameTitle = JobFileName[i] + FrameTitle;
    }
    //
    CFrameWnd *pFrame = new CMainFrame(FrameTitle);
    m_pMainWnd = pFrame;
    if (!pFrame->LoadFrame(IDR_MAINFRAME))
        return FALSE;
    //     The main window has been initialized, so show and update it.
    pFrame->ShowWindow(m_nCmdShow);
    pFrame->UpdateWindow();
    try
    {
        // Get the list of the "Profiles" to run.
        for (std::string profile : GetProfilesList(fs::path((string)JobFileName)))
        {
            // Now run the model.
            RunTheModel(profile.c_str());
        }
    }
    catch (Cotton2KException e)
    {
        AfxMessageBox(e.what());
    }
    //     End of simulation of all profiles
    AfxMessageBox(" Simulation Ended. \n\n To Exit - close the Job window. ");
    return TRUE;
}

/////////////////////////////////////////////////////////////////////////
CString C2KApp::GetJobFile()
//     Function GetJobFile() is called from InitInstance() when the name of the Job file
//  has not been obtained from parsing the command line. It enables the user to input
//  this name by a standard "open file" dialog.
//     It Checks the file extension and returns the name of the Job file (with
//  full directory path).
//
{
    //     Activate the extended "open file" dialog for windows
    COpenDlg dlg(TRUE);
    if (dlg.DoModal() != IDOK)
        return "No File found";
    //     Get file name with path
    CString strFileName = dlg.GetPathName();
    CString Exten = _T(".JOB");
    //     Check file extension and define Job name
    if (strFileName.Right(4) == Exten)
        return strFileName;
    else
    {
        int strlen = strFileName.GetLength();
        int newlen = strlen;
        for (int ii = 4; ii > 0; ii--)
        {
            char Dummy = strFileName[strlen - ii];
            if (Dummy == '.')
            {
                newlen = strlen - ii;
                return strFileName.Left(newlen) + Exten;
            }
        }
        return strFileName.Left(newlen) + Exten;
    }
}

/////////////////////////////////////////////////////////////////////////
std::vector<std::string> C2KApp::GetProfilesList(const fs::path &JobFileName)
//     Function GetProfilesList() opens the "JOB" file, reads it, and gets
//  from it the profile list. It is called from InitInstance().
//     Input argument: name of the JOB file (with path)
//
{
    //     Check if the file exists
    if (!fs::exists(JobFileName))
        throw FileNotExists(JobFileName);
    //     Open the selected Job file for input.
    ifstream DataFile(JobFileName, ios::in);
    if (DataFile.fail())
        //     When the file can not be opened:
        throw FileNotOpened(JobFileName);
    char m_TempString[100];
    //     Skip the first line of the file.
    DataFile.getline(m_TempString, 99);
    //     Read the simulation profiles from the next lines of the file, and store them
    //  in the StringArray ProfileArray.
    std::vector<std::string> ProfileArray;
    int readlength = 99;
    while (readlength > 0)
    {
        //     Create array of Profile names
        DataFile.get(m_TempString, 21);
        if (DataFile.eof() == 1)
            break;
        std::string m_String = m_TempString;
        readlength = m_String.length();
        if (readlength > 4)
            ProfileArray.push_back(m_String.substr(0, readlength - 4));
        DataFile.getline(m_TempString, 99);
    }
    DataFile.close();
    return ProfileArray;
}

/////////////////////////////////////////////////////////////////////////////
void C2KApp::RunTheModel(const char *profile)
//     This function calls the following functions for each profile:
//          ReadInput(), DailySimulation() and  DataOutput()
//
{
    // Read the input data for this simulation
    Simulation sim = ReadInput(profile);
    // Create a modeless dialog with progress control
    pdlg = new CProgCtrlDlg;
    pdlg->m_uiRangeTo = sim.day_finish - sim.day_start + 1;
    pdlg->m_ProfileName = profile;
    pdlg->m_Running = "Running the Simulation";
    pdlg->Create();
    // Do daily simulations
    DailySimulation(sim);
    //     Write output data
    pdlg->m_Running = "Writing Output Files";
    DataOutput(sim);
    pdlg->EndDialog(0);
    delete pdlg; //  check if needed
}

////////////////////////////////////////////////////////////////////////////////
void C2KApp::DailySimulation(Simulation &sim)
//     This function controls the dynamic phase of the simulation.
//     It calls the functions:
//        SimulateThisDay(), WriteStateVariables()
//
//     The following global variable are referenced:   Kday.
//
{
    try
    {
        for (int i = 0; i < sim.day_finish - sim.day_start + 1; i++)
        {
            pdlg->ProgressStepit();
            if (i > 0)
            {
                memcpy(&sim.states[i], &sim.states[i - 1], sizeof(State));
            }
            else
            {
                State &state0 = sim.states[0];
                state0.number_of_layers_with_root = 7;
                state0.plant_height = 4.0;
                state0.abscised_fruit_sites = 0;
                state0.abscised_leaf_weight = 0;
                state0.cumulative_nitrogen_loss = 0;
                state0.water_stress = 1;
                state0.carbon_stress = 1;
                state0.extra_carbon = 0;
                state0.number_of_vegetative_branches = 1;
                state0.number_of_squares = 0;
                state0.number_of_green_bolls = 0;
                state0.number_of_open_bolls = 0;
                for (int k = 0; k < 3; k++)
                {
                    state0.vegetative_branches[k].number_of_fruiting_branches = 0;
                    for (int l = 0; l < 30; l++)
                    {
                        state0.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes = 0;
                        state0.vegetative_branches[k].fruiting_branches[l].delay_for_new_node = 0;
                        state0.vegetative_branches[k].fruiting_branches[l].main_stem_leaf = {0, 0, 0, 0, 0, 0};
                        for (int m = 0; m < 5; m++)
                            state0.site[k][l][m] = {0, 0, 0, Stage::NotYetFormed, {0, 0, 0, 0}, {0, 0}, {0, 0, 0}, {0, 0}, {0, 0}};
                    }
                }
            }
            strcpy(sim.states[i].date, DoyToDate(sim.day_start + i, sim.year));
            SimulateThisDay(sim, i);
        }
    }
    catch (SimulationEnd)
    {
    }
}

//////////////////////////////////////////////////
void C2KApp::SimulateThisDay(Simulation &sim, const int &u)
//     This function executes all the simulation computations in a day. It is called from
//  DailySimulation().   It calls the following functions:
//     DoyToDate(), ColumnShading(), DayClim(), SoilTemperature(), SoilProcedures(),
//     SoilNitrogen(), SoilSum(), PhysiologicalAge(), Defoliate(), Stress(),
//     GetNetPhotosynthesis(), PlantGrowth(), CottonPhenology(), PlantNitrogen(),
//     CheckDryMatterBal(), DailyOutput();
//
//     The following global variables are referenced here:
//  iyear, Kday, LastDayWeatherData, LeafAreaIndex, OutIndex.
//
//     The following global variables are set here:
//  isw, Kday.
//
{
    double rracol[20]; // the relative radiation received by a soil column, as affected by shading by plant canopy.
    // days from emergence
    if (sim.day_emerge <= 0)
        Kday = 0;
    else
        Kday = sim.day_start - sim.day_emerge + u + 1;
    if (Kday < 0)
        Kday = 0;
    //     The following functions are executed each day (also before emergence).
    ColumnShading(sim, u, rracol);   // computes light interception and soil shading.
    DayClim(sim, u);                 // computes climate variables for today.
    SoilTemperature(sim, u, rracol); // executes all modules of soil and canopy temperature.
    SoilProcedures(sim, u);          // executes all other soil processes.
    SoilNitrogen(sim, u);            // computes nitrogen transformations in the soil.
    SoilSum(sim);                    // computes totals of water and N in the soil.
                                     //     The following is executed each day after plant emergence:
    if (sim.day_start + u >= sim.day_emerge && isw > 0)
    {
        //     If this day is after emergence, assign to isw the value of 2.
        isw = 2;
        sim.states[u].day_inc = PhysiologicalAge(sim.states[u].hours); // physiological days increment for this day. computes physiological age
        Defoliate(sim, u);                                                                           // effects of defoliants applied.
        Stress(sim, u);                                                                              // computes water stress factors.
        GetNetPhotosynthesis(sim, u, sim.states[u].day_length);                                      // computes net photosynthesis.
        PlantGrowth(sim, u, 3, sim.states[u].day_length);                                            // executes all modules of plant growth.
        CottonPhenology(sim, u);                                                                     // executes all modules of plant phenology.
        PlantNitrogen(sim, u);                                                                       // computes plant nitrogen allocation.
        CheckDryMatterBal(sim.profile_name, sim.states[u].date, sim.states[u].abscised_leaf_weight); // checks plant dry matter balance.
                                                                                                     //     If the relevant output flag is not zero, compute soil nitrogen balance and soil
                                                                                                     //  nitrogen averages by layer, and write this information to files.
    }
    //     Call DailyOutput for reporting some simulation data for this day.
    DailyOutput(sim, u);
    //     Check if the date to stop simulation has been reached, or if this is the last day
    //  with available weather data. Simulation will also stop when no leaves remain on the plant.
    //
    if (sim.day_start + u >= LastDayWeatherData || (Kday > 10 && LeafAreaIndex < 0.0002))
        throw SimulationEnd();
}
/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
    CAboutDlg();

    enum
    {
        IDD = IDD_ABOUTBOX
    };

protected:
    virtual void DoDataExchange(CDataExchange *pDX); // DDX/DDV support
    DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange *pDX)
{
    CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
//
void C2KApp::OnAppAbout()
{
    CAboutDlg aboutDlg;
    aboutDlg.DoModal();
}
