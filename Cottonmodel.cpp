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
//    DoAdjustments()
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
#include "PlantAdjustment.h"
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
        m_String;
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
    pdlg->m_uiRangeTo = sim.number_of_states;
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
//     This function controls the dynamic phase of the simulation, allowing
//  for in-run adjustments when there is an input of plant map adjustments.
//     It calls the functions:
//        DoAdjustments(), SimulateThisDay(), WriteStateVariables()
//
//     The following global variable are referenced:   Kday.
//
{
    string Date;

    int Daynum = sim.day_start - 1;
    int NumLayersWithRoots = 7; // number of soil layers with roots. Initialized with 7
    double PlantHeight = 4.0;
    double AbscisedFruitSites = 0; // total number of abscised fruit sites, per plant.
    double AbscisedLeafWeight = 0; // weight of abscissed leaves, g per plant.
    double WaterStress = 1;        // general water stress index (0 to 1).
    try
    {
        for (int DayOfSimulation = 1; DayOfSimulation <= sim.number_of_states; DayOfSimulation++)
        {
            BOOL bAdjustToDo;
            tie(bAdjustToDo, Date, Daynum, NumLayersWithRoots,
                PlantHeight, AbscisedFruitSites, AbscisedLeafWeight, WaterStress) = DoAdjustments(sim, Date,
                                                                                                  Daynum,
                                                                                                  DayOfSimulation,
                                                                                                  NumLayersWithRoots,
                                                                                                  PlantHeight,
                                                                                                  AbscisedFruitSites,
                                                                                                  AbscisedLeafWeight,
                                                                                                  WaterStress);
            pdlg->ProgressStepit();
            Daynum++;
            tie(Date, NumLayersWithRoots, PlantHeight,
                AbscisedFruitSites, AbscisedLeafWeight, WaterStress) = SimulateThisDay(sim, Daynum,
                                                                                       NumLayersWithRoots, PlantHeight,
                                                                                       AbscisedFruitSites,
                                                                                       AbscisedLeafWeight, WaterStress);
            // If there are pending plant adjustments, call WriteStateVariables() to write
            // state variables of this day in a scratch file.
            if (bAdjustToDo)
                WriteStateVariables(TRUE, Date, Daynum, DayOfSimulation, sim.first_bloom, sim.first_square, NumLayersWithRoots,
                                    PlantHeight, AbscisedFruitSites, AbscisedLeafWeight, WaterStress, sim.root_weight,
                                    sim.root_age, sim.climate);
        }
    }
    catch (SimulationEnd)
    {
    }
}

///////////////////////////////////////////////////////////////////////////////
tuple<BOOL, string, int, int, double, double, double, double> C2KApp::DoAdjustments(
    Simulation &sim,
    string Date,
    int Daynum,
    int DayOfSimulation,
    int NumLayersWithRoots,
    double PlantHeight,
    double AbscisedFruitSites,
    double AbscisedLeafWeight,
    double WaterStress)
//     This function is called from DailySimulation(). It checks if plant adjustment data
//  are available for this day and calls the necessary functions to compute adjustment.
//  It calls PlantAdjustments(), SimulateThisDay(), WriteStateVariables()
//
//     The following global variable are referenced:  Kday,
//     The following global variable are set:   MapDataDate, nadj, NumAdjustDays.
//
{
    //     Check if plant map data are available for this day. If there are no more adjustments, return.
    static int kprevadj = 0; // day after emergence of the previous plant map adjustment.
    int sumsad = 0;          // sum for checking if any adjustments are due
    for (int i = 0; i < 30; i++)
        sumsad += MapDataDate[i];
    if (sumsad <= 0)
        return make_tuple(FALSE, Date, Daynum, NumLayersWithRoots,
                          PlantHeight, AbscisedFruitSites, AbscisedLeafWeight, WaterStress);
    //     Loop for all adjustment data, and check if there is an adjustment for this day.
    for (int i = 0; i < 30; i++)
    {
        if (Daynum == MapDataDate[i])
        {
            //     Compute NumAdjustDays, the number of days for retroactive adjustment. This can not be more
            //  than 12 days, limited by the date of the previous adjustment.
            NumAdjustDays = Kday - kprevadj;
            if (NumAdjustDays > 12)
                NumAdjustDays = 12;
            //     Loop for six possible adjustments. On each iteration call first PlantAdjustments(), which
            //  will assign TRUE to nadj(jj) if adjustment is necessary, and compute the necessary parameters.
            for (int jj = 0; jj < 5; jj++)
            {
                tie(Date, Daynum, sim.first_bloom, sim.first_square, NumLayersWithRoots, PlantHeight, AbscisedFruitSites,
                    AbscisedLeafWeight, WaterStress) = PlantAdjustments(i, jj, sim.profile_name, Date, Daynum,
                                                                        DayOfSimulation, sim.first_bloom, sim.first_square,
                                                                        NumLayersWithRoots, PlantHeight,
                                                                        AbscisedFruitSites, AbscisedLeafWeight,
                                                                        WaterStress, sim.root_weight, sim.root_age);
                //     If adjustment is necessary, rerun the simulation for the previous NumAdjustDays (number
                //  of days) and call WriteStateVariables() to write state variables in scratch structure.
                if (nadj[jj])
                    for (int j1 = 0; j1 < NumAdjustDays; j1++)
                    {
                        Daynum++;
                        tie(Date, NumLayersWithRoots, PlantHeight,
                            AbscisedFruitSites, AbscisedLeafWeight, WaterStress) = SimulateThisDay(sim, Daynum,
                                                                                                   NumLayersWithRoots,
                                                                                                   PlantHeight,
                                                                                                   AbscisedFruitSites,
                                                                                                   AbscisedLeafWeight,
                                                                                                   WaterStress);
                        if (Kday > 0)
                            WriteStateVariables(TRUE, Date, Daynum, DayOfSimulation, sim.first_bloom, sim.first_square,
                                                NumLayersWithRoots, PlantHeight, AbscisedFruitSites, AbscisedLeafWeight,
                                                WaterStress, sim.root_weight, sim.root_age, sim.climate);
                    } // end for j1, and if nadj
            }         // end for jj
                      //     After finishing this adjustment date, set kprevadj (date of previous adjustment, to
                      //  be used for next adjustment), and assign zero to the present msadte, and to array nadj[].
            kprevadj = MapDataDate[i] - sim.day_emerge + 1;
            MapDataDate[i] = 0;
            for (int jj = 0; jj < 5; jj++)
                nadj[jj] = FALSE;
        } // end if Daynum
    }     // end do i
    return make_tuple(TRUE, Date, Daynum, NumLayersWithRoots,
                      PlantHeight, AbscisedFruitSites, AbscisedLeafWeight, WaterStress);
}

//////////////////////////////////////////////////
tuple<string, int, double, double, double, double> C2KApp::SimulateThisDay(
    Simulation &sim,
    const int &Daynum,
    int NumLayersWithRoots,
    double PlantHeight,
    double AbscisedFruitSites,
    double AbscisedLeafWeight,
    double WaterStress)
//     This function executes all the simulation computations in a day. It is called from
//  DailySimulation(), and DoAdjustments().   It calls the following functions:
//     DoyToDate(), ColumnShading(), DayClim(), SoilTemperature(), SoilProcedures(),
//     SoilNitrogen(), SoilSum(), PhysiologicalAge(), Pix(), Defoliate(), Stress(),
//     GetNetPhotosynthesis(), PlantGrowth(), CottonPhenology(), PlantNitrogen(),
//     CheckDryMatterBal(), PlantNitrogenBal(), SoilNitrogenBal(), SoilNitrogenAverage(),
//     DailyOutput();
//
//     The following global variables are referenced here:
//  iyear, Kday, LastDayWeatherData, LeafAreaIndex, OutIndex, pixday.
//
//     The following global variables are set here:
//  DayInc, DayOfSimulation, isw, Kday.
//
{
    //    Compute Daynum (day of year), Date, and DayOfSimulation (days from start of simulation).
    string Date = DoyToDate(Daynum, iyear);
    double DayLength; // day length, hours.
    int DayOfSimulation = Daynum - sim.day_start + 1;
    double rracol[20]; // the relative radiation received by a soil column, as affected by shading by plant canopy.
                       //    Compute Kday (days from emergence).
    if (sim.day_emerge <= 0)
        Kday = 0;
    else
        Kday = Daynum - sim.day_emerge + 1;
    if (Kday < 0)
        Kday = 0;
    //     The following functions are executed each day (also before emergence).
    ColumnShading(sim, Daynum, sim.day_emerge, PlantHeight, rracol); // computes light interception and soil shading.
    tie(DayLength) = DayClim(sim.profile_name, Date, Daynum, DayOfSimulation, sim.day_start, sim.day_finish, sim.latitude,
                             sim.longitude, sim.climate); // computes climate variables for today.
    tie(sim.day_emerge) = SoilTemperature(sim, sim.profile_name, Daynum, DayOfSimulation, sim.day_emerge, sim.day_start, sim.day_finish, sim.day_plant,
                                          sim.day_start_mulch, sim.day_end_mulch, sim.mulch_indicator, sim.mulch_transmissivity_short_wave, sim.mulch_transmissivity_long_wave, PlantHeight,
                                          rracol, sim.climate); // executes all modules of soil and canopy temperature.
    SoilProcedures(sim, Daynum, DayOfSimulation, NumLayersWithRoots, WaterStress); // executes all other soil processes.
    SoilNitrogen(Daynum, sim.day_start);                                           // computes nitrogen transformations in the soil.
    SoilSum();                                                                     // computes totals of water and N in the soil.
                                                                                   //     The following is executed each day after plant emergence:
    if (Daynum >= sim.day_emerge && isw > 0)
    {
        //     If this day is after emergence, assign to isw the value of 2.
        isw = 2;
        double DayInc = PhysiologicalAge(AirTemp); // physiological days increment for this day. computes physiological age
        if (pixday[0] > 0)
            Pix();                                                 // effects of pix applied.
        Defoliate(sim.profile_name, Date, Daynum, sim.day_emerge); // effects of defoliants applied.
        tie(WaterStress) = Stress(sim.profile_name, PlantHeight,
                                  NumLayersWithRoots); // computes water stress factors.
        GetNetPhotosynthesis(Daynum, sim.day_emerge, sim.day_start_co2, sim.day_end_co2, sim.co2_enrichment_factor,
                             DayLength, sim.climate); // computes net photosynthesis.
        tie(NumLayersWithRoots, PlantHeight) = PlantGrowth(
            sim,
            Date,
            Daynum,
            DayOfSimulation,
            3, // the number of root classes defined in the model.
            NumLayersWithRoots,
            PlantHeight,
            DayInc,
            DayLength,
            WaterStress); // executes all modules of plant growth.
        tie(sim.first_bloom, sim.first_square, AbscisedFruitSites, AbscisedLeafWeight) = CottonPhenology(Daynum, sim.day_emerge,
                                                                                                         sim.first_bloom, sim.first_square,
                                                                                                         DayInc, WaterStress,
                                                                                                         AbscisedLeafWeight, sim.climate); // executes all modules of plant phenology.
        PlantNitrogen(sim.profile_name, Daynum, sim.day_emerge);                                                                           // computes plant nitrogen allocation.
        CheckDryMatterBal(sim.profile_name, Date, AbscisedLeafWeight);                                                                     // checks plant dry matter balance.
                                                                                                                                           //     If the relevant output flag is not zero, compute soil nitrogen balance and soil
                                                                                                                                           //  nitrogen averages by layer, and write this information to files.
        if (OutIndex[20] > 0)
        {
            PlantNitrogenBal(sim.profile_name);    // checks plant nitrogen balance.
            SoilNitrogenBal(sim.profile_name);     // checks soil nitrogen balance.
            SoilNitrogenAverage(sim.profile_name); // computes average soil nitrogen by layers.
        }
    }
    //     Call DailyOutput for reporting some simulation data for this day.
    DailyOutput(sim.profile_name, Date, Daynum, DayOfSimulation, sim.day_emerge, sim.day_finish, sim.first_bloom, sim.first_square,
                NumLayersWithRoots, PlantHeight, AbscisedFruitSites, AbscisedLeafWeight, WaterStress, sim.root_weight,
                sim.root_age, sim.climate);
    //     Check if the date to stop simulation has been reached, or if this is the last day
    //  with available weather data. Simulation will also stop when no leaves remain on the plant.
    //
    if (Daynum >= sim.day_finish || Daynum >= LastDayWeatherData || (Kday > 10 && LeafAreaIndex < 0.0002))
        throw SimulationEnd();
    return make_tuple(Date, NumLayersWithRoots, PlantHeight,
                      AbscisedFruitSites, AbscisedLeafWeight, WaterStress);
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
