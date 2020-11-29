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
#include "CottonSimulation.h"
#include "GeneralFunctions.h"

///////////////////////////////////////////////////////////////////////////////
//    Class C2KApp
///////////////////////////////////////////////////////////////////////////////
BEGIN_MESSAGE_MAP(C2KApp, CWinApp)
	ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
	ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
END_MESSAGE_MAP()
/////////////////////////////////////////////////////////////////////////////
// C2KApp construction
C2KApp::C2KApp()
{
}
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
    _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
//     Parse command line for standard shell commands.
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);
    CString JobFileName;     // Name of the Job file (with full directory).
//     The name of the "Job" file is usually taken from the command line.
    JobFileName = cmdInfo.m_strFileName; 
//     If the "Job" file is not defined in the command line, call GetJobFile().
	if (JobFileName.GetLength() == 0)
        JobFileName = GetJobFile();
//     Build and set title of the main frame.
    FrameTitle = "";
    CString Slash = "\\";
    int nLen = JobFileName.GetLength();
//
    for (int i = nLen-1; i >= 0; i--)
    {
        if (JobFileName[i] == Slash) 
            break;
        FrameTitle = JobFileName[i] + FrameTitle;
    }
//
    CFrameWnd* pFrame = new CMainFrame;
	m_pMainWnd = pFrame;
	if (!pFrame->LoadFrame(IDR_MAINFRAME))
		return FALSE;
//     The main window has been initialized, so show and update it.
	pFrame->ShowWindow(m_nCmdShow);
	pFrame->UpdateWindow();
//     Get the list of the "Profiles" to run.
	GetProfilesList(JobFileName);
//     Now run the model.
    try
    {
        RunTheModel();
    }
    catch (Cotton2KException e)
    {
        AfxMessageBox(e.what());
    }
	return TRUE;
}
/////////////////////////////////////////////////////////////////////////////
int C2KApp::ExitInstance() 
{
	return CWinApp::ExitInstance();
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
        for ( int ii = 4; ii > 0 ; ii-- )
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
void C2KApp::GetProfilesList(CString JobFileName)
//     Function GetProfilesList() opens the "JOB" file, reads it, and gets
//  from it the profile list. It is called from InitInstance().
//     Input argument: name of the JOB file (with path)
//    
{
//     Check if the file exists
    CFile file;
    CFileStatus status;
    if (!file.GetStatus(JobFileName, status))
    {
        CString strMessage;
        AfxFormatString1(strMessage, IDS_FILE_NOT_EXISTS, JobFileName);
        AfxMessageBox(strMessage);
		return;
    }
//     Open the selected Job file for input.
    ifstream DataFile(JobFileName, ios::in);
    if ( DataFile.fail() )
    {
//     When the file can not be opened:
       DataFile.close();
       CString strMessage;
       AfxFormatString1(strMessage, IDS_FILE_NOT_OPENED, JobFileName);
       AfxMessageBox(strMessage);
       return;
	}
	char m_TempString[100];
//     Skip the first line of the file.
    DataFile.getline(m_TempString, 99);
//     Read the simulation profiles from the next lines of the file, and store them
//  in the StringArray ProfileArray.
    ProfileArray.RemoveAll();
	int readlength = 99;
	while (readlength > 0)
    {
//     Create array of Profile names
           DataFile.get(m_TempString,21);
           if (DataFile.eof() == 1) 
               break;
           CString m_String = (CString) m_TempString;
           m_String.Remove(' ');
	       readlength = m_String.GetLength();
		   if (readlength > 4)
		       ProfileArray.Add( m_String.Left(readlength - 4) );
           DataFile.getline(m_TempString, 99);
	}
    DataFile.close();
}
/////////////////////////////////////////////////////////////////////////////
void C2KApp::RunTheModel()
//     This function calls the following functions for each profile: 
//          ReadInput(), DailySimulation() and  DataOutput()
//     Global variables set: ProfileName 
//     Global variables referenced:  DayFinish, DayStart 
//
{
    string ProfileName; // name of input file with profile data (without the extension ".PRO").
    string Date;        // date string formatted as "dd-MMM-yyyy", for example 25-JUN-2003
    for (int i = 0; i < ProfileArray.GetSize(); i++)
    {
		ProfileName = ProfileArray.GetAt(i);
//     Read the input data for this simulation
        ReadInput(ProfileName);
//     Create a modeless dialog with progress control
        int range = DayFinish - DayStart + 1;
        pdlg = new CProgCtrlDlg;
        pdlg->m_uiRangeTo = range;
        pdlg->m_ProfileName = ProfileName.c_str();
        pdlg->m_Running = "Running the Simulation";
        pdlg->Create();
//     Do daily simulations
        tie(Date) = DailySimulation(ProfileName, Date);
//     Write output data
        pdlg->m_Running = "Writing Output Files";
        tie(Date) = DataOutput(ProfileName, Date);
        pdlg->EndDialog(i);
        delete pdlg;  //  check if needed
    }
//     End of simulation of all profiles
    AfxMessageBox(" Simulation Ended. \n\n To Exit - close the Job window. " );
}
////////////////////////////////////////////////////////////////////////////////
tuple<string> C2KApp::DailySimulation(string ProfileName, const string& Date)
//     This function controls the dynamic phase of the simulation, allowing
//  for in-run adjustments when there is an input of plant map adjustments.
//     It calls the functions:
//        DoAdjustments(), SimulateThisDay(), WriteStateVariables()
//
//     The following global variable are referenced:   DayStart, Kday.
//     The following global variable are set:   bEnd, Daynum.
//
{
      string date = Date;
      Daynum = DayStart - 1;
	  bEnd = FALSE;
//     Start the daily loop. If variable bEnd has been assigned a value
//  of TRUE end simulation.
	  while (! bEnd)
	  {
         BOOL bAdjustToDo;
         tie(bAdjustToDo, date) = DoAdjustments(ProfileName, date);
         pdlg->ProgressStepit();
//     Execute simulation for this day.
         tie(date) = SimulateThisDay(ProfileName);
//     If there are pending plant adjustments, call WriteStateVariables() to write
//  state variables of this day in a scratch file.
         if (bAdjustToDo)
		        WriteStateVariables(TRUE, date);  
	  } // end while
      return make_tuple(date);
}
///////////////////////////////////////////////////////////////////////////////
tuple<BOOL, string> C2KApp::DoAdjustments(string ProfileName, const string& Date)
//     This function is called from DailySimulation(). It checks if plant adjustment data
//  are available for this day and calls the necessary functions to compute adjustment.
//  It calls PlantAdjustments(), SimulateThisDay(), WriteStateVariables()
//
//     The following global variable are referenced:  DayEmerge, Daynum, Kday, 
//     The following global variable are set:   MapDataDate, nadj, NumAdjustDays.
//
{
//     Check if plant map data are available for this day. If there are no more adjustments, return.
	  static int kprevadj = 0;  // day after emergence of the previous plant map adjustment.
      int sumsad = 0;  // sum for checking if any adjustments are due 
      string date = Date;
      for (int i = 0; i < 30; i++)
          sumsad += MapDataDate[i];
      if (sumsad <= 0)
          return make_tuple(FALSE, date);
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
                tie(date) = PlantAdjustments(i, jj, ProfileName, date);
//     If adjustment is necessary, rerun the simulation for the previous NumAdjustDays (number
//  of days) and call WriteStateVariables() to write state variables in scratch structure.
                if (nadj[jj])
                   for(int j1 = 0; j1 < NumAdjustDays; j1++)
                   {
                       tie(date) = SimulateThisDay(ProfileName);
                       if (Kday > 0) 
                           WriteStateVariables(TRUE, date); 
                   }     // end for j1, and if nadj 
            }            // end for jj
//     After finishing this adjustment date, set kprevadj (date of previous adjustment, to 
//  be used for next adjustment), and assign zero to the present msadte, and to array nadj[].
            kprevadj = MapDataDate[i] - DayEmerge + 1;
            MapDataDate[i] = 0;
            for(int jj = 0; jj < 5; jj++)
                 nadj[jj] = FALSE;
            continue;
         }// end if Daynum
      }// end do i
      return make_tuple(TRUE, date);
}
//////////////////////////////////////////////////
tuple<string> C2KApp::SimulateThisDay(string ProfileName)
//     This function executes all the simulation computations in a day. It is called from
//  DailySimulation(), and DoAdjustments().   It calls the following functions:
//     DoyToDate(), ColumnShading(), DayClim(), SoilTemperature(), SoilProcedures(),
//     SoilNitrogen(), SoilSum(), PhysiologicalAge(), Pix(), Defoliate(), Stress(),
//     GetNetPhotosynthesis(), PlantGrowth(), CottonPhenology(), PlantNitrogen(),
//     CheckDryMatterBal(), PlantNitrogenBal(), SoilNitrogenBal(), SoilNitrogenAverage(),
//     DailyOutput();
//     
//     The following global variables are referenced here:  DayEmerge, DayFinish, 
//  DayStart, iyear, Kday, LastDayWeatherData, LeafAreaIndex, OutIndex, pixday.
//
//     The following global variables are set here:
//  bEnd, Date, DayInc, Daynum, DayOfSimulation, isw, Kday.
//
{
//    Compute Daynum (day of year), Date, and DayOfSimulation (days from start of simulation).
         Daynum++;
	     string Date = DoyToDate(Daynum, iyear);
         DayOfSimulation = Daynum - DayStart + 1; 
//    Compute Kday (days from emergence).
      if (DayEmerge <= 0)
         Kday = 0;
      else
         Kday = Daynum - DayEmerge + 1;
      if (Kday < 0) 
          Kday = 0;
//     The following functions are executed each day (also before emergence).
      ColumnShading();      // computes light interception and soil shading.
      DayClim(ProfileName, Date);            // computes climate variables for today.
      SoilTemperature(ProfileName);    // executes all modules of soil and canopy temperature.
      SoilProcedures(ProfileName);     // executes all other soil processes.
      SoilNitrogen();       // computes nitrogen transformations in the soil.
      SoilSum();            // computes totals of water and N in the soil.
//     The following is executed each day after plant emergence:
      if (Daynum >= DayEmerge && isw > 0)
	  {
//     If this day is after emergence, assign to isw the value of 2.
         isw = 2;
         DayInc = PhysiologicalAge();    // computes physiological age
         if(pixday[0] > 0)
             Pix();        // effects of pix applied.
         Defoliate(ProfileName, Date);         // effects of defoliants applied.
         Stress(ProfileName);            // computes water stress factors.
         GetNetPhotosynthesis();         // computes net photosynthesis.
         PlantGrowth(ProfileName, Date);       // executes all modules of plant growth.
         CottonPhenology();              // executes all modules of plant phenology.
         PlantNitrogen(ProfileName);     // computes plant nitrogen allocation.
         CheckDryMatterBal(ProfileName, Date); // checks plant dry matter balance.
//     If the relevant output flag is not zero, compute soil nitrogen balance and soil
//  nitrogen averages by layer, and write this information to files.
         if ( OutIndex[20] > 0 )
		 {
	        PlantNitrogenBal(ProfileName);          // checks plant nitrogen balance.
	        SoilNitrogenBal(ProfileName);           // checks soil nitrogen balance.
	        SoilNitrogenAverage(ProfileName);       // computes average soil nitrogen by layers.
		 }
      }
//     Call DailyOutput for reporting some simulation data for this day.
      DailyOutput(ProfileName, Date);
//     Check if the date to stop simulation has been reached, or if this is the last day
//  with available weather data. Simulation will also stop when no leaves remain on the plant.
//  bEnd = TRUE  indicates stopping this simulation.
//
      if ( Daynum >= DayFinish || Daynum >= LastDayWeatherData 
           || (Kday > 10 && LeafAreaIndex < 0.0002) )
		  bEnd = TRUE;
      return make_tuple(Date);
}
/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

	enum { IDD = IDD_ABOUTBOX };
protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
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
