// Dialogs.cpp : implementation file for dialog classes
//               InputDatesDlg, COpenDlg, CProgCtrlDlg 
//
#include "stdafx.h"
#include "cottonmodel.h"
#include "Dialogs.h"

////////////////////////////////////////////////////////////////////////
//   class InputDatesDlg dialog
////////////////////////////////////////////////////////////////////////
//     This dialogs enables the user to define some optional outputs: type of output,
//  and starting and ending date for this output.
//     It is used by functions DayClim() and SoilTemperatureInit().
//
InputDatesDlg::InputDatesDlg(CWnd* pParent /*=NULL*/)
	: CDialog(InputDatesDlg::IDD, pParent)
{
	m_OutputType = _T("");
	m_StartDate = 0;
	m_EndDate = 0;
}

void InputDatesDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, m_OutputType);
	DDV_MaxChars(pDX, m_OutputType, 100);
	DDX_Text(pDX, IDC_EDIT2, m_StartDate);
	DDV_MinMaxInt(pDX, m_StartDate, 1, 730);
	DDX_Text(pDX, IDC_EDIT3, m_EndDate);
	DDV_MinMaxInt(pDX, m_EndDate, 1, 730);
}

BEGIN_MESSAGE_MAP(InputDatesDlg, CDialog)
END_MESSAGE_MAP()

void InputDatesDlg::OnOK() 
{
	CDialog::OnOK();
}

void InputDatesDlg::OnCancel() 
{
	CDialog::OnCancel();
}
/////////////////////////////////////////////////////////////////////////////
//     class   COpenDlg dialog
/////////////////////////////////////////////////////////////////////////////
//  This dialog is used to open a "JOB" file.
//     It is used by function GetJobFile().
//
IMPLEMENT_DYNAMIC(COpenDlg, CFileDialog)
//
COpenDlg::COpenDlg(BOOL bOpenFileDialog, LPCTSTR lpszDefExt, LPCTSTR lpszFileName,
        DWORD dwFlags, LPCTSTR lpszFilter, CWnd* pParentWnd) :
          CFileDialog(bOpenFileDialog, lpszDefExt, lpszFileName, dwFlags,
          lpszFilter, pParentWnd)
{
    m_ofn.lpstrInitialDir = "JOBS";
    m_ofn.Flags = OFN_EXPLORER | OFN_ENABLEHOOK | OFN_HIDEREADONLY  | OFN_NOCHANGEDIR;
//
    m_ofn.lpstrTitle = "Open a Job File for Cotton2K";
    m_ofn.lpstrFilter = "COTTON2K Job Files (*.JOB)\0 *.JOB\0\0";
}

void COpenDlg::DoDataExchange(CDataExchange* pDX)
{
      CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(COpenDlg, CDialog)
END_MESSAGE_MAP()
/////////////////////////////////////////////////////////////////////////////
BOOL COpenDlg::OnInitDialog()
{
    VERIFY(CDialog::OnInitDialog());
    return TRUE;
}
/////////////////////////////////////////////////////////////////////////////
//  class  CProgCtrlDlg dialog
/////////////////////////////////////////////////////////////////////////////
//     This dialog is used to show graphically the progress of the run of the simulation.
//     It is used by function RunTheModel()
CProgCtrlDlg::CProgCtrlDlg(CWnd* pParent /*=NULL*/)
    : CDialog(CProgCtrlDlg::IDD, pParent)
{
	m_iDelta = 0;
	m_uiStep = 1;
	m_uiRangeFrom = 0;
	m_uiPos = 0;
	m_bVertical = FALSE;
	m_bSmooth = FALSE;
	m_ProfileName = _T("");
	m_Running = _T("");

	m_nID = CProgCtrlDlg::IDD;
    m_pParent = pParent;
}
///////////////////////////////////////////////////////////////
void CProgCtrlDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDV_MaxChars(pDX, m_ProfileName, 30);
	DDX_Text(pDX, IDC_EDIT_RUNNING, m_Running);
	DDV_MaxChars(pDX, m_Running, 50);
}
////////////////////////////////////////////////////////////////
BOOL CProgCtrlDlg::Create()
{
	return CDialog::Create(m_nID, m_pParent);
}
////////////////////////////////////////////////////////////////
BOOL CProgCtrlDlg::OnInitDialog()
{
	CDialog::OnInitDialog();
    CString m_Title = "Simulation of  ";
    m_Title += m_ProfileName;
    SetWindowText(m_Title);
//	
	CRect rect(20, 30, 520, 45);
// Initialise controls
	m_Progress.Create( WS_VISIBLE | PBS_SMOOTH, rect, this, IDC_PROG_HORZPOS );
	m_Progress.SetRange( static_cast<short>(m_uiRangeFrom), static_cast<short>(m_uiRangeTo) );
	m_Progress.SetStep( m_uiStep );

    return TRUE;
}

BEGIN_MESSAGE_MAP(CProgCtrlDlg, CDialog)
END_MESSAGE_MAP()
/////////////////////////////////////////////////////////////////////////////
void CProgCtrlDlg::ProgressStepit()
{
	m_Progress.StepIt();
	UpdateData();
}
/////////////////////////////////////////////////////////////////////////////
// m_Progress

m_Progress::m_Progress()
{
}

m_Progress::~m_Progress()
{
}

BEGIN_MESSAGE_MAP(m_Progress, CProgressCtrl)
END_MESSAGE_MAP()

void CProgCtrlDlg::OnChangeEditRunning(CString Text) 
{
    m_Running = Text;
    UpdateData(FALSE);
}
