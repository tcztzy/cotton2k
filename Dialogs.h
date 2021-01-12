// Dialogs.h : header file
// declarations for classes InputDatesDlg, COpenDlg, CProgCtrlDlg, m_Progress.
//
#pragma once

#include "resource.h"

/////////////////////////////////////////////////////////////////////////////
// InputDatesDlg dialog
//
class InputDatesDlg : public CDialog {
public:
    InputDatesDlg(CWnd *pParent = NULL);   // standard constructor

    enum {
        IDD = IDD_DIALOG1
    };
    CString m_OutputType;
    int m_StartDate;
    int m_EndDate;
protected:
    virtual void DoDataExchange(CDataExchange *pDX);    // DDX/DDV support
    virtual void OnOK();

    virtual void OnCancel();

DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
//    class  COpenDlg dialog
/////////////////////////////////////////////////////////////////////////////
//
class COpenDlg : public CFileDialog {
DECLARE_DYNAMIC(COpenDlg)

public:
    COpenDlg(BOOL bOpenFileDialog, // TRUE for FileOpen, FALSE for FileSaveAs
             LPCTSTR lpszDefExt = NULL,
             LPCTSTR lpszFileName = NULL,
             DWORD dwFlags = OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
             LPCTSTR lpszFilter = NULL,
             CWnd *pParentWnd = NULL);

    BOOL OnInitDialog();

    enum {
        IDD = ID_FILE_OPEN
    };
protected:
    virtual void DoDataExchange(CDataExchange *pDX);    // DDX/DDV support
DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
class CProgCtrlDlg : public CDialog {
public:
    CProgCtrlDlg(CWnd *pParent = NULL);

    BOOL Create();

    virtual void ProgressStepit();

    enum {
        IDD = IDD_PROGCTRL
    };
    int m_iDelta;
    UINT m_uiStep;
    UINT m_uiRangeFrom;
    UINT m_uiPos;
    UINT m_uiRangeTo;
    BOOL m_bVertical;
    BOOL m_bSmooth;
    CString m_ProfileName;
    CString m_Running;
    CProgressCtrl m_Progress;

    virtual BOOL OnInitDialog();

protected:
    virtual void DoDataExchange(CDataExchange *pDX);

    CWnd *m_pParent;
    int m_nID;
public:
    afx_msg void OnChangeEditRunning(CString Text);

DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
// m_Progress window
//
class m_Progress : public CProgressCtrl {
public:
    m_Progress();

    virtual ~m_Progress();

protected:
DECLARE_MESSAGE_MAP()
};
