/***********************************************************************************
 *
 *   ChildProcess.cpp
 *
 *   Erik Persson and Björn Nilsson, 2000-
 */

#include <afx.h>
#include <afxwin.h>
#include "ChildProcess.h"
#include <system/TextOutput.h>
#include <system/FileSystem.h>
#include <system/stdout.h>

/********************************* DoChildProcess *********************************/

#define BUFSIZE 200

#define UNDOCUMENTED
// #define FUSKINPUT
// #define INPUT

#ifdef UNDOCUMENTED
/* msdos.h is not available but here are the flags */
#define FOPEN           0x01    /* file handle open */
#define FEOFLAG         0x02    /* end of file has been encountered */
#if defined (_M_M68K) || defined (_M_MPPC)
#define FWRONLY         0x04    /* file handle associated with write only file */
#define FLOCK           0x08    /* file has been successfully locked at least once */
#else  /* defined (_M_M68K) || defined (_M_MPPC) */
#define FCRLF           0x04    /* CR-LF across read buffer (in text mode) */
#define FPIPE           0x08    /* file handle refers to a pipe */
#endif  /* defined (_M_M68K) || defined (_M_MPPC) */
#ifdef _WIN32
#define FNOINHERIT      0x10    /* file handle opened _O_NOINHERIT */
#else  /* _WIN32 */
#define FRDONLY         0x10    /* file handle associated with read only file */
#endif  /* _WIN32 */

#define FAPPEND         0x20    /* file handle opened O_APPEND */
#define FDEV            0x40    /* file handle refers to device */
#define FTEXT           0x80    /* file handle is in text mode */
#endif

#include <windows.h>
#include <io.h>

/*************************************/

void ReportErr(const char *ach)
{
	TRACE("%s\n",ach);
	//::MessageBox(NULL, ach, "Error", MB_OK | MB_ICONINFORMATION );
}

HANDLE MakeUninheritable(HANDLE h)
{
    HANDLE hDup;
	BOOL fSuccess = DuplicateHandle(GetCurrentProcess(), h,
        GetCurrentProcess(), &hDup , 0,
        FALSE,
        DUPLICATE_SAME_ACCESS);
    if( !fSuccess )
		ReportErr("MakeUninheritable error");
    CloseHandle(h);
	return hDup;
}

#ifdef FUSKINPUT
BOOL CALLBACK onWnd(  HWND hwnd,      // handle to parent window
  LPARAM lParam)   // application-defined value); 
{
    DWORD procId;
    CChildJob *pJob= (CChildJob *) lParam;

    if (GetWindowThreadProcessId(hwnd, &procId) ==  pJob->m_dwChildThreadId)
    {
        pJob->m_hChildWnd= hwnd;
        return FALSE; // Stop enumeration
    }
    return TRUE;
}
#endif

/***********************************************************************************/

class CChildJob
{
public:
	HANDLE m_hChildStdoutRd;
	DWORD m_dwChildThreadId;
	CWinThread *m_pThreadFromChild;
	CString m_sOut;

	bool DoChildProcess(const char *aCmd, int *pProcessRetVal);
	void TransferFromChild();
};

/***********************************************************************************/

void CChildJob::TransferFromChild()
{
/*   BOOL  fConnected = ConnectNamedPipe(m_hChildStdoutRd, NULL) ? 
         TRUE : 
         (GetLastError() == ERROR_PIPE_CONNECTED); 

	if (!fConnected) return;
	*/

	bool bProceed= true;
	CHAR aBuf[BUFSIZE]; 
	DWORD dwBytes;

	while (bProceed)
	{
		if (ReadFile( m_hChildStdoutRd, aBuf, BUFSIZE-1, &dwBytes, NULL))
		{	// && dwBytes
			// Read success
			aBuf[dwBytes]= 0;

			// Put information in buffer
			m_sOut += aBuf;
		}
		else 
		{
			ReportErr("Read error"); // GetLastError();
			bProceed= false;
		}
	}
}

UINT CChildJob_ThreadFromChild( LPVOID pParam )
{
	//Controlling function
	CChildJob *pS= (CChildJob *) pParam;
	pS->TransferFromChild();

	return 0; // means OK
}

/***********************************************************************************/

bool CChildJob::DoChildProcess(const char *aCmd, int *pProcessRetVal)
{ 
	// Error can be retreived by GetLastError()
#ifdef INPUT
	HANDLE hChildStdinRd, hChildStdinWr;
#endif
	HANDLE hChildStdoutRd, hChildStdoutWr;
	DWORD retval;
	bool bOk= false;

	m_sOut= "";

	SECURITY_ATTRIBUTES pipe_attr; 

	// First clear
#ifdef INPUT
	m_hChildStdinWr= INVALID_HANDLE_VALUE;
	m_pThreadToChild= NULL;
#endif
	m_hChildStdoutRd= INVALID_HANDLE_VALUE;
	m_pThreadFromChild= NULL;
	
	// If bInheritHandle is set, pipe handles are inherited,
	// and the pipe won't work properly 
 
	pipe_attr.nLength = sizeof(SECURITY_ATTRIBUTES); 
	pipe_attr.bInheritHandle = TRUE; // FALSE;
	pipe_attr.lpSecurityDescriptor = NULL; 

#ifdef INPUT
#ifdef FUSKINPUT
	// no stdin pipe
	hChildStdinRd= INVALID_HANDLE_VALUE;
	hChildStdinWr= INVALID_HANDLE_VALUE;
#else
	// Create stdin pipe
	if (!CreatePipe(&hChildStdinRd, &hChildStdinWr, &pipe_attr, 0)) 
	{	ReportErr("Stdin pipe creation failed\n"); 
		goto cleanup0;
	}
#endif
#endif
  
	// Create stdout pipe
	if (! CreatePipe(&hChildStdoutRd, &hChildStdoutWr, &pipe_attr, 0)) 
	{	ReportErr("Stdout pipe creation failed\n"); 
		goto cleanup0;
	}
	
	// Save handles
#ifdef INPUT
#ifndef FUSKINPUT
	m_hChildStdinWr= MakeUninheritable(hChildStdinWr);
#endif
#endif
	m_hChildStdoutRd= MakeUninheritable(hChildStdoutRd);

/*	{
		HANDLE hSaveStdout;

		hSaveStdout= GetStdHandle(STD_OUTPUT_HANDLE);
		(void) SetStdHandle( STD_OUTPUT_HANDLE, hChildStdoutWr);

		hOut = CreateFile("CONOUT$", GENERIC_WRITE, FILE_SHARE_WRITE,
				NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
		FILE_FLAG_NO_BUFFERING
	}
*/
	PROCESS_INFORMATION proc_info; 
	STARTUPINFO startup_info; 

	// Set up members of STARTUPINFO structure. 

	ZeroMemory( &startup_info, sizeof(STARTUPINFO) );
	startup_info.cb = sizeof(STARTUPINFO); 
#ifdef UNDOCUMENTED
	startup_info.hStdInput= INVALID_HANDLE_VALUE;
	startup_info.hStdOutput= INVALID_HANDLE_VALUE;
	startup_info.hStdError= INVALID_HANDLE_VALUE;
#else
#ifdef FUSKINPUT
	startup_info.hStdInput= INVALID_HANDLE_VALUE; //hChildStdinRd; 
#else	
	startup_info.hStdInput= hChildStdinRd; 
#endif
	startup_info.hStdOutput= hChildStdoutWr;
	startup_info.hStdError= hChildStdoutWr; // INVALID_HANDLE_VALUE;
#endif
	startup_info.dwFlags= STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW;
	startup_info.wShowWindow= SW_HIDE; // SW_SHOWMINNOACTIVE;

#ifdef UNDOCUMENTED
	// Use undocumented redirection mechanism
	{
		int nh= 3;
        char *aFileFlag;
        UNALIGNED long *ahFile;

        startup_info.cbReserved2 = (WORD)(sizeof( int ) + (nh *
                                  (sizeof( char ) + sizeof( long ))));

        startup_info.lpReserved2 = (unsigned char *) malloc( startup_info.cbReserved2);

        *((UNALIGNED int *)(startup_info.lpReserved2)) = nh;

        aFileFlag = (char *)(startup_info.lpReserved2 + sizeof( int ));

        ahFile = (UNALIGNED long *)(startup_info.lpReserved2 + sizeof( int ) +
                  (nh * sizeof( char )));

		ahFile[0]= (long) INVALID_HANDLE_VALUE;
		aFileFlag[0]= (char) (FPIPE | FOPEN | FTEXT);
#ifdef INPUT
#ifndef FUSKINPUT
		ahFile[0]= (long) INVALID_HANDLE_VALUE;
		aFileFlag[0]= (char) (FPIPE | FOPEN | FTEXT);
#endif
#endif
		ahFile[0]= (long) INVALID_HANDLE_VALUE;
		aFileFlag[0]= (char) (FPIPE | FOPEN | FTEXT);

		// Stdout: use illegal flag combination. (both pipe and device!)
		ahFile[1]= (long) hChildStdoutWr;
		aFileFlag[1]= (char) (FPIPE | FDEV | FOPEN | FTEXT);

		// Stderr: use illegal flag combination. (both pipe and device!)
		ahFile[2]= (long) hChildStdoutWr; // INVALID_HANDLE_VALUE
		aFileFlag[2]= (char) (FPIPE | FDEV | FOPEN | FTEXT);
	}
#endif

	// Create the child process. 
	if (!CreateProcess(NULL, 
		(char *) aCmd,  // pointer to command line string 
		NULL,          // process security attributes 
		NULL,          // primary thread security attributes 
		TRUE,          // handles are inherited 
		CREATE_NEW_CONSOLE, // 0 // creation flags 
		NULL,          // use parent's environment 
		NULL,          // use parent's current directory 
		&startup_info,  // STARTUPINFO pointer 
		&proc_info))  // receives PROCESS_INFORMATION 
	{	ReportErr("Create process failed"); 
		goto cleanup1;
	}

#ifdef UNDOCUMENTED
	// Free allocated structure
    free( startup_info.lpReserved2 );
#endif

#ifdef INPUT
#ifndef FUSKINPUT
    (void) CloseHandle(hChildStdinRd);
#endif
#endif
	(void) CloseHandle(hChildStdoutWr);

	// Save thread ID
	m_dwChildThreadId= proc_info.dwThreadId;

#ifdef INPUT
#ifdef FUSKINPUT
    // Find window
    
    WaitForInputIdle(proc_info.hProcess,10000);

    m_hChildWnd= NULL;
    (void) EnumWindows(onWnd,(LPARAM) this);
#endif
#endif

#ifdef INPUT
	m_pThreadToChild= AfxBeginThread(CChildJob_ThreadToChild, (LPVOID) this);
	//m_pThreadToChild->m_bAutoDelete= false;
#endif
	
	m_pThreadFromChild= AfxBeginThread(CChildJob_ThreadFromChild, (LPVOID) this);
	//m_pThreadFromChild->m_bAutoDelete= false;

	WaitForSingleObject(proc_info.hProcess, (DWORD)(-1L));
	//Sleep(20000);

    /* return termination code and exit code -- note we return
       the full exit code */
	bOk= (GetExitCodeProcess(proc_info.hProcess, &retval) != 0);

	CloseHandle(proc_info.hProcess);

	//if (m_pThreadToChild) delete m_pThreadToChild;
	//if (m_pThreadFromChild) delete m_pThreadFromChild;

cleanup1:
#ifdef INPUT
	if (m_hChildStdinWr != INVALID_HANDLE_VALUE)
		(void) CloseHandle(m_hChildStdinWr);
	m_hChildStdinWr= INVALID_HANDLE_VALUE;
#endif
	if (m_hChildStdoutRd != INVALID_HANDLE_VALUE)
		(void) CloseHandle(m_hChildStdoutRd);
	m_hChildStdoutRd= INVALID_HANDLE_VALUE;

cleanup0:
	*pProcessRetVal= retval;

	printf(m_sOut);
	
	return bOk;
}

bool DoChildProcess(const char *aCmd, int *pProcessRetVal)
{
	// Return false if process couldn't be created
	// Otherwise returns a process return value in *pProcessRetVal

	CChildJob job;

	return job.DoChildProcess(aCmd,pProcessRetVal);
	/*

	STARTUPINFO si = {0};
	si.cb= sizeof(si);

	PROCESS_INFORMATION proc_info;
	
	if (!CreateProcess(NULL,
		  (char *) aCmd, // TODO: varför är dena inte konstant?
		  NULL,	
		  NULL,
		  FALSE,	// bInheritHandles
		  0,		// dwCreationFlags
		  NULL,NULL,
		  &si,&proc_info))
	{	
		// Error can be retreived by GetLastError()
		return false;
	}
	else
	{
		WaitForSingleObject(proc_info.hProcess, (DWORD)(-1L));
	
		DWORD retval= 0;
		(void) GetExitCodeProcess(proc_info.hProcess, &retval);
		*pProcessRetVal= retval;

		CloseHandle(proc_info.hProcess);	
		return true;
	}
	*/
}
