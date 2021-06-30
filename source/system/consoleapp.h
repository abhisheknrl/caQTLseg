/*******************************************************************************
 *
 *   consoleapp.h -- 
 *
 *   Björn Nilsson, 2004-2008
 */

#ifndef CONSOLEAPP_H
#define CONSOLEAPP_H

#include "types/vector.h"
#include "types/string.h"
#include "assert.h"
#include "filehelpers.h"

class CConsoleApp
{
protected:
	CVector<CString> m_vNonSwitchArgs;

	virtual bool OnSwitch(const char *ach, CString &sSwitch) { return true; };
	virtual const char *ScanSwitch(const char *arg, CString &sSwitch);
	virtual bool ParseCommandLine(int argc, TCHAR* argv[], TCHAR* envp[]);
	virtual ~CConsoleApp() {};

public:
	virtual int Main(int argc, TCHAR* argv[], TCHAR* envp[])=0;
};

#endif
