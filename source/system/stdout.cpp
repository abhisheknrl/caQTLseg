/***********************************************************************************
 *
 *   stdout.cpp -- Standard Output using ITextOutput
 *
 */

#include <stdio.h>
// #include "stdout.h"  
// NOTE: Our stdout.h header must not be included here. If it 
// is, the OnPut code will not choose the standard printf etc.
#include "textoutput.h"

#ifdef SYSTEM_STDOUT_H
#error stdout.h has been included
#endif

/***********************************************************************************/
// CMyNilOutput: default output object

class CMyNilOut: public ITextOutput
{
	void OnPut(const char *ach, size_t len) 
	{ 
		// This must choose standard printf
		printf("%s", (const char *)CString(ach, len)); 
	}
};

static CMyNilOut g_NilOut;

/***********************************************************************************/

static ITextOutput *g_pStdout= &g_NilOut;

ITextOutput *SetStdout(ITextOutput *pOut)
{
	// Set the standard output and return the previous setting
	ASSERT(pOut);
	ITextOutput *pOld= g_pStdout;
	g_pStdout= pOut;
	return pOld;
}

void system_printf(const char *ach, ...)
{
	va_list args;
	va_start(args,ach);
	g_pStdout->WriteFmtVA(ach,args);
	va_end(args);
}

void system_putchar(int c)
{
	g_pStdout->Write(c);
}

void system_puts(const char *ach)
{
	// puts, differently from fputs, automatically appends newline
	g_pStdout->Write(ach);
	g_pStdout->Write('\n');
}
