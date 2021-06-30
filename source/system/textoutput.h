/*******************************************************************************
 *
 *   TextOutput.h -- Text Output Class
 *
 *   Erik Persson and Björn Nilsson, 2000-
 */

#ifndef TEXTOUTPUT_H
#define TEXTOUTPUT_H

#include <stdarg.h>          // for va_list
#include "types/stringfun.h" // for FormatVA

// Interface for text output 
class ITextOutput
{
public:
	// The single interface function
	virtual void OnPut(const char *aText, size_t TextLen) = 0;
	
	// Inline helper functions
	void Write(const char *ach)
	{	
		OnPut(ach, strlen(ach));
	}

	void Write(char c)
	{	
		OnPut(&c, 1);
	}

	void WriteFmtVA(const char *ach, va_list argptr)
	{
		Write( FormatVA(ach, argptr) );
	}

	void WriteFmt(const char *ach, ...)
	{
		va_list args;
		va_start(args,ach);
		WriteFmtVA(ach,args);
		va_end(args);
	}
};

#endif // TEXTOUTPUT_H
