/***********************************************************************************
 *
 *   stderr.h -- Standard Error using ITextOutput
 *
 */

#ifndef SYSTEM_STDERR_H
#define SYSTEM_STDERR_H

class ITextOutput;
ITextOutput *SetStderr(ITextOutput *pOut);

void eprintf(const char *ach, ...);
void eputchar(char c);
void eputs(const char *ach);

#endif
