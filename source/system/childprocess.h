/***********************************************************************************
 *
 *   ChildProcess.h
 *
 *   Erik Persson and Bj�rn Nilsson, 2000-
 */

#ifndef CHILDPROCESS_H
#define CHILDPROCESS_H

class ITextOutput;

bool DoChildProcess(const char *aCmd, int *pProcessRetVal);

#endif // CHILDPROCESS_H
