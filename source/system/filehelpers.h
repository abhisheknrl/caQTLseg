//
// filehelpers.h -- Facilitates working with file and reporting errors
//

#ifndef FILEHELPERS_H
#define FILEHELPERS_H

#include <stdio.h>

bool CheckLicence(const int nDays);
void ReportError(const char *aMsg, const char *aContext, int nVerbosity);
void ReportError(const char *aMsg, const char *aContext);
void ReportWarning(const char *aMsg, const char *aContext, int g_Verbosity);
void ReportWarning(const char *aMsg, const char *aContext);
FILE *OpenInputFile(const char *aFilename);
FILE *OpenInputFile(const char *aFilename, int nVerbosity);
FILE *OpenOutputFile(const char *aFilename);
FILE *OpenOutputFile(const char *aFilename, int nVerbosity);
FILE *OpenOutputFile(const char *aFilename, int nVerbosity, const char *aMessage);
bool CloseFile(FILE *pf);

int GetVerbosity();
int SetVerbosity(int v);

bool IsLittleEndian();
bool IsBigEndian();

#endif
