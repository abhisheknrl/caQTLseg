/***********************************************************************************
 *
 *   TimeHelpers.h -- CTime Helper Functions
 *
 *   Erik Persson and Björn Nilsson, 2000-
 */

#ifndef TIMEHELPERS_H
#define TIMEHELPERS_H

#include <time.h>

inline CString GetTimeString(const CTime& t)
{
	time_t ltime= t.GetTime();
	return CString(ctime( &ltime ));
}

inline CTime GetTime(void)
{
	time_t osBinaryTime;  // C run-time time (defined in <time.h>)
	time( &osBinaryTime ); 
	return CTime(osBinaryTime);
}

#endif // TIMEHELPERS_H
