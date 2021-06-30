/*******************************************************************************
 *
 *   IPropertyList.h -- Property List Interface
 *
 *   Erik Persson, 2001
 */

#ifndef IPROPERTYLIST_H
#define IPROPERTYLIST_H

#include "types/string.h"

class IPropertyList
{
	// Abstract Methods
public:
	virtual CString GetProperty(const char *aProp) const = 0;
	virtual bool SetProperty(const char *aProp, const char *aVal)= 0;	
};

#endif
