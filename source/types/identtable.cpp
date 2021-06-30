/***********************************************************************************
 *
 *   identtable.cpp -- Identifier Table
 *
 *   Erik Persson and Björn Nilsson, 1999-
 */

#include <afx.h>
#include "identtable.h"

/***********************************************************************************/

class CIdentData
{
public:
	CString m_s;

	CIdentData(const char *ach)
	{
		m_s= ach;
	}
};

/***********************************************************************************/

const char *GetIdentString(HIDENT hIdent)
{
	ASSERT(hIdent != INVALID_IDENT);
	return hIdent->m_s;
}

/***********************************************************************************/

CIdentTable::CIdentTable()
{
}

CIdentTable::~CIdentTable()
{
	void *p;
	POSITION pos= m_Map.GetStartPosition();
	CString key;
	while (pos)
	{
		m_Map.GetNextAssoc(pos,key,p);

		delete (CIdentData *) p;
	}
}

HIDENT CIdentTable::GetIdentHandle(const char *ach)
{
	void *p;

	if (m_Map.Lookup(ach,p))
		return ((CIdentData *) p);

	CIdentData *pData= new CIdentData(ach);

	m_Map.SetAt(ach,(void *) pData);

	return (HIDENT) pData;
}

const char *CIdentTable::GetString(HIDENT hIdent)
{
	return hIdent->m_s;
}
