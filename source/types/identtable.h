/***********************************************************************************
 *
 *   identtable.h -- Identifier Table
 *
 *   Erik Persson and Björn Nilsson, 1999-
 */

#ifndef IDENTTABLE_H
#define IDENTTABLE_H

#include <types/stringref.h>
#include <afxcoll.h>

class CIdentData;

typedef CIdentData *HIDENT;
#define INVALID_IDENT (NULL)

const char *GetIdentString(HIDENT hIdent);

class CIdentTable
{
	CMapStringToPtr m_Map;

public:
	CIdentTable();
	~CIdentTable();

	HIDENT GetIdentHandle(const char *ach);
	const char *GetString(HIDENT hIdent);
};

#endif // IDENTTABLE_H
