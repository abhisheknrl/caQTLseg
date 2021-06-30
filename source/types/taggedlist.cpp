/*******************************************************************************
 *
 *   taggedlist.cpp -- Tagged list
 *
 *   Erik Persson and Björn Nilsson, 2001-
 */

#include <stdio.h>
#include "types/string.h"
#include "types/taggedlist.h"
#include "types/stringfun.h"

CTaggedList::CTaggedList()
{
}

CTaggedList::CTaggedList(const CTaggedList& other)
{
	m_vsKey= other.m_vsKey;
	m_vsValue= other.m_vsValue;
}

CTaggedList::~CTaggedList()
{
}

int CTaggedList::FindKey(const char *aProp) const
{
	// Return index of property or -1 for failure
	for (int i= 0; i<m_vsValue.GetSize(); i++)
		if (!strcmp(m_vsKey[i],aProp))
			return i;
	return -1; // not found
}

bool CTaggedList::AppendAssoc(const char *aProp, const char *aVal)
{
	try
	{
		m_vsKey.Add(aProp);
		m_vsValue.Add(aVal);
		return true;
	}
	catch(...) // CMemoryException
	{
		return false;
	}
}

CString CTaggedList::GetProperty(const char *aProp) const
{
	int n= FindKey(aProp);
	if (n>=0)
	{
		ASSERT( n<m_vsValue.GetSize() );
		return m_vsValue[n];
	}
	return "";
}

bool CTaggedList::SetProperty(const char *aProp, const char *aVal)
{
	try
	{
		int n= FindKey(aProp);
		if (n>=0)
		{
			ASSERT( n<m_vsValue.GetSize() );
			m_vsValue.SetAt(n, aVal);
			return true;
		}
		else
			return AppendAssoc(aProp, aVal);
	}
	catch(...) // CMemoryException
	{
		return false;
	}
}

int CTaggedList::GetCount() const
{
	return m_vsValue.GetSize();
}

CString CTaggedList::GetKey(int n) const
{
	ASSERT(n>=0 && n<m_vsKey.GetSize());
	return m_vsKey[n];
}

CString CTaggedList::GetValue(int n) const
{
	ASSERT(n>=0 && n<m_vsValue.GetSize());
	return m_vsValue[n];
}

void CTaggedList::SetFloat(const char *aKey, const float val)
{
	SetProperty(aKey, ::Format("%f", val));
}

bool CTaggedList::GetFloat(const char *aKey, float &val) const
{
	CString s= GetProperty(aKey);
	if (s.IsEmpty())
		return false; // Fail because empty is illegal format
	return sscanf(s, "%f", &val)>0;
}

bool CTaggedList::GetDouble(const char *aKey, double &val) const
{
	CString s= GetProperty(aKey);
	if (s.IsEmpty())
		return false; // Fail because empty is illegal format
	val= atof(s);
	// TODO: error check...
	return true; // sscanf(s, "%f", &val)>0;
}

void CTaggedList::SetInt(const char *aKey, const int val)
{
	SetProperty(aKey, ::Format("%d", val));
}

bool CTaggedList::GetInt(const char *aKey, int &val) const
{
	CString s= GetProperty(aKey);
	if (s.IsEmpty())
		return false; // Fail because empty is illegal format
	return sscanf(s, "%d", &val)>0;
}

void CTaggedList::SetBool(const char *aKey, const bool val)
{
	SetProperty(aKey, ::Format("%d", int(val)));
}

bool CTaggedList::GetBool(const char *aKey, bool &val) const
{
	CString s= GetProperty(aKey);
	if (s.IsEmpty())
		return false; // Fail because empty is illegal format
	int i;
	if (sscanf(s, "%d", &i)==0)
		return false;
	val= i!=0;
	return true;
}

void CTaggedList::SetString(const char *aKey, const CString &val)
{
	(void) SetProperty(aKey, val);
}

bool CTaggedList::GetString(const char *aKey, CString &val) const
{
	val= GetProperty(aKey);
	return true; // Cannot fail as any format is allowed
}
