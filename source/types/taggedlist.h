/*******************************************************************************
 *
 *   taggedlist.h -- Tagged list
 *
 *   Erik Persson and Björn Nilsson, 2001-
 */

#ifndef TAGGEDLIST_H
#define TAGGEDLIST_H

#include "ipropertylist.h"
#include "types/vector.h"
/*
#include "types/ref.h"
#include "types/claimable.h"
*/

class CTaggedList : public IPropertyList
{
	CVector<CString> m_vsKey;
	CVector<CString> m_vsValue;

public:
	CTaggedList();
	CTaggedList(const CTaggedList& other);
	virtual ~CTaggedList();

	// IPropertyList Interface
public:
	CString GetProperty(const char *aProp) const;
	bool SetProperty(const char *aProp, const char *aVal);	

	// CTaggedList methods
public:
	void Clear()
	{
		m_vsKey.SetSize(0);
		m_vsValue.SetSize(0);
	}

	int GetCount() const;
	int FindKey(const char *aKey) const;
	CString GetKey(int n) const;
	CString GetValue(int n) const;

	void SetFloat(const char *aKey, const float val);
	bool GetFloat(const char *aKey, float &val) const;
	bool GetDouble(const char *aKey, double &val) const;
	void SetInt(const char *aKey, const int val);
	bool GetInt(const char *aKey, int &val) const;
	void SetBool(const char *aKey, const bool val);
	bool GetBool(const char *aKey, bool &val) const;
	void SetString(const char *aKey, const CString &val);
	bool GetString(const char *aKey, CString &val) const;

	bool AppendAssoc(const char *aProp, const char *aVal);
};

#endif // TAGGEDLIST_H
