/***********************************************************************************
 *
 *   unsignedset.h -- Set of unsigned
 *
 *   Erik Persson and Björn Nilsson, 1999-
 */

#ifndef UNSIGNEDSET_H
#define UNSIGNEDSET_H

#include "types/vector.h"

class CUnsignedSet
{
	CVector<bool> m_vbIncluded;
	int m_Magic;

private:
	void AssureSize(int size);

public:
	CUnsignedSet()
	{ m_Magic= 7777;}

	CUnsignedSet(const CUnsignedSet& other)
	{  m_Magic= 7777;
		*this= other; // assignment is implemented
	}

	~CUnsignedSet()
	{ ASSERT(m_Magic==7777);
	  m_Magic= 2222;
	}

	void Clear();
	void Include(int n, bool bIncluded);
	bool Includes(int n) const;
	void Invert(int nSize);

	const CUnsignedSet&
		operator=(const CUnsignedSet& other);
	const CUnsignedSet& operator&=
	  (const CUnsignedSet& other);
	const CUnsignedSet& operator|=
	  (const CUnsignedSet& other);
	const CUnsignedSet& operator^=
		(const CUnsignedSet& other);

	int GetFirst() const;
	int GetNext(int n) const;
	int GetCount() const;
};

#endif // UNSIGNEDSET_H
