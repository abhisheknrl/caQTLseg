/***********************************************************************************
 *
 *   unsignedset.cpp -- Set of unsigned
 *
 *   Erik Persson and Björn Nilsson, 1999-
 */

#include "unsignedset.h"
#include "minmax.h"

/********************************/

void CUnsignedSet::AssureSize(int size)
{
	ASSERT(m_Magic==7777);
	
	int oldsize= m_vbIncluded.GetSize();
	if (oldsize >= size)
		return;

	int newsize= oldsize;
	if (newsize<32) newsize= 32;
	while (newsize<size)
		newsize= (newsize<<1);
	m_vbIncluded.SetSize(newsize);

	// Clear new bits
	for (int i= oldsize; i<newsize; i++)
		m_vbIncluded.SetAt(i, false);
}

void CUnsignedSet::Clear()
{
	int i, n= m_vbIncluded.GetSize();
	for (i= 0; i<n; i++)
		m_vbIncluded.SetAt(i, false);
}

void CUnsignedSet::Include(int n, bool bIncluded)
{
	AssureSize(n+1);
	ASSERT(n >= 0 && n < m_vbIncluded.GetSize());
	m_vbIncluded.SetAt(n, bIncluded);
}

bool CUnsignedSet::Includes(int n) const
{
	ASSERT(m_Magic==7777);
	if (n >= 0 && n < m_vbIncluded.GetSize())
		return m_vbIncluded[n];
	return false;
}

const CUnsignedSet&
CUnsignedSet::operator= (const CUnsignedSet& other)
{
	m_vbIncluded= other.m_vbIncluded;
	return (*this);
}

const CUnsignedSet&
CUnsignedSet::operator&= (const CUnsignedSet& other)
{
	int i, imax= m_vbIncluded.GetSize(), imax2= other.m_vbIncluded.GetSize();
	for (i= 0; i<imax; i++)
		if (i >= imax2 || !other.m_vbIncluded[i])
			m_vbIncluded.SetAt(i, false);
	return *this;
}

const CUnsignedSet&
CUnsignedSet::operator|= (const CUnsignedSet& other)
{
	AssureSize(other.m_vbIncluded.GetSize());
	int i, imax= __min(m_vbIncluded.GetSize(), other.m_vbIncluded.GetSize());
	for (i= 0; i<imax; i++)
		if (other.m_vbIncluded[i])
			m_vbIncluded.SetAt(i, true);
	return *this;
}

const CUnsignedSet&
CUnsignedSet::operator^= (const CUnsignedSet& other)
{
	AssureSize(other.m_vbIncluded.GetSize());

	int i, imax= __min(m_vbIncluded.GetSize(), other.m_vbIncluded.GetSize());
	for (i= 0; i<imax; i++)
		if (other.m_vbIncluded[i])
			m_vbIncluded.SetAt(i, !m_vbIncluded[i]);
	return *this;
}

void CUnsignedSet::Invert(int nSize)
{
	// Just as in mathematical set theory, a universe (nSize) is required to form the complement
	AssureSize(nSize);
	for (int i= 0; i<nSize; i++)
		m_vbIncluded.SetAt(i, !m_vbIncluded[i]);
}

int CUnsignedSet::GetFirst() const
{
	int i, imax= m_vbIncluded.GetSize();
	i= 0;
	while (i<imax)
	{	
		if (m_vbIncluded[i]) 
			return i;
		i++;
	}

	return -1; // not found
}

int CUnsignedSet::GetNext(int n) const
{
	// -1 tolerant

	int i= n, imax= m_vbIncluded.GetSize();
	i++;
	while (i<imax)
	{	
		if (m_vbIncluded[i]) 
			return i;
		i++;
	}

	return -1; // not found
}

int CUnsignedSet::GetCount() const
{
	// Slow implementation
	// TODO: Add a count field that makes this O(1)

	int imax= m_vbIncluded.GetSize();

	int cnt= 0;
	for (int i= 0; i<imax; i++)
	{
		if (m_vbIncluded[i])
			cnt++;
	}

	return cnt;
}

