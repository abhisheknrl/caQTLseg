/*******************************************************************************
 *
 *   RangedSet.cpp -- Implementation of CRangedSet
 *
 *   Björn Nilsson, 2004
 */

#include "rangedset.h"
#include <limits.h>

int CRangedSet::GetLargestMin(int nValue) const
{
	// TODO: Convert to binary search for better performance.
	int i=-1;
	while (i+1<m_vSubsets.GetSize() && m_vSubsets[i+1].m_nMin<=nValue)
		i++;
	return i;
}

int CRangedSet::GetInterval(int nValue) const
{
	// Returns index of the interval that contains nValue, 
	// or -1 if nValue is not included in the set.
	int i= GetLargestMin(nValue);
	if (i>=0 && nValue<m_vSubsets[i].m_nMax)
		return i;
	return -1;
}

int CRangedSet::GetMin() const 
{
	// Get smallest item (infinity if empty)
	if (m_vSubsets.IsEmpty())
		return INT_MAX;
	return m_vSubsets[0].m_nMin;
}

int CRangedSet::GetMax() const 
{
	// Get largest item (-infinity if empty)
	int i= m_vSubsets.GetSize()-1;
	if (i<0)
		return INT_MIN;
	return m_vSubsets[i].m_nMax-1;
}

int CRangedSet::GetCount() const 
{
	// Get number of elements in the set.
	int n= 0;
	for (int i=m_vSubsets.GetSize();--i>=0;)
		n += m_vSubsets[i].m_nMax - m_vSubsets[i].m_nMin;
	return n;
}

void CRangedSet::Include(int nMin, int nMax)
{
	if (nMax<nMin)
		return;

	CVector<subset> vNew;
	int n= m_vSubsets.GetSize();

	// Keep preceeding intervals.
	int i=0;
	while (i<n && m_vSubsets[i].m_nMax<nMin)
		vNew.Add(m_vSubsets[i++]);
	if (i>=0 && i<n && m_vSubsets[i].m_nMin<nMin)
		nMin= m_vSubsets[i].m_nMin;

	// Skip contained intervals.
	int j= i;
	while (j<n && nMax>=m_vSubsets[j].m_nMax)
		j++; 
	if (j<n && nMax>=m_vSubsets[j].m_nMin)
	{
		nMax= m_vSubsets[j].m_nMax;
		j++;
	}

	// Add interval containing new inclusion.
	subset sNew;
	sNew.m_nMin= nMin;
	sNew.m_nMax= nMax;
	vNew.Add(sNew);

	// Keep succeeding intervals.
	while (j<n)
		vNew.Add(m_vSubsets[j++]);
	m_vSubsets= vNew;
}

void CRangedSet::Exclude(int nMin, int nMax)
{
	ASSERT(false); // TODO: Implement Exclude.
}

void CRangedSet::Invert()
{
	subset s;
	CVector<subset> vNew;

	if (IsEmpty())
	{
		s.m_nMin= INT_MIN;
		s.m_nMax= INT_MAX;
		vNew.Add(s);
	}
	else
	{
		if (m_vSubsets[0].m_nMin!=INT_MIN)
		{
			s.m_nMin= INT_MIN;
			s.m_nMax= m_vSubsets[0].m_nMin;
			vNew.Add(s);
		}

		int n= m_vSubsets.GetSize();
		for (int i=0;i<n-1;)
		{
			s.m_nMin= m_vSubsets[i].m_nMax;
			i++;
			s.m_nMax= m_vSubsets[i].m_nMin;
			vNew.Add(s);
		}

		if (m_vSubsets[n-1].m_nMax!=INT_MAX)
		{
			s.m_nMin= m_vSubsets[n-1].m_nMax;
			s.m_nMax= INT_MAX;
			vNew.Add(s);
		}		
	}

	m_vSubsets= vNew;
}

