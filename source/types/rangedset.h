/*******************************************************************************
 *
 *   Rangedset.h -- Specifies subsets of integers using half-open 
 *					subsets of the form [min, max)
 *
 *   Björn Nilsson, 2004
 */

#ifndef RANGEDSET_H
#define RANGEDSET_H

#include "types/vector.h"
#include "types/assert.h"

class CRangedSet
{
protected:
	struct subset 
	{ 
		int m_nMin; // Closed bound
		int m_nMax; // Open bound
	};
	CVector<subset> m_vSubsets;

	int GetLargestMin(int nValue) const;
	int GetInterval(int nValue) const;
public:
	CRangedSet() {};
	virtual ~CRangedSet() {};

	void Include(int nMin, int nMax);
	void Exclude(int nMin, int nMax);
	void Invert();
	
	void Empty() { m_vSubsets.SetSize(0); }
	bool IsEmpty() const { return m_vSubsets.IsEmpty(); }
	int GetMin() const;
	int GetMax() const;
	int GetCount() const;
	bool IsIncluded(int nValue) const { return GetInterval(nValue)>=0; }

	int GetSubsetCount() { return m_vSubsets.GetSize(); }
	void GetSubset(int nIndex, int &nMin, int &nMax)
	{ 
		ASSERT(nIndex>=0 && nIndex<GetSubsetCount());
		nMin= m_vSubsets[nIndex].m_nMin;
		nMax= m_vSubsets[nIndex].m_nMax;
	}
};

#endif
