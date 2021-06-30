/*******************************************************************************
 *
 *   dataset.cpp -- 
 *
 *   Björn Nilsson, Medical Mathematics, 2005-2006
 */

#include "dataset.h"
#include "stringfun.h"

/******************************************************************************/
// Construction

/*
CRef<CDataset> CDataset::CreateCopy() const
{
	ASSERT(this);

	CRef<CDataset> rCopy;
	try
	{
		rCopy.CreateNew();
	}
	catch(...)
	{
		return NULL;
	}

	CDataset &copy = rCopy.GetExclusive();
	copy.m_pData= m_pData;
	copy.m_vAnnotRows= m_vAnnotRows;
	copy.m_vAnnotCols= m_vAnnotCols;
	copy.m_nDisplayAnnotCol= m_nDisplayAnnotCol;
	copy.m_nDisplayAnnotRow= m_nDisplayAnnotRow;

	return rCopy;
}

bool CDataset::IsValid() const
{
	if (m_pData)
	{
		// Has data
		int w= m_Data.GetWidth();
		for (int i=0;i<m_vAnnotRows.GetSize();i++)
			if (m_vAnnotRows[i].m_vValues.GetSize()!=w)
				return false; // Annotation row has wrong length
		int h= m_Data.GetHeight();
		for (i=0;i<m_vAnnotCols.GetSize();i++)
			if (m_vAnnotCols[i].m_vValues.GetSize()!=h)
				return false; // Annotation col has wrong height 
	}
	else
	{
		// No data means only annotation columns allowed
		if (m_vAnnotRows.GetSize()>0)
		{
			int w= m_vAnnotRows[0].m_vValues.GetSize();
			for (int i=1;i<m_vAnnotRows.GetSize();i++)
				if (m_vAnnotRows[i].m_vValues.GetSize()!=w)
					return false; // Annotation rows do not have equal lengths
		}
		if (m_vAnnotCols.GetSize()>0)
			return false;
	}

	return true;
}
*/

/******************************************************************************/
// Missing and negative data

int CDataset::GetNaNCount() const
{
	// Returns the no of NaN floats (missing data) in the matrix
	int n= 0;
	if (HasData())
		for (int r=0;r<m_Data.GetHeight();r++)
			for (int c=0;c<m_Data.GetWidth();c++)
				if (IsNaN(r,c))
					n++;
	return n;
}

int CDataset::GetNegativeCount() const
{
	// Returns the no of negative elements in the matrix
	int n= 0;
	if (HasData())
		for (int r=0;r<m_Data.GetHeight();r++)
			for (int c=0;c<m_Data.GetWidth();c++)
				if (IsN(r,c) && m_Data.GetAt(r,c)<0)
					n++;
	return n;
}

/******************************************************************************/
// Display

int CDataset::GetDisplayWidth() const
{
	if (HasData())
		return m_Data.GetWidth();
	else if (m_vAnnotRows.GetSize()>0)
		return m_vAnnotRows[0].m_vValues.GetSize();
	return 0;
}

int CDataset::GetDisplayHeight() const
{
	if (HasData())
		return m_Data.GetHeight();
	else if (m_vAnnotCols.GetSize()>0)
		return m_vAnnotCols[0].m_vValues.GetSize();
	return 0;
}

int CDataset::FindAnnotRow(const CString &sName)
{
	for (int i=0;i<m_vAnnotRows.GetSize();i++)
		if (m_vAnnotRows[i].m_sTitle==sName)
			return i;
	return -1;
}

int CDataset::FindAnnotCol(const CString &sName)
{
	for (int i=0;i<m_vAnnotCols.GetSize();i++)
		if (m_vAnnotCols[i].m_sTitle==sName)
			return i;
	return -1;
}

int CDataset::ResolveAnnotCol(CString &sName)
{
	// Returns: Index of column sName or -1. Translates sName if necessary.
	int c;
	const char *ach0= sName;
	const char *ach1= Scan_UnsignedInt(sName, &c);
	if (ach1==ach0+sName.GetLength() && c>=0 && c<m_vAnnotCols.GetSize())
	{
		sName= m_vAnnotCols[c].m_sTitle;
		return c;
	}
	return FindAnnotCol(sName);
}

int CDataset::ResolveAnnotRow(CString &sName)
{
	// Returns: Index of row sName or -1. Translates sName if necessary.
	int r;
	const char *ach0= sName;
	const char *ach1= Scan_UnsignedInt(sName, &r);
	if (ach1==ach0+sName.GetLength() && r>=0 && r<m_vAnnotRows.GetSize())
	{
		sName= m_vAnnotRows[r].m_sTitle;
		return r;
	}
	return FindAnnotRow(sName);
}
