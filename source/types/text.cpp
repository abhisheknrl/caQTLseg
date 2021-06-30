/*******************************************************************************
 *
 *   Text.cpp -- 
 *
 *   Björn Nilsson, 2002
 */

#include "text.h"

CText::CText()
{
	m_nSelStart= 0;
	m_nSelEnd= 0;
	m_nDirtyMin= -1;
	m_nDirtyMax= -1;
	(void) MakeNonEmpty();
}

CText::~CText()
{
}

int CText::NormalizeLine(int nLine)
{
	// Kernel function that splits a line at all '\n' except terminating ones.
	// Returns index of last line.
	CVector<CLine> vNew;
	const CString &s= m_vLines[nLine].m_sText;

	int nBaseOffs= m_vLines[nLine].m_nOffs;

	const char *aText= s;
	const char *aStart= aText;
	const char *aEnd= aStart + s.GetLength() - 1; // -1 in order to skip the last '\n'
	while (aText<aEnd)
	{
		if (*aText=='\n')
		{
			CLine l0;
			
			int nLen= aText-aStart+1;
			l0.m_sText= CString(aStart, nLen);
			// TokenizeLine(l0);
			vNew.Add(l0);

			aStart= aText + 1;
		}

		aText++;
	}

	if (vNew.GetSize()>0)
	{
		CLine l0;
		l0.m_sText= CString(aStart, (LPCTSTR)s + s.GetLength() - aStart);

		// TokenizeLine(l0);
		m_vLines.SetAt(nLine, l0);
		m_vLines.InsertVector(nLine, vNew);
	}
	else
	{
		CLine l0= m_vLines[nLine];
		// TokenizeLine(l0);
		m_vLines.SetAt(nLine, l0);
	}

	return nLine + vNew.GetSize();
}

int CText::LineFromOffset(int nOffs) const
{
	// Finds the line containing nOffs using binary search
	int first= 0;
	int last= m_vLines.GetSize()-1;

	if (first>last) 
		return 0;
	if (nOffs>m_vLines[last].m_nOffs + m_vLines[last].GetLength())
		return last;

	while (first<last)
	{
		int i= (first+last)>>1; // ranges from first..last-1
		const CLine &l= m_vLines[i];
		if (l.m_nOffs + l.GetLength() > nOffs)
			last= i;
		else
			first= i+1;
	}

	return first;
}

int CText::GetTextSize() const
{
	int n= m_vLines.GetSize()-1;
	if (n>=0)
		return m_vLines[n].m_nOffs + m_vLines[n].GetLength();
	else
		return 0;
}

bool CText::MakeNonEmpty()
{
	// Helper function to ReplaceSel
	// Returns true if text was empty.
	if (m_vLines.IsEmpty())
	{
		CLine l;
		l.m_nOffs= 0;
		m_vLines.Add(l);
		return true;
	}

	return false;
}

void CText::UpdateOffset(int nFromLine)
{
	int nOffs= 0;
	if (nFromLine>0)
		nOffs= m_vLines[nFromLine-1].m_nOffs + m_vLines[nFromLine-1].GetLength();

	int L= m_vLines.GetSize();
	CLine *pLine= m_vLines.GetBuffer(L);
	CLine *pLineEnd= &pLine[L];
	pLine = &pLine[nFromLine];
	while (pLine<pLineEnd)
	{
		pLine->m_nOffs= nOffs;
		nOffs += pLine->GetLength();
		pLine++;
	}
}

CString CText::GetString() const
{
	// Retrieve the whole text as a CString
	
	int nSize= GetTextSize();
	CString s;	
	char *aText= s.GetBuffer(nSize);

	int i= 0;
	int I= m_vLines.GetSize();
	while (i<I)
	{
		const CLine &line= m_vLines[i];
		(void) strncpy(aText, line.m_sText, line.GetLength());
		aText += line.GetLength();
		i++;
	}

	s.ReleaseBuffer(nSize);
	return s;
}

void CText::ClipSelection(int &nStart, int &nEnd)
{
	if (nStart<0)
		nStart=0;
	if (nEnd>GetTextSize())
		nEnd= GetTextSize();
}

bool CText::ReplaceSel(CString sText)
{
	ClipSelection(m_nSelStart, m_nSelEnd);
	if (m_nSelStart>m_nSelEnd)
	{
		ASSERT(false);
		return false;
	}

	MakeNonEmpty();

	int nOldSize= GetTextSize();
	int l0= LineFromOffset(m_nSelStart);
	int l1= LineFromOffset(m_nSelEnd);

	bool bOk= false;
	try
	{
		// Replace selection
		CLine line0= m_vLines[l0];
		int n= m_nSelStart-m_vLines[l0].m_nOffs;
		if (l0==l1)
		{
			if (m_nSelEnd>m_nSelStart)
				line0.Delete(n, m_nSelEnd-m_nSelStart);
			if (!sText.IsEmpty())
				line0.Insert(n, sText);
		}
		else
		{
			// l1!=l0 => m_nSelEnd>m_nSelStart
			const CLine &line1= m_vLines[l1];
			line0.m_sText= line0.Left(n) + sText + line1.Mid(m_nSelEnd-line1.m_nOffs);
		}
		m_vLines.SetAt(l0, line0);

		// Delete emptied lines
		if (l1>l0)
			m_vLines.Delete(l0 + 1, l1 - l0);
		ASSERT(m_vLines.GetSize()>l0);

		bOk= true;
	}
	catch(...)
	{
	}

	// Normalize line.
	m_nDirtyMin= l0;
	m_nDirtyMax= NormalizeLine(l0);

	// Update offsets
	UpdateOffset(l0);

	// Adjust selection
	m_nDirtyDelta= GetTextSize() - nOldSize;
	SetSel(m_nSelStart, m_nSelEnd + m_nDirtyDelta);

	return bOk;
}

void CText::CleanEol(CString &s)
{
	// Helper function for converting 13-10 EOLs to 10-EOLs.
	int L= s.GetLength();
	char *aRead= s.GetBuffer(L);
	char *aWrite= aRead;
	char *aEnd= aRead + L;
	while (aRead<aEnd)
	{
		if (*aRead==13)
			L--;
		else
			*aWrite++= *aRead;

		aRead++;
	}

	s.ReleaseBuffer(L);
}

void CText::GetTokens(int nL0, int nL1, CVector<unsigned long> &vTokens)
{
	int i=nL0;
	while (i<=nL1)
	{
		if (m_vLines[i].m_vTokens.GetSize())
			vTokens.Append(m_vLines[i].m_vTokens);
		i++;
	}
}
