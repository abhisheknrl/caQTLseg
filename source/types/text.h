/*******************************************************************************
 *
 *   Text.h -- Container class for huge texts.
 *
 *   Björn Nilsson, Sketchware, 2002
 */

#ifndef TEXTCONTAINER_H
#define TEXTCONTAINER_H

#include <afxwin.h>
#include "vector.h"

class CText  
{
public:
	CText();
	virtual ~CText();
	
	// Physical "lines", i.e. eol-delimited paragraphs.
	class CLine 
	{
	public:
		int m_nOffs;
		CString m_sText;
		// loword=token offset, hiword= token value
		CVector<unsigned long> m_vTokens; 

		CLine() {};
		CLine(const CLine &l)
		{
			m_nOffs= l.m_nOffs;
			m_sText= l.m_sText;
			m_vTokens= l.m_vTokens;
		};
		bool IsEmpty() const
		{	return m_sText.IsEmpty()!=0; }
		int GetLength() const
		{	return m_sText.GetLength(); }
		int Delete( int nIndex, int nCount = 1 )
		{	return m_sText.Delete(nIndex, nCount); }
		int Insert( int nIndex, LPCTSTR pstr )
		{	return m_sText.Insert(nIndex, pstr); }
		CString Mid( int nFirst, int nCount ) const
		{	return m_sText.Mid(nFirst, nCount);	}
		CString Mid( int nFirst ) const
		{	return m_sText.Mid(nFirst); }
		CString Left( int nCount ) const
		{	return m_sText.Left(nCount); }
		void AddToken(int nValue, int nOffset)
		{	m_vTokens.Add((nValue << 16) + (nOffset & 0xffff)); }
	};

protected:
	CVector<CLine> m_vLines;
	int m_nSelStart;
	int m_nSelEnd;
	int m_nDirtyMin;
	int m_nDirtyMax;
	int m_nDirtyDelta;
public:
	bool MakeNonEmpty();
	int LineFromOffset(int nOffs) const;
	int NormalizeLine(int nLine);
	int GetLineCount() const
	{	return m_vLines.GetSize(); }
	int GetTextSize() const;
	CLine *GetBuffer(int nSize)
	{	return m_vLines.GetBuffer(nSize); }
	const CLine& GetLine(int nLine) const
	{	return m_vLines[nLine]; }
	void SetAt(int nPos, const CLine &line)
	{	m_vLines.SetAt(nPos, line); }
	void UpdateOffset(int nFromLine);
	CString GetString() const;

	void GetSel(int &nStart, int &nEnd) const
	{
		nStart= m_nSelStart;
		nEnd= m_nSelEnd;
	}
	void SetSel(int nStart, int nEnd)
	{
		m_nSelStart= max(0, nStart);
		m_nSelEnd= min(nEnd, GetTextSize());
	}
	void ClipSelection(int &nStart, int &nEnd);
	void GetDirtyRange(int &nMin, int &nMax) const
	{	nMin= m_nDirtyMin; nMax= m_nDirtyMax; }
	int GetDirtyDelta()
	{	return m_nDirtyDelta; }
	void ClearDirty()
	{	m_nDirtyMin= m_nDirtyMax= -1; m_nDirtyDelta= 0; }
	void GetTokens(int nL0, int nL1, CVector<unsigned long> &vTokens);

	// Virtual so that undo functionality can be incorporated at low level.
	virtual bool ReplaceSel(CString sText);
	void CleanEol(CString &s);
};

#endif
