/*******************************************************************************
 *
 *   dataset.h -- 
 *
 *   Björn Nilsson, Medical Mathematics, 2005-2006
 */

#ifndef DATASET_H
#define DATASET_H

#include "types/claimable.h"
#include "types/ref.h"
#include "types/vector.h"
#include "types/matrix.h"
#include "annotation.h"
#include "types/nan.h"

class CDataset
{
public:
	CMatrix<float> m_Data;
	CVector<CAnnotation> m_vAnnotRows;
	CVector<CAnnotation> m_vAnnotCols;
	int m_nDisplayAnnotRow;
	int m_nDisplayAnnotCol;

public:
	// Construction etc
	CDataset() { m_nDisplayAnnotRow= 0; m_nDisplayAnnotCol= 0; }
	CRef<CDataset> CreateCopy() const;
	~CDataset() {}
	bool IsValid() const;

	// Data matrix
	const CMatrix<float> &GetDataRef() const { return *(const CMatrix<float> *)&m_Data; } 	
	bool HasData() const { return !m_Data.IsEmpty(); } //return m_pData!=NULL; }
	int GetDataWidth() const { return m_Data.GetWidth(); }
	int GetDataHeight() const { return m_Data.GetHeight(); }
	int GetWidth() const { return m_Data.GetWidth(); }
	int GetHeight() const { return m_Data.GetHeight(); }
	void GetDataSize(int &w, int &h) const { w= m_Data.GetWidth(); h= m_Data.GetHeight(); }
	void GetSize(int &w, int &h) const { w= m_Data.GetWidth(); h= m_Data.GetHeight(); }
	
	// Missing data
	bool IsNaN(int r, int c) const { return ::IsNaN(m_Data.GetAt(r,c)); }
	bool IsN(int r, int c) const { return ::IsN(m_Data.GetAt(r,c)); }
	int GetNaNCount() const;
	int GetNegativeCount() const;

	// Annotations
	void SetAnnotRows(const CVector<CAnnotation> &vRows) { m_vAnnotRows= vRows; }
	void SetAnnotCols(const CVector<CAnnotation> &vCols) { m_vAnnotCols= vCols; }
	const CVector<CAnnotation> GetAnnotCols() const { return m_vAnnotCols; }
	const CVector<CAnnotation> GetAnnotRows() const { return m_vAnnotRows; }
	const CAnnotation &GetAnnotCol(int nCol) const { return m_vAnnotCols[nCol]; }
	const CAnnotation &GetAnnotRow(int nRow) const { return m_vAnnotRows[nRow]; }
	int GetAnnotRowCount() const { return m_vAnnotRows.GetSize(); }
	int GetAnnotColCount() const { return m_vAnnotCols.GetSize(); }

	// Display
	int GetDisplayAnnotCol() const { return min(m_vAnnotCols.GetSize()-1, m_nDisplayAnnotCol); }
	int GetDisplayAnnotRow() const { return min(m_vAnnotRows.GetSize()-1, m_nDisplayAnnotRow); }
	int SetDisplayAnnotCol(int nCol) { m_nDisplayAnnotCol= nCol; }
	int SetDisplayAnnotRow(int nRow) { m_nDisplayAnnotRow= nRow; }
	int GetDisplayWidth() const;
	int GetDisplayHeight() const;
	void GetDisplaySize(int &w, int &h) const { w= GetDisplayWidth(); h= GetDisplayHeight(); }

	// Column and row identification 
	int FindAnnotRow(const CString &sName);
	int FindAnnotCol(const CString &sName);
	int ResolveAnnotCol(CString &sName);
	int ResolveAnnotRow(CString &sName);
};

#endif