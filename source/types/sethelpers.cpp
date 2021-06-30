//
//	sethelpers.cpp	--
//

#include "sethelpers.h"
#include "system/environment.h"
#include "system/filehelpers.h"
#include "types/tablefun.h"
#include "types/stringfun.h"

const char *ScanPathwayEntry_Panther(const char *ach0, const char **aref, int *nref)
{
	const char *ach1= ach0;
	while (*ach1 && *ach1!=';' && !(ach1[0]=='-' && ach1[1]=='>'))
		ach1++;
	*aref= ach0;
	*nref= ach1-ach0;

	if (*ach1==';')
		return ach1+1;

	if (ach1[0]=='-' && ach1[1]=='>')
	{
		// Skip trailing descriptor
		ach1 += 2;
		while (*ach1 && *ach1!=';')
			ach1++;
		while (*ach1==';')
			ach1++;
	}

	return ach1;
}

bool Sets_ScanCAT(const CTableIndex &idx_sco, int nAnnotRows_sco, int nKeyCol_sco,
				  const CTableIndex &idx_set, int nAnnotRows_set, int nKeyCol_set, const CVector<Table_SortItem_double> &vSI_set,
				  CVector<CString> &vGeneSetName, CVector< CVector<int> > &vGeneSet, int c)
{
	// Scan gene set definitions in the original CAT/TXT Rendercat format 

	const int C= idx_set.GetMinColCount();

	CVector<Table_SortItem_double> vSI;
	for (int r0=nAnnotRows_sco;r0<idx_sco.GetRowCount();r0++)
	{
		CString sID= Table_GetAt(idx_sco, r0, nKeyCol_sco);

		CVector<CString> vID;
		SplitAtChar(sID, '|', vID); // allow multi-id gene annotation

		int r1= -1;
		for (int i=0;r1<0 && i<vID.GetSize();i++)
			r1= Table_FindInIndex_double(vID[i], vSI_set);

		if (r1>=0)
		{
			CString s= Table_GetAt(idx_set, r1 + nAnnotRows_set, c);
			const char *ach0= s;
			int i= vSI.GetSize();

			while (true)
			{
				const char *aref; int nref;
				ach0= ScanPathwayEntry_Panther(ach0, &aref, &nref);
				if (nref<=0)
					break;

				// Gene already assigned to this group?
				// Corrects bug discovered by Petra Håkansson
				int j=i;
				for (;j<vSI.GetSize();j++)
					if (strncmp(aref, vSI[j].m_sKey, nref)==0)
						break;

				if (j==vSI.GetSize())
				{
					Table_SortItem_double si;
					si.m_sKey= CString(aref, nref);
					si.m_nIndex= r0 - nAnnotRows_sco; // Index w.r.t. data matrix
					ASSERT(si.m_nIndex>=0);
					vSI.Add(si);
				}
			}
		}
	}
	Table_SortItems_double(vSI, true, false);

	// Decode vSI
	for (int i=0;i<vSI.GetSize();)
	{
		CVector<int> vSet;
		vSet.Add(vSI[i].m_nIndex);
		CString sName= vSI[i].m_sKey;
		int i_end= i+1;
		for (;i_end<vSI.GetSize() && vSI[i_end].m_sKey==sName;i_end++) 
		{
			if (vSet.Find(vSI[i_end].m_nIndex)<0)
				vSet.Add(vSI[i_end].m_nIndex); // Add if not already member -- ensures uniqueness
		}

		if (vSet.GetSize()>1) // Robustness thing, protect against 1-gene groups (not proper groups)
		{
			vGeneSet.Add(vSet);
			vGeneSetName.Add(sName);
		}
		i= i_end;
	}

	return true;
}

CString GetTakitDir()
{
	// TODO: Should be useful in many apps. Move to ConsoleApp or Environment?
	CString sDir= Environment::GetEnv("TAKIT_HOME");
	if (sDir.IsEmpty())
		ReportWarning("Environment variable TAKIT_HOME has not been set", NULL);
	return sDir;
}

bool LoadSymbolData(CLoadBuf* &pLB, CTableIndex &ti)
{
	CString sFilename= g_FileSystem.GetCombinedName(GetTakitDir(), "external\\geneinfo\\Homo_sapiens_symbols_and_synonyms.txt");
	CString sErr;
	pLB= Table_LoadTabDataset(sFilename, ti, sErr);
	if (pLB!=NULL && ti.IsRectangular())
		return true;
	ReportError("Unable to open Homologene", sFilename);
	return false;
}

CString GeneIDfromSymbol(const CTableIndex &tiSymb, CString s, int cSymb, int cID)
{
	// This is slow, but returns the first listed hit in homo_sapiens_symbols_and_synonyms.
	// The synonyms part sometimes contains duplicates entries for the same symbol (higher, incorrect GeneIDs).
	// If faster version is implemented, must ensure the lowest GeneID is chosen.
	int r=0;
	for (;r<tiSymb.GetRowCount();r++)
	{
		int n= tiSymb.GetLastPointer(r, cSymb)-tiSymb.GetFirstPointer(r, cSymb);
		if (n==s.GetLength() && strncmp(s, tiSymb.GetFirstPointer(r, cSymb), n)==0)
			break;
	}

	if (r<tiSymb.GetRowCount())
		return Table_GetAt(tiSymb, r, cID);
	return "";
}

bool Sets_ScanAndTranslateGMT(const CTableIndex &idx_sco, int nAnnotRows_sco,
							  const CTableIndex &idx_set, int nAnnotRows_set,
							  CVector<CString> &vGeneSetName, CVector< CVector<int> > &vGeneSet)
{
	bool bOk= false;
	CLoadBuf *pLB;
	CTableIndex tiSymb;
	if (LoadSymbolData(pLB, tiSymb))
	{
		// Scan gene set definitions in Broad Insitute .GMT format
		const int cSymb= Table_FindAnnotCol(tiSymb, "Symbol", tiSymb.GetColCount(0));
		const int cID= Table_FindAnnotCol(tiSymb, "GeneID", tiSymb.GetColCount(0));
		if (cSymb>=0 && cID>=0)
		{
			for (int r=0;r<idx_set.GetRowCount();r++)
			{
				int C= idx_set.GetColCount(r);
				CString sSetName= Table_GetAt(idx_set, r, 0);
				if (!sSetName.IsEmpty())
				{
					CVector<int> vSet;
					printf("Translating set '%s'\n", (const char *)sSetName);
					for (int c=1;c<C;c++)
					{
						CString sSymb= Table_GetAt(idx_set, r, c);
						if (!sSymb.IsEmpty())
						{
							CString sID= GeneIDfromSymbol(tiSymb, sSymb, cSymb, cID); // TODO: This is SLOW! -- replace
							if (!sID.IsEmpty())
							{
								int i= Table_FindAnnotRow(idx_sco, sID, idx_sco.GetRowCount()); // TODO: This is SLOW! -- replace
								i -= nAnnotRows_sco; // indices refer to data matrix
								if (i>=0 && vSet.Find(i)<0)
								{
									// if (GetVerbosity()>=2)
									//	printf("Translated: %s --> %s\n", sSymb, sID);
									vSet.Add(i); // Add if valid and unique
								}
							}
							else
								printf("Not translated: %s\n", (const char *)sSymb);
						}
					}

					if (vSet.GetSize()>0)
					{
						vGeneSetName.Add(sSetName);
						vGeneSet.Add(vSet);
						// printf("Set added (n=%d): %s\n", vSet.GetSize(), sSetName);
					}
					else
						ReportWarning("Set not added", sSetName);
				}
			}

			bOk= true;
		}
		else
			ReportError("Invalid Symbol-GeneID translation table", NULL);

		g_FileSystem.FreeBuf(pLB);
	}
	return bOk;
}

bool Sets_ScanGMT(const CTableIndex &idx_sco, int nAnnotRows_sco, int nKeyCol_sco,
				  const CTableIndex &idx_set, int nAnnotRows_set, int nKeyCol_set,
				  CVector<CString> &vGeneSetName, CVector< CVector<int> > &vGeneSet)
{
	// Scan gene set definitions in Broad Insitute .GMT format
	printf("Scanning .gmt file...\n");

	// Prepare sorted search index
	CVector<Table_SortItem_double> vSI;
	for (int r=nAnnotRows_sco;r<idx_sco.GetRowCount();r++)
	{
		Table_SortItem_double si;
		si.m_sKey= Table_GetAt(idx_sco, r, nKeyCol_sco);
		si.m_nIndex= r;
		vSI.Add(si);
	}
	Table_SortItems_double(vSI, true, false);

	// 
	for (int r=0;r<idx_set.GetRowCount();r++)
	{
		int C= idx_set.GetColCount(r);
		CString sSetName= Table_GetAt(idx_set, r, 0);
		if (!sSetName.IsEmpty())
		{
			CVector<int> vSet;
			for (int c=1;c<C;c++)
			{
				int i=-1;
				CString s = Table_GetAt(idx_set, r, c);
				if (false)
				{
					// Linear search -- very slow, obsolete
					i= Table_FindAnnotRow(idx_sco, s, idx_sco.GetRowCount()); // TODO: This is SLOW! -- replace
				}
				else
				{
					// Binary search -- faster
					if (!s.IsEmpty())
					{
						int p_left= 0;
						int p_right= vSI.GetSize();
						while (p_left<p_right && i<0)
						{
							int p_mid= (p_left+p_right) >> 1;
							// printf("p_mid= %d, vSI[p_mid]= %s\n", p_mid, vSI[p_mid]);
							int cmp= strcmp(s, vSI[p_mid].m_sKey);
							if (cmp<0)
								p_right= p_mid;
							else if (cmp>0)
								p_left= p_mid+1;
							else
								i= vSI[p_mid].m_nIndex;
						}
					}
				}

				if (i>=0)
				{
					i -= nAnnotRows_sco; // Ensure i refers to data matrix, not table row
					if (vSet.Find(i)<0) // Ensure members remain unique
						vSet.Add(i);
				}
				else
				{
					if (!s.IsEmpty())
						printf("Not found: %s\n", (const char *)Table_GetAt(idx_set, r,c));
				}
			}

			// if (vSet.GetSize()>0)
			{
				vGeneSetName.Add(sSetName);
				vGeneSet.Add(vSet);
			}
		}
	}

	return true;
}


bool Sets_CheckIdentifiers(const CTableIndex &idx_sco, int nAnnotRows_sco, int nKeyCol_sco)
{
	// Assert that the gene identifiers in the score list are unique

	CVector<Table_SortItem_double> vSI(idx_sco.GetRowCount()-nAnnotRows_sco);
	for (int r=nAnnotRows_sco;r<idx_sco.GetRowCount();r++)
	{
		Table_SortItem_double si;
		si.m_sKey= Table_GetAt(idx_sco, r, nKeyCol_sco);
		si.m_nIndex= r;
		vSI.SetAt(r-nAnnotRows_sco, si);
	}
	Table_SortItems_double(vSI, true, false);

	for (int r=0;r<vSI.GetSize();)
	{
		int r_end;
		for (r_end=r+1;r_end<vSI.GetSize() && vSI[r_end].m_sKey==vSI[r].m_sKey;r_end++) {}

		if (r_end-r>1)
		{
			CString s= ::Format("Identifier '%s' occurs more than once (%d times) in gene list", vSI[r].m_sKey, r_end-r);
			ReportWarning(s, NULL);
		}
		r= r_end;
	}
	return true;
}

/******************************************************************************/

bool Sets_SaveGMT(const char *aFilename, 
				  const CTableIndex &idx_sco, int nAnnotRows_sco, int nKeyCol_sco, 
				  const CVector<CString> &vGeneSetName, const CVector< CVector<int> > &vGeneSet)
{
	FILE *pf= OpenOutputFile(aFilename);
	if (!pf)
		return false;

	bool bOk= true;
	for (int i=0;bOk && i<vGeneSetName.GetSize();i++)
	{
		fprintf(pf, vGeneSetName[i]);
		for (int j=0;j<vGeneSet[i].GetSize();j++)
			fprintf(pf, "\t%s", (const char *)Table_GetAt(idx_sco,  nAnnotRows_sco + vGeneSet[i][j], nKeyCol_sco));
		fprintf(pf, "\n");
	}

	(void) CloseFile(pf);
	return bOk;
}

void Sets_CreateKey(const CTableIndex &idx, int nAnnotRows, int nKeyCol, CVector<Table_SortItem_double> &vSI)
{
	ASSERT(vSI.GetSize()==0);
	(void) Table_GetKeyColumn(idx, nKeyCol, nAnnotRows, vSI);
	Table_SortItems_double(vSI, true, false);

	// Deduct number of annotation rows because internal 
	// routines use matrix indices, not table indices.
	Table_SortItem_double *pSI= vSI.GetBuffer(vSI.GetSize());
	for (int j=0;j<vSI.GetSize();j++)
		pSI[j].m_nIndex -= nAnnotRows;
}

bool Sets_Init_sco(const char *aFilename,
				   const CTableIndex &idx_sco, int &nAnnotRows_sco, int &nAnnotCols_sco, 
				   CString sKeyCol_sco, int &nKeyCol_sco,
				   CString sScoreCol_sco, int &nScoreCol_sco, CMatrix<double> &m_Data_sco)
{
	if (idx_sco.GetMinColCount()<2)
	{
		ReportError("Score table must contain at least two columns", aFilename);
		return false;
	}

	if (idx_sco.IsEmpty())
	{
		ReportError("Score table must be non-empty", aFilename);
		return false;
	}
	if (!idx_sco.IsRectangular())
	{
		ReportError("Score table must be rectangular", aFilename);
		return false;
	}

	// TODO: Merge InitIDX0 and InitDatasets. Works as is, but looks is structurally suboptimal
	nAnnotCols_sco= -1;
	nAnnotRows_sco= -1;
	CMatrix<short> mFlags;
	if (!Table_DetectAnnotation(idx_sco, nAnnotRows_sco, nAnnotCols_sco, mFlags))
	{
		ReportError("Unable to auto-detect number of annotations columns/rows", NULL);
		return false;
	}
	if (nAnnotCols_sco==idx_sco.GetMinColCount() || nAnnotRows_sco==idx_sco.GetRowCount())
	{
		ReportError("Score table only contains annotations, no numerical data", aFilename);
		return false;
	}

	if (!Table_ScanMatrix(idx_sco, nAnnotRows_sco, nAnnotCols_sco, idx_sco.GetRowCount()-1, idx_sco.GetMinColCount()-1, m_Data_sco, NULL))
	{
		ReportError("Illegal doubleing point format", aFilename);
		return false;
	}

	nKeyCol_sco= -1;
	if (sKeyCol_sco.IsEmpty())
	{
		// Nothing said, try 'GeneID'
		CString sGeneID= "GeneID";
		nKeyCol_sco= Table_ResolveCol(idx_sco, sGeneID);
		if (nKeyCol_sco>=0)
			ReportWarning("Assuming index column 'GeneID'", 0);
		else
		{
			ReportError("Index column not specified", sKeyCol_sco);
			return false;
		}
	}
	else
	{
		nKeyCol_sco= Table_ResolveCol(idx_sco, sKeyCol_sco); // TODO: Add switch to specify col name
		if (nKeyCol_sco<0)
		{
			ReportError("Column not found", sKeyCol_sco);
			return false;
		}
	}

	nScoreCol_sco= -1;
	if (sScoreCol_sco.IsEmpty())
	{
		if (nAnnotCols_sco>=idx_sco.GetMinColCount())
		{
			ReportError("Score table must contain at least one data column", 0);
			return false;
		}
		nScoreCol_sco= nAnnotCols_sco;
	}
	else
	{
		nScoreCol_sco= Table_ResolveCol(idx_sco, sScoreCol_sco);
		if (nScoreCol_sco<0)
		{
			ReportError("Column not found", sScoreCol_sco);
			return false;
		}
	}

	if (!Sets_CheckIdentifiers(idx_sco, nAnnotRows_sco, nKeyCol_sco))
		return false;

	return true;
}

CLoadBuf *Sets_Load_sco(const char *aFilename, 
						CTableIndex &idx_sco, int &nAnnotRows_sco, int &nAnnotCols_sco, 
						CString sKeyCol_sco, int &nKeyCol_sco, 
						CString sScoreCol_sco, int &nScoreCol_sco, 
						CVector<Table_SortItem_double> &vKey_sco, CMatrix<double> &m_Data_sco)
{
	// Load score table
	CString sErr;
	CLoadBuf *pLB= Table_LoadTabDataset(aFilename, idx_sco, sErr);
	if (!pLB)
	{
		ReportError(sErr, aFilename);
		return NULL;
	}
	if (!Sets_Init_sco(aFilename, idx_sco, nAnnotRows_sco, nAnnotCols_sco, sKeyCol_sco, nKeyCol_sco, sScoreCol_sco, nScoreCol_sco, m_Data_sco))
	{
		g_FileSystem.FreeBuf(pLB);
		return NULL;
	}

	Sets_CreateKey(idx_sco, nAnnotRows_sco, nKeyCol_sco, vKey_sco);
	return pLB;
}

CLoadBuf *Sets_Load_set(const char *aFilename,
						CTableIndex &idx_set, int &nAnnotRows_set, int &nAnnotCols_set, 
						int &nKeyCol_set, CVector<Table_SortItem_double> &vKey_set)
{
	// Load category definition table (RenderCat format or Broad .GMT format)
	CString sErr;
	CLoadBuf *pLB= Table_LoadTabDataset(aFilename, idx_set, sErr);
	if (!pLB)
	{
		ReportError(sErr, aFilename);
		return NULL;
	}
	nKeyCol_set= 0; // first column per default (both .cat and .gmt format)
	nAnnotCols_set= idx_set.GetMinColCount();
	nAnnotRows_set;
	CString sExt= g_FileSystem.GetFileExtension(aFilename);
	if (strcmp_nocase(sExt, "GMT")==0)
		nAnnotRows_set= 0;
	else 
		nAnnotRows_set= 1;

	// Check format
	bool bOk= false;
	if (idx_set.GetMinColCount()<1)
		ReportWarning("Gene set table must contain at least two columns", aFilename); // used to be error, now warning to allow empty gene sets
	// else
		bOk= true;

	// Unload if error
	if (!bOk)
	{
		g_FileSystem.FreeBuf(pLB);
		return NULL;
	}

	Sets_CreateKey(idx_set, nAnnotRows_set, nKeyCol_set, vKey_set);
	return pLB;
}