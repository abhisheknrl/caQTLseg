/*******************************************************************************
 *
 *   caQTLseg.exe - Segmentation of chromatin accessibility (caQTL) data
 *
 *   Björn Nilsson and Abhishek Niroula, 2020-
 *
 */

//
// input: 
//
//		-	filelist: files with cut site coordinates in the region-of-interest
//		-	filefolder: if you wish to replace pre-pend a different folder than what is in the file list
//		-	genotypefile: file with genotypes for each sample
//		-	rsID
//		-	minbp, maxbp: defines the interval to segment (any cutsite outside it will be excluded)\
//		-	subsampling
//		-	lambda: penalty
//

#include "types/tablefun.h"
#include "system/consoleapp.h"
#include "system/filesystem.h"
#include "types/stringfun.h"
#include "math/statfun.h"
#include "math/variance.h"
#include <float.h>
#include "types/minmax.h"
#include "math/rand.h"
#include "system/linereader.h"
#include <stdint.h>

struct SegItem
{
	int m_HeapPos;
	double m_cost; // best cost so far
	int m_state; // 
	int m_PreviousIndex; // index of previous segment in vector 
	int m_leftbp; // bp position on the left side of this segitem (cell), this is to simplify downsampling
};

int g_nRefs= 0;

/******************************************************************************/

class CApp : public CConsoleApp
{
protected:
	bool m_bHelpShown;

	//
	CString m_sFileFolder;
	CString m_sRSID;
	int m_nMinBP, m_nMaxBP;
	
	//
	CString m_sCHR;
	int m_nBinSpacing;
	int m_nBinSize;

	CVector<double> m_vLambda1, m_vLambda2;
	CVector<CString> m_vsLambda1, m_vsLambda2;

	int m_nMaxSegmentLength_state1;

	bool m_bMulti;
	int m_nMultiWW0, m_nMultiWW1;
	bool m_bMultiQuantileNormalize;
	bool m_bMultiOverwriteCache;

	int m_nPermutations;
	bool m_bExportIGV;
	bool m_bExportSEG;
	bool m_bExportPI0;
	CString m_sNaN;

	bool Optimize(const double lambda1, const double lambda2, CVector<double>& vG, CMatrix<int>& mCounts, int nMinBP, int nMaxBP, int nBinSpacing, const CString& sBasename, int& nBPstate1);

	//
	CTableIndex m_idx_genotypes;

	void DisplayHelp();
	void SetDefaultParameters();
	bool CheckParameters();
	virtual bool OnSwitch(const char *ach, CString &sSwitch);

	//
	bool ParseGenotypes(const char *aFilename, CVector<double>& vG, CString &sA1, CString &sA2, CVector<CString>& vSampleID, int& n00, int& n01, int& n11, int& nA);
	bool ParseGenotypes_multi(const char* aFilename, CVector<CString>& vSampleID, CMatrix<double>& mG, CVector<CString>& vRSID, CVector<CString>& vCHR, CVector<int>& vBP, CVector<CString>& vA1, CVector<CString>& vA2, CVector<int>& vn00, CVector<int>& vn01, CVector<int>& vn11, CVector<int>& vnNA);

	bool LoadCutSiteFiles(const char* aFolderName, const CVector<CString>& vSampleID, CMatrix<int>& mCounts, const CString& sCHR, int nMinBP, int nMaxBP, int nBinSpacing);
	bool LoadCutSiteFiles_multi(const char* aFolderName, const CVector<CString>& vSampleID, const CVector<CString>& vSNP_CHR, CVector<int>& vSNP_BP, const int& ww0, const int& ww1, CMatrix<double>& mData_out, CMatrix<double>& mG);
	bool BED2BIN(const char* aFilenameBED, const char* aFilenameBIN, const CVector<CString>& vsCHRdictionary, const bool bOverwriteCache);

	bool DumpCounts(const char* aFilename, CMatrix<double>& mCounts_norm, CString& sChr, int nMinBP, int nMaxBP, int nBinSpacing);
	bool ExportIGV(const CString& sBasename, CMatrix<double>& mCounts_norm, const double* pG, const double* pX, const double* pY);
	bool ExportPI0(const CString& sFilename, CVector<double>& vLambda1, CVector<double>& vLambda2, CMatrix<double>& mPI0, CMatrix<double>& mH1proportion);

public:
	virtual int Main(int argc, TCHAR* argv[], TCHAR* envp[]);
};

//
// CCommand overrides

void CApp::DisplayHelp()
{
	m_bHelpShown=true;

	CString sName= g_FileSystem.GetExecutableFileName();
	sName= g_FileSystem.GetNameWithoutExtension(g_FileSystem.GetStrippedName(sName));
	printf("Usage: %s <genotype file> <cut files folder> <output basename> [options]\n\n", (const char *)sName);
	printf("Options:\n\n");
	printf("-rsID\t\tSet genetic marker.\n");
	printf("-chr\t\tSet chromosome of region-of-interest.\n");
	printf("-minbp\t\tSet start of region-of-interest in base pairs.\n");
	printf("-maxbp\t\tSet end of region-of-interest in base pairs.\n");
	printf("-binspacing\tSet bin spacing in base pairs.\n");
	printf("-binwidth\tSet bin size in base pairs.\n");
	printf("-l1\t\tSet lambda1 penalty. Use comma to separate mutiple values.\n");
	printf("-l2\t\tSet lambda2 penalty. Use comma to separate mutiple values.\n");
	printf("-p\t\tSet number of permutations used to estimate pi0.\n");
	printf("-seg\t\tExport .seg file containing a list of segments with allele-dependent accessibility.\n");
	printf("-pi0\t\tExport .pi0 file with pi0 estimates for each combination of (lambda1, lambda2) values. Useful if multiple lambda1 and lambda2 are being used.\n");
	printf("-igv\t\tExport .igv file for Integrative Genomics Viewer.\n");
	printf("-?\t\tDisplay help text.\n\n");
}

bool CApp::OnSwitch(const char *ach, CString &sSwitch)
{
	if (sSwitch == "?")
	{
		DisplayHelp();
		return true;
	}
	else if (sSwitch == "minbp")
	{
		if (Table_CheckIntFormat_complete(ach) != 0)
		{
			m_nMinBP = atoi(ach);
			if (m_nMinBP < 0)
			{
				ReportError("-minbp must be non-negative", 0);
				return false;
			}
			return true;
		}
	}
	else if (sSwitch == "maxbp")
	{
		if (Table_CheckIntFormat_complete(ach) != 0)
		{
			m_nMaxBP = atoi(ach);
			if (m_nMaxBP < 0)
			{
				ReportError("-maxbp must be non-negative", 0);
				return false;
			}
			return true;
		}
	}
	else if (sSwitch == "rsid" || sSwitch=="rsID" || sSwitch=="snp")
	{
		m_sRSID = ach;
		return true;
	}
	else if (sSwitch == "l1")
	{
		CVector<CString> v;
		SplitAtChar(ach, ',', v);
		int i = 0;
		for (; i < v.GetSize(); i++)
		{
			if (Table_CheckFloatFormat_complete(v[i]) == 0)
				break;
			m_vLambda1.Add(atof(v[i]));
			m_vsLambda1.Add(v[i]);
		}
		if (i == v.GetSize())
			return true;
	}
	else if (sSwitch== "l2")
	{
		CVector<CString> v;
		SplitAtChar(ach, ',', v);
		int i = 0;
		for (; i < v.GetSize(); i++)
		{
			if (Table_CheckFloatFormat_complete(v[i]) == 0)
				break;
			m_vLambda2.Add(atof(v[i]));
			m_vsLambda2.Add(v[i]);
		}
		if (i == v.GetSize())
			return true;
	}
	else if (sSwitch == "binspacing")
	{
		if (Table_CheckIntFormat_complete(ach) != 0)
		{
			m_nBinSpacing = atoi(ach);
			if (m_nBinSpacing < 10)
			{
				ReportError("Please specify a reasonable bin spacing, check -binspacing switch", 0);
				return false;
			}
			return true;
		}
	}
	else if (sSwitch == "binsize")
	{
		if (Table_CheckIntFormat_complete(ach) != 0)
		{
			m_nBinSize = atoi(ach);
			if (m_nBinSize < 10)
			{
				ReportError("Please specify a reasonable bin spacing, check -binsize switch", 0);
				return false;
			}
			return true;
		}
	}
	else if (sSwitch == "chr")
	{
		m_sCHR = ach;
		return true;
	}
	else if (sSwitch == "multi")
	{
		m_bMulti = true;
		return true;
	}
	else if (sSwitch == "ww0")
	{		
		if (Table_CheckIntFormat_complete(ach) != 0)
		{
			m_nMultiWW0 = atoi(ach);
			return true;
		}
	}
	else if (sSwitch == "ww1")
	{
		if (Table_CheckIntFormat_complete(ach) != 0)
		{
			m_nMultiWW1 = atoi(ach);
			return true;
		}
	}
	else if (sSwitch == "n")
	{
		m_bMultiQuantileNormalize = true;
		return true;
	}
	else if (sSwitch == "o")
	{
		m_bMultiOverwriteCache = true;
		return true;
	}
	else if (sSwitch == "p")
	{
		if (Table_CheckIntFormat_complete(ach) != 0)
		{
			m_nPermutations = atoi(ach);
			if (m_nPermutations < 10 || m_nPermutations>10000)
			{
				ReportError("Please specifiy a reasonable number of permutations, check -p switch", 0);
				return false;
			}
			return true;
		}
	}
	else if (sSwitch == "igv")
	{
		// generate .igv files (one per l1,l2 pair)
		m_bExportIGV = true;
		return true;
	}
	else if (sSwitch == "seg")
	{
		// generate .seg files (one per l1,l2 pair)
		m_bExportSEG = true;
		return true;
	}
	else if (sSwitch == "pi0")
	{
		// generate .pi0 file with the estimated pi0 values for all l1,p2 pairs (one l1 per row, one l2 per column; requires -p set)
		m_bExportPI0 = true;
		return true;
	}
	else if ((sSwitch == "NaN") || (sSwitch = "nan") || (sSwitch == "nanstring") || (sSwitch == "NaNstring"))
	{
		// text used in output files for NaN values
		m_sNaN= ach;
		return true;
	}

	CString s= sSwitch;
	if (ach && *ach)
		s += ":" + CString(ach);
	ReportError("Illegal switch", (const char *)s);
	return false;
}

void CApp::SetDefaultParameters()
{
	m_bHelpShown= false;
	m_nMinBP = -1;
	m_nMaxBP = -1;
	m_nBinSpacing = -1;
	m_nBinSize = -1;
	m_nMaxSegmentLength_state1 = 10000;
	m_nPermutations = -1;
	m_bMulti = false;
	m_nMultiWW0 = -1;
	m_nMultiWW1 = -1;
	m_bMultiQuantileNormalize = false;
	m_bMultiOverwriteCache = false;
	m_bExportIGV = false;
	m_bExportSEG = false;
	m_bExportPI0 = false;
	m_sNaN = "n/a";
}

bool CApp::CheckParameters()
{
	if (m_bHelpShown)
		return false;

	if (!m_bMulti)
	{
		// single-region (standard) mode
		if (m_vLambda1.IsEmpty())
		{
			ReportError("-l1 has not been set", 0);
			return false;
		}
		for (int i=0;i<m_vLambda1.GetSize();i++)
		{
			if (m_vLambda1[i] < 1.0)
			{
				ReportError("-l1 must be greater than 1.0", ::Format("%g", m_vLambda1[i]));
				return false;
			}
		}
		if (m_vLambda2.IsEmpty())
		{
			ReportError("-l2 has not been set", 0);
			return false;
		}
		for (int i = 0; i < m_vLambda2.GetSize(); i++)
		{
			if (m_vLambda2[i] < 0.0)
			{
				ReportError("-l2 must be non-negative", ::Format("%g", m_vLambda2[i]));
				return false;
			}
		}
		if (m_bExportPI0 && m_nPermutations < 0)
		{
			ReportError("to use -pi0, you need to specify the number of permutation used to estimate the null distribution, use -p switch", 0);
			return false;
		}
		if (m_sCHR.IsEmpty())
		{
			ReportError("-chr has not been set", 0);
			return false;
		}
		if (m_nMinBP < 0)
		{
			ReportError("-minbp has not been set", 0);
			return false;
		}
		if (m_nMaxBP < 0)
		{
			ReportError("-minbp has not been set", 0);
			return false;
		}
		if (m_nMinBP >= m_nMaxBP)
		{
			ReportError("-minbp must be less than -maxbp", 0);
			return false;
		}
		if (m_nBinSpacing < 0)
		{
			ReportError("-binspacing not set", 0);
			return false;
		}
		if (m_nBinSize < 0)
		{
			ReportError("-binsize not set", 0);
			return false;
		}
		if ((m_nBinSize % m_nBinSpacing) != 0 || ((m_nBinSize/m_nBinSpacing) % 2==0))
		{
			ReportError("-binsize must be an odd multiple of -binSpacing",0);
			return false;
		}
		if (m_sRSID.IsEmpty())
		{
			ReportError("-rsid has not been set", 0);
			return false;
		}
	}
	else
	{
		// multi-region mode
		if (m_nMultiWW0 < 0)
		{
			ReportError("-ww0 not set", 0);
			return false;
		}
		if (m_nMultiWW1 < 0)
		{
			ReportError("-ww1 not set", 0);
			return false;
		}
		if (m_nMultiWW0 > m_nMultiWW1)
		{
			ReportError("-ww0 must be smaller than -ww1", 0);
			return false;
		}
	}

	return true;
}

//
// helper functions

void PQfromZ(const CVector<double> &vZ, double *pZ0, const int &nZ0, CVector<double> &vP, CVector<double> &vQ)
{
	const int nZ= vZ.GetSize();
	
	if (pZ0 == 0)
	{
		// no null distribution
		vP.SetSize(nZ);
		vP.Fill(Stat_GetNaN_double());
		vQ.SetSize(nZ);
		vQ.Fill(Stat_GetNaN_double());
		return;
	}

	// sort z scores
	CVector<Stat_SortItem_double> vSI(nZ);
	Stat_SortItem_double *pSI= vSI.GetBuffer();
	for (int i=0;i<nZ;i++)
	{
		Stat_SortItem_double si;
		si.m_Value= __abs(vZ[i]);
		si.m_nIndex= i;
		vSI.SetAt(i, si);
	}
	Stat_QuicksortIndexed_double(pSI, nZ, true);

	// sort null distribution of z scores
	for (int i=0;i<nZ0;i++)
		pZ0[i]= __abs(pZ0[i]);
	Stat_Quicksort_double(pZ0, nZ0, true);

	// calculate false discovery rate (q-value) for each z
	vQ.SetSize(nZ);
	vP.SetSize(nZ);
	double pi0= 1.0; // no pi0 estimation implemented at the moment (but close to 1.0 for most applications)
	double qmin= pi0; 
	int i0= 0;
	for (int i=0;i<nZ;i++)
	{
		while (i0<nZ0 && pZ0[i0]<=pSI[i].m_Value)
			i0++;

		int j= pSI[i].m_nIndex;
		double A1= double(nZ-i)/nZ;
		double A0= double(nZ0-i0)/nZ0;
		double q= pi0*A0/A1;
		if (q<qmin)
			qmin= q;
		vQ.SetAt(j, qmin);

		double p= double(nZ0-i0+1)/(nZ0+1);
		vP.SetAt(j, p);
	}
}

CVector<double> QfromP(CVector<double> vP, bool bSortQ, bool bEstimatePi0)
{
	// Sort P's
	const int N= vP.GetSize();
	CVector<Stat_SortItem_double> vSI;
	vSI.SetSize(N);
	for (int i=0;i<N;i++)
	{
		Stat_SortItem_double si;
		si.m_Value= vP[i];
		si.m_nIndex= i;
		vSI.SetAt(i, si);
	}
	Stat_SortItem_double *pSI= vSI.GetBuffer(N);
	Stat_QuicksortIndexed_double(pSI, N, true);

	/*
	// Override Benjamini-Yuketeli correction, i.e. assume independent or 
	// positively dependent p's which yields the standard, and slightly 
	// less convervative Benjamini-Hochberg correction.

	double benj_yuke= 0.0;
	for (i=0;i<N;i++)
		benj_yuke += 1.0/(i+1);
	*/

	// Estimate pi0
	double pi0= 1.0;
	if (bEstimatePi0 && N>1000)
	{
		printf("Estimating pi0\n");
		int nHigh= 0;
		double fLambda= 0.6f; // Cutoff
		for (int i=0;i<N;i++)
		{
			if (vP[i]>fLambda)
				nHigh++;
		}
		pi0= (double(nHigh)/(1.0f-fLambda))/N;
		if (pi0>1.0)
			pi0= 1.0;
	}
	//else
	//	printf("Skipping pi0 estimation\n");

	//if (GetVerbosity()>1)
	//	printf("Type of multiple-testing correction: False discovery rate (pi0=%.2f).\n", pi0);

	CVector<double> vQ;
	vQ.SetSize(N);
	double fmin= pi0;
	for (int i=N;--i>=0;)
	{
		int j= vSI[i].m_nIndex;
		double f= pi0*vSI[i].m_Value/(double(i+1)/N);

		// Constrained
		if (f<fmin)
			fmin= f;

		if (bSortQ)
			vQ.SetAt(i, fmin);
		else
			vQ.SetAt(j, fmin);
	}

	return vQ;
}

int FindInIndex_alphabetic(const CVector<Table_SortItem_int> &vSI, CString &sValue)
{
	// Assumes vSI is sorted alphabetically in ascending order.
	// Returns position of sValue in vSI.
	int i0= 0;
	int i1= vSI.GetSize();
	while (i0<i1)
	{
		int i_mid= (i0+i1)/2;

		if (vSI[i_mid].m_sKey==sValue)
			return i_mid;

		if (strcmp(vSI[i_mid].m_sKey, sValue)<0)
			i0= i_mid+1;
		else
			i1= i_mid;		
	}

	return -1;
}

CVector<double> GetCorrectedP(const CVector<double> &v_p)
{
	double p_mid= Stat_GetMedian_double(v_p);
	double lambda= 0.5/p_mid;
	printf("Assay inflation factor= %g\n", lambda);

	const int n=v_p.GetSize();
	CVector<double> v_out(n);
	for (int i=0;i<n;i++)
		v_out.SetAt(i, __min(1.0, v_p[i]*lambda));

	return v_out;
}

void FitModel_state0(const double* pA, const double* pG, const CMatrix<double>& mCounts_norm, SegItem* pS, const int s0, const int s1, double& L2)
{
	// fits an allele-independent model

	const int n_samples = mCounts_norm.GetHeight();

	L2 = 0;
	for (int i = s0; i < s1; i++)
	{
		double e = 0;
		for (int j = 0; j < n_samples; j++)
		{
			double z = mCounts_norm.GetAt(j, i);
			e += z;
		}
		double x = e / n_samples;

		for (int j = 0; j < n_samples; j++)
		{
			double residual = x - mCounts_norm.GetAt(j, i);
			L2 += residual * residual;
		}
	}
}

void FitModel_state0_init(double* pA, double *p_state0_cumulator, const CMatrix<double>& mCounts_norm)
{
	const int n_seg = mCounts_norm.GetWidth();
	const int n_samples = mCounts_norm.GetHeight();

	for (int i = 0; i < n_seg; i++)
		pA[i] = 0;
	for (int j = 0; j < n_samples; j++)
	{
		const double* p_counts = mCounts_norm.GetPointer(j, 0);
		for (int i = 0; i < n_seg; i++)
			pA[i] += p_counts[i];
	}
	for (int i = 0; i < n_seg; i++)
		pA[i] /= n_samples;

	for (int i = 0; i < n_seg; i++)
	{
		double L2_cumulator = 0;
		for (int j = 0; j < n_samples; j++)
		{
			const double* p_counts = mCounts_norm.GetPointer(j, 0);
			double res = p_counts[i] - pA[i];
			L2_cumulator += res * res;
		}
		p_state0_cumulator[i] = L2_cumulator;
	}
}

void FitModel_state0_cumulate(const double* pA, const double* pG, const CMatrix<double>& mCounts_norm, const int s, double& L2_cumulator)
{
	const int n_samples = mCounts_norm.GetHeight();
	for (int j = 0; j < n_samples; j++)
	{
		const double* p_counts = mCounts_norm.GetPointer(j, 0);
		double res = p_counts[s] - pA[s];
		L2_cumulator += res * res;
	}
}

void FitModel_state1_init(double* pC, double* pE, double* pX, double* pY, const double* pG, const CMatrix<double>& mCounts_norm, const double s_g, const double ss_g)
{
	// initializes pC, pE, pX, pY given pG, mCounts_norm, s_g and ss_g
	const int n_seg = mCounts_norm.GetWidth();
	const int n_samples = mCounts_norm.GetHeight();

	double a, b, d;
	a = s_g;
	b = ss_g;
	d = n_samples;
	double mul = 1.0 / (a * a - b * d);

	for (int i = 0; i < n_seg; i++)
	{
		pC[i] = 0;
		pE[i] = 0;
	}

	for (int j = 0; j < n_samples; j++)
	{
		const double* p_counts = mCounts_norm.GetPointer(j, 0);
		for (int i = 0; i < n_seg; i++)
		{
			pE[i] += p_counts[i];
			pC[i] += pG[j] * p_counts[i];
		}
	}

	for (int i = 0; i < n_seg; i++)
	{
		pX[i] = (a * pC[i] - b * pE[i]);
		pY[i] = (a * pE[i] - pC[i] * d);
	}
}

void FitModel_state1_cumulate(const double *pC, const double *pE, const double *pX, const double *pY, const double* pG, const CMatrix<double>& mCounts_norm, const int s, const double s_g, const double ss_g, double& L2_cumulator)
{
	// fits an allele-dependent model
	const int n_samples = mCounts_norm.GetHeight();

	// solve ax+by=c, dx+ay=e
	double a, b, d;
	a = s_g;
	b = ss_g;
	d = n_samples;
	double mul = 1.0 / (a * a - b * d);

	for (int j = 0; j < n_samples; j++)
	{
		const double* p_counts = mCounts_norm.GetPointer(j, 0);
		double res = (pX[s] + pG[j] * pY[s])*mul - p_counts[s];
		L2_cumulator += res * res;
	}
}

void FitModel_state1(double* pC, double* pE, double* pX, double* pY, double* pG, const CMatrix<double>& mCounts_norm, SegItem* pS, const int s0, const int s1, const double s_g, const double ss_g, double& L2)
{
	// fits an allele-dependent model
	const int n_samples = mCounts_norm.GetHeight();

	// solve ax+by=c, dx+ay=e
	double a, b, d;
	a = s_g;
	b = ss_g;
	d = n_samples;
	double mul = 1.0 / (a * a - b * d);

	L2 = 0;
	for (int j = 0; j < n_samples; j++)
	{
		const double* p_counts = mCounts_norm.GetPointer(j, 0);
		for (int i = s0; i < s1; i++)
		{
			double res = (pX[i] + pG[j] * pY[i]) * mul - p_counts[i];
			L2 += res * res;
		}
	}
}

void UpHeap(int* pH, int h0, SegItem* pS)
{
	int s = pH[h0];
	h0 = pS[s].m_HeapPos;
	const double c = pS[s].m_cost;

	int h1;
	while (pS[pH[h1 = (h0 >> 1)]].m_cost > c)
	{
		pH[h0] = pH[h1];
		pS[pH[h0]].m_HeapPos = h0;
		h0 = h1;
	}
	pH[h0] = s;
	pS[pH[h0]].m_HeapPos = h0;

}

void CheckHeap(int* pH, int nH, SegItem* pS)
{
#ifdef _DEBUG
	for (int i = 1; i < nH; i++)
	{
		if (pS[pH[i]].m_cost < pS[pH[i >> 1]].m_cost)
		{
			printf("nH= %d\n", nH);
			printf("%d\n", i);
			printf("pH[%d]= %d\n", i, pH[i]);
			printf("pS[%d].m_cost= %g", pH[i], pS[pH[i]].m_cost);
			i = i >> 1;
			printf("%d\n", i);
			printf("pH[%d]= %d\n", i, pH[i]);
			printf("pS[%d].m_cost= %g", pH[i], pS[pH[i]].m_cost);

			ASSERT(false);
		}
	}
#endif
}

bool CApp::DumpCounts(const char *aFilename, CMatrix<double>& mCounts_norm, CString &sChr, int nMinBP, int nMaxBP, int nBinSpacing)
{
	FILE* pf = OpenOutputFile(aFilename);
	if (!pf)
		return false;

	fprintf(pf, "CHR\tBP");
	for (int c = 0; c < mCounts_norm.GetWidth(); c++)
		fprintf(pf, "\tS%d", c + 1);
	fprintf(pf, "\n");

	for (int c=0;c<mCounts_norm.GetWidth();c++)
	{
		fprintf(pf, "%s\t%d", (const char *)sChr, nMinBP+nBinSpacing*c);
		for (int r = 0; r < mCounts_norm.GetHeight(); r++)
			fprintf(pf, "\t%g", mCounts_norm.GetAt(r,c));
		fprintf(pf, "\n");
	}
	CloseFile(pf);
	return true;
}

void UpdateCost(int *pH, int &nH, SegItem *pS, int s1, int s_root, double delta_cost, int state_new)
{
 	double cost_new = delta_cost + pS[s_root].m_cost;
	if (cost_new < pS[s1].m_cost)
	{
		pS[s1].m_cost = cost_new;
		pS[s1].m_PreviousIndex = s_root;
		pS[s1].m_state = state_new;
		if (pS[s1].m_HeapPos < 0)
		{
			nH++;
			pH[nH] = s1;
			pS[s1].m_HeapPos = nH;
		}
		UpHeap(pH, pS[s1].m_HeapPos, pS);
	}
}

void UpdateCost2(int* pH, int& nH, SegItem* pS, int s1, int s_root, double cost_new, int state_new)
{
	// assumes the commented lines are checked before the call, to avoid unnecessary stack ops
	// double cost_new = delta_cost + pS[s_root].m_cost;
	// if (cost_new < pS[s1].m_cost)
	// {
		pS[s1].m_cost = cost_new;
		pS[s1].m_PreviousIndex = s_root;
		pS[s1].m_state = state_new;
		if (pS[s1].m_HeapPos < 0)
		{
			nH++;
			pH[nH] = s1;
			pS[s1].m_HeapPos = nH;
		}
		UpHeap(pH, pS[s1].m_HeapPos, pS);
	// }
}

double FitModel_AutoSelectLambda1(const double* pA, CVector<double>& vG, CMatrix<double>& mCounts_norm,
	const double s_g, const double ss_g, const double lambda_2)
{
	// permute genotypes
	CVector<double> vG_perm = vG;

	const int n_perm = 20;
	CMatrix<double> mF(n_perm, mCounts_norm.GetWidth());

	CVector<double> vC_perm(mCounts_norm.GetWidth());
	CVector<double> vE_perm(mCounts_norm.GetWidth());
	CVector<double> vX_perm(mCounts_norm.GetWidth());
	CVector<double> vY_perm(mCounts_norm.GetWidth());
	double* pC_perm = vC_perm.GetBuffer();
	double* pE_perm = vE_perm.GetBuffer();
	double* pX_perm = vX_perm.GetBuffer();
	double* pY_perm = vY_perm.GetBuffer();

	//double L2_0_tot = 0;
	//double L2_1_tot = 0;
	for (int p = 0; p < n_perm; p++)
	{
		if (true)
		{
			for (int i = 0; i < vG_perm.GetSize(); i++)
			{
				int a = rand_int31() % vG_perm.GetSize();
				int b = rand_int31() % vG_perm.GetSize();
				double g = vG_perm.GetAt(a);
				vG_perm.SetAt(a, vG_perm.GetAt(b));
				vG_perm.SetAt(b, g);
			}
		}
		double* pG_perm = vG_perm.GetBuffer();

		FitModel_state1_init(pC_perm, pE_perm, pX_perm, pY_perm, pG_perm, mCounts_norm, s_g, ss_g);

		//double L2_0_tot_p = 0;
		//double L2_1_tot_p = 0;
		for (int i = 0; i < mCounts_norm.GetWidth(); i++)
		{
			double L2_0 = 0;
			FitModel_state0_cumulate(pA, pG_perm, mCounts_norm, i, L2_0);
			double L2_1 = 0;
			FitModel_state1_cumulate(pC_perm, pE_perm, pX_perm, pY_perm,
				pG_perm, mCounts_norm, i, s_g, ss_g, L2_1);

			//L2_0_tot_p += L2_0;
			//L2_1_tot_p += L2_1;

			double f = L2_0 / L2_1;
			if (IsNaN_double(f))
				f = 1.0;
			mF.SetAt(p, i, f);
		}

		//L2_0_tot += L2_0_tot_p;
		//L2_1_tot += L2_1_tot_p;
	}

	//printf("tot ratio= %g\n", L2_0_tot / L2_1_tot);

	double* pF = mF.GetPointer(0, 0);
	Stat_Quicksort_double(pF, mF.GetWidth() * mF.GetHeight(), true);
	if (false)
	{
		for (double q = 0.95; q <= 1.00001; q += 0.001)
		{
			double f = Stat_GetQuantile_double(pF, mF.GetWidth() * mF.GetHeight(), q);
			printf("%g\t%g\n", q, f);
		}
	}

	return Stat_GetQuantile_double(pF, mF.GetWidth() * mF.GetHeight(), 0.95);
}

bool CApp::ExportIGV(const CString& sBasename, CMatrix<double>& mCounts_norm, const double* pG, const double* pX, const double* pY)
{
	FILE* pf = OpenOutputFile(sBasename + ".igv");
	if (!pf)
		return false;

	const int n_bins = mCounts_norm.GetWidth();
	CVector<int> v_n_00(n_bins);
	CVector<int> v_n_01(n_bins);
	CVector<int> v_n_11(n_bins);
	CVector<double> v_mu_00(n_bins);
	CVector<double> v_mu_01(n_bins);
	CVector<double> v_mu_11(n_bins);
	CVector<double> v_mu_00_norm(n_bins);
	CVector<double> v_mu_01_norm(n_bins);
	CVector<double> v_mu_11_norm(n_bins);
	CVector<double> v_r(n_bins);
	CVector<double> v_p(n_bins);
	CVector<double> v_q(n_bins);

	const int n_samples = mCounts_norm.GetHeight();
	CVector<double> vX(n_samples);

	fprintf(pf, "Chromosome\tStart\tEnd\tFeature\tn00\tn01\tn11\taverage00\taverage01\taverage11\taverage00norm\taverage01norm\taverage11norm\tr\tp\tlog10p\tq\tlog10q\n");
	for (int i = 0; i < n_bins; i++)
	{
		double mu_00 = 0;
		double mu_01 = 0;
		double mu_11 = 0;
		int n_00 = 0;
		int n_01 = 0;
		int n_11 = 0;
		for (int r = 0; r < n_samples; r++)
		{
			if (pG[r] < 0.5)
			{
				mu_00 += mCounts_norm.GetAt(r, i);
				n_00++;
			}
			else if (pG[r] > 1.5)
			{
				mu_11 += mCounts_norm.GetAt(r, i);
				n_11++;
			}
			else
			{
				mu_01 += mCounts_norm.GetAt(r, i);
				n_01++;
			}
			vX.SetAt(r, mCounts_norm.GetAt(r, i));
		}

		double mu = (mu_00 + mu_01 + mu_11) / (n_00 + n_01 + n_11);

		mu_00 /= n_00;
		mu_01 /= n_01;
		mu_11 /= n_11;
		v_mu_00.SetAt(i, mu_00);
		v_mu_01.SetAt(i, mu_01);
		v_mu_11.SetAt(i, mu_11);
		if (mu != 0.0)
		{
			v_mu_00_norm.SetAt(i, mu_00 / mu);
			v_mu_01_norm.SetAt(i, mu_01 / mu);
			v_mu_11_norm.SetAt(i, mu_11 / mu);
		}
		else
		{
			v_mu_00_norm.SetAt(i, 1.0);
			v_mu_01_norm.SetAt(i, 1.0);
			v_mu_11_norm.SetAt(i, 1.0);
		}

		double corr = Stat_GetCorrelation_double(vX.GetBuffer(), pG, vX.GetSize());
		
		double t = corr * sqrt(vX.GetSize() - 2) / sqrt(1.0 - corr * corr);
		double p = 2*(1.0-Stat_Tcdf(__abs(t), vX.GetSize()-2));

		v_n_00.SetAt(i, n_00);
		v_n_01.SetAt(i, n_01);
		v_n_11.SetAt(i, n_11);
		v_r.SetAt(i, corr);
		v_p.SetAt(i, p);
	}
	v_q = QfromP(v_p, false, false);

	for (int i = 0; i < n_bins; i++)
	{
		int nStart, nEnd;
		if (true)
		{
			nStart = m_nMinBP + m_nBinSpacing * i;
			nEnd = nStart + m_nBinSpacing;
		}
		else
		{
			nStart = m_nMinBP + ((2 * i + 1) * m_nBinSpacing) / 2;
			nEnd = nStart + 1;
		}

		fprintf(pf, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
			(const char*)m_sCHR, nStart, nEnd, "density", 
			v_n_00[i], v_n_01[i], v_n_11[i], 
			v_mu_00[i], v_mu_01[i], v_mu_11[i], 
			v_mu_00_norm[i], v_mu_01_norm[i], v_mu_11_norm[i], 
			v_r[i], 
			v_p[i], -log10(v_p[i]), v_q[i], -log10(v_q[i]));
	}

	CloseFile(pf);
	return true;
}

bool CApp::ExportPI0(const CString &sFilename, CVector<double>& vLambda1, CVector<double>& vLambda2, CMatrix<double>& mPI0, CMatrix<double>& mH1proportion)
{
	if (vLambda1.GetSize() != mPI0.GetHeight())
	{
		ReportError("Internal error: vLambda1.GetSize() != mPI0.GetHeight()", 0);
		return false;
	}
	if (vLambda2.GetSize() != mPI0.GetWidth())
	{
		ReportError("Internal error: vLambda2.GetSize() != mPI0.GetWidth()", 0);
		return false;
	}

	FILE* pf = OpenOutputFile(sFilename);
	if (!pf)
		return false;

	for (int i1 = -1; i1 < vLambda1.GetSize(); i1++)
	{
		if (i1 < 0)
		{
			// header row
			fprintf(pf, "%g", 0.0); // dummy value
			for (int i2 = 0; i2 < vLambda2.GetSize(); i2++)
				fprintf(pf, "\t%g", vLambda2[i2]);
		}
		else
		{
			// body row

			// H1 proportion with PI0 within parenthesis
			fprintf(pf, "%g", vLambda1[i1]);
			for (int i2 = 0; i2 < vLambda2.GetSize(); i2++)
			{
				fprintf(pf, "\t%.2g%% ", 100*mH1proportion.GetAt(i1, i2));
				if (IsN_double(mPI0.GetAt(i1, i2)))
					fprintf(pf, "(%.3g)", mPI0.GetAt(i1, i2));
				else
					fprintf(pf, "(%s)", (const char *)m_sNaN);
			}
		}
		fprintf(pf, "\n");
	}

	CloseFile(pf);
	return true;
}

bool CApp::Optimize(const double lambda1, const double lambda2, CVector<double>& vG, CMatrix<int>& mCounts, int nMinBP, int nMaxBP, int nBinSpacing, const CString &sBasename, int &nBPstate1)
{
	if (mCounts.GetWidth() > 100000)
	{
		ReportError("Too many dynamic programming cells, increase -BinSpacing or change -minbp, -maxbp", 0);
		return false;
	}

	// Segmentation voxels
	int n_items = mCounts.GetWidth();
	CVector<SegItem> vS(n_items + 1); // +1 for downheap anchor with infinite cost
	SegItem* pS = vS.GetBuffer();
	for (int i = 0; i < n_items; i++)
	{
		pS[i].m_cost = Stat_GetPosInfty_double();
		pS[i].m_HeapPos = -1;
		pS[i].m_PreviousIndex = -1;
		pS[i].m_state = -1;
		pS[i].m_leftbp = nMinBP + nBinSpacing * i;
	}
	pS[n_items].m_cost = 1.0e100;
	pS[n_items].m_HeapPos = -1;
	pS[n_items].m_leftbp = -1;
	pS[n_items].m_PreviousIndex = -1;
	pS[n_items].m_state = -1;

	// Initialize count data and buffers
	CMatrix<double> mCounts_norm;
	const int dBP = nMaxBP - nMinBP;
	mCounts_norm.ReInit(mCounts.GetHeight(), mCounts.GetWidth());

	// code for non-cumulative initiation and normalization
	int dc = m_nBinSize / m_nBinSpacing / 2;
	for (int r = 0; r < mCounts.GetHeight(); r++)
	{
		int s = 0;
		for (int c = 0; c < mCounts.GetWidth(); c++)
			s += mCounts.GetAt(r, c);
		for (int c = 0; c < mCounts.GetWidth(); c++)
		{
			int sc = 0;
			int n = 0;
			for (int i = c - dc; i <= c + dc; i++)
			{
				if (i >= 0 && i < mCounts.GetWidth())
				{
					sc += mCounts.GetAt(r, i);
					n++;
				}
			}
			mCounts_norm.SetAt(r, c, double(sc) / n / s * mCounts.GetWidth()); //  1000); // TODO: Change 1000 to Counts.GetWidth(), retune lambda1 and lambda2
		}
	}

	// initiate model fitting buffers for state 0 
	CVector<double> vA(mCounts_norm.GetWidth());
	double* pA = vA.GetBuffer();
	CVector<double> v_state0_cumulator(mCounts_norm.GetWidth());
	double* p_state0_cumulator= v_state0_cumulator.GetBuffer();
	FitModel_state0_init(pA, p_state0_cumulator, mCounts_norm);

	// initiate model fitting buffers for state 1
	CVector<double> vC(mCounts_norm.GetWidth());
	CVector<double> vE(mCounts_norm.GetWidth());
	CVector<double> vX(mCounts_norm.GetWidth());
	CVector<double> vY(mCounts_norm.GetWidth());
	double* pC = vC.GetBuffer();
	double* pE = vE.GetBuffer();
	double* pX = vX.GetBuffer();
	double* pY = vY.GetBuffer();

	double* pG = vG.GetBuffer(); // genotype buffer (input)
	double s_g = 0;
	double ss_g = 0;
	for (int i = 0; i < vG.GetSize(); i++)
	{
		s_g += pG[i];
		ss_g += pG[i] * pG[i];
	}
	FitModel_state1_init(pC, pE, pX, pY, pG, mCounts_norm, s_g, ss_g);
	CVector<double> v_state1_cumulator(mCounts_norm.GetWidth());
	double* p_state1_cumulator = v_state1_cumulator.GetBuffer();
	for (int i = 0; i < mCounts_norm.GetWidth(); i++)
	{
		double L2_cumulator = 0;
		FitModel_state1_cumulate(pC, pE, pX, pY, pG, mCounts_norm, i, s_g, ss_g, L2_cumulator);
		p_state1_cumulator[i] = L2_cumulator;
	}

	double Lambda1_adjusted = lambda1;
	double Lambda2_adjusted = lambda2 * mCounts_norm.GetHeight();	// make effect of lambda2 independent of sample size

	//
	// The optimization is dynamic programming, implemented one using Dijkstraa's shortest-path 
	// algorithm with heap-sorted queue and backpointers (=computationally efficient)
	//
	CVector<int> vH(2 * (n_items + 1));
	vH.Fill(n_items); // point to downheap anchor element
	int* pH = vH.GetBuffer();
	pH[0] = 0;
	pS[0].m_cost = 0.0;
	pS[0].m_PreviousIndex = -1;

	pH[1] = 0;
	pS[0].m_HeapPos = 1;
	int nH = 1;

	while (nH > 0)
	{
		int s_root = pH[1];
		pS[s_root].m_HeapPos = -1;

		pH[1] = pH[nH];
		pH[nH] = n_items; // point to downheap anchor element
		pS[pH[1]].m_HeapPos = 1;

		nH--;

		int h0 = 1;
		int h1_left = h0 << 1;
		int h1_right = h1_left + 1;
		while ((pS[pH[h1_left]].m_cost < pS[pH[h0]].m_cost) || (pS[pH[h1_right]].m_cost < pS[pH[h0]].m_cost))
		{
			if (pS[pH[h1_left]].m_cost <= pS[pH[h1_right]].m_cost)
			{
				// swap left
				int s = pH[h0];
				pH[h0] = pH[h1_left];
				pH[h1_left] = s;

				pS[pH[h0]].m_HeapPos = h0;
				pS[pH[h1_left]].m_HeapPos = h1_left;

				h0 = h1_left;
			}
			else
			{
				// swap right
				int s = pH[h0];
				pH[h0] = pH[h1_right];
				pH[h1_right] = s;

				pS[pH[h0]].m_HeapPos = h0;
				pS[pH[h1_right]].m_HeapPos = h1_right;

				h0 = h1_right;
			}
			h1_left = h0 << 1;
			h1_right = h1_left + 1;
		}

		// expand rightward paths from root element
		int state_new;
		if (pS[s_root].m_state == 0)
		{
			// last state 0 --> fit state 1
			state_new = 1;
			double L2_cumulator_state1 = 0;
			for (int s1 = s_root + 1; s1 < n_items; s1++)
			{
				int seglen = pS[s1].m_leftbp - pS[s_root].m_leftbp;
				if (state_new == 1 && seglen > m_nMaxSegmentLength_state1) // limit max seg length in state1 
					break;
				L2_cumulator_state1 += p_state1_cumulator[s1 - 1];

				double delta_cost = Lambda1_adjusted * L2_cumulator_state1 + Lambda2_adjusted;
				double cost_new = delta_cost + pS[s_root].m_cost;
				if (cost_new < pS[s1].m_cost)
					UpdateCost2(pH, nH, pS, s1, s_root, cost_new, state_new);
			}
		}
		else if (pS[s_root].m_state == 1)
		{
			// last state 1 --> fit state 0
			state_new = 0;
			double L2_cumulator_state0 = 0;
			for (int s1 = s_root + 1; s1 < n_items; s1++)
			{
				L2_cumulator_state0 += p_state0_cumulator[s1 - 1];

				double delta_cost = L2_cumulator_state0 + Lambda2_adjusted;
				double cost_new = delta_cost + pS[s_root].m_cost;
				if (cost_new < pS[s1].m_cost)
					UpdateCost2(pH, nH, pS, s1, s_root, cost_new, state_new);
			}
		}
		else
		{
			// last state undefined --> fit best state
			double L2_cumulator_state0 = 0;
			double L2_cumulator_state1 = 0;
			for (int s1 = s_root + 1; s1 < n_items; s1++)
			{
				int seglen = pS[s1].m_leftbp - pS[s_root].m_leftbp;

				FitModel_state0_cumulate(pA, pG, mCounts_norm, s1-1, L2_cumulator_state0);
				if (seglen < m_nMaxSegmentLength_state1)
					FitModel_state1_cumulate(pC, pE, pX, pY, pG, mCounts_norm, s1 - 1, s_g, ss_g, L2_cumulator_state1);
				else
					L2_cumulator_state1 = 1.0e100;

				if (L2_cumulator_state0 < Lambda1_adjusted*L2_cumulator_state1)
				{
					state_new = 0;
					UpdateCost(pH, nH, pS, s1, s_root, L2_cumulator_state0 + Lambda2_adjusted, state_new);
				}
				else
				{
					state_new = 1;
					UpdateCost(pH, nH, pS, s1, s_root, Lambda1_adjusted*L2_cumulator_state0 + Lambda2_adjusted, state_new);
				}
			}
		}
	}

	//
	// dynamic programming back-tracking step
	int s = n_items - 1;
	int n_seg = 0;
	int n_right = __min(m_nMaxBP, pS[s].m_leftbp + nBinSpacing);
	int state = pS[s].m_state;
	int s_prev = s;
	s = pS[s].m_PreviousIndex;

	//
	// export solution
	FILE* pf_seg = 0;
	if (m_bExportSEG && !sBasename.IsEmpty())
	{
		pf_seg = OpenOutputFile(sBasename + ".seg");
		if (!pf_seg)
			return false;
	}

	nBPstate1 = 0;
	if (!sBasename.IsEmpty() && (GetVerbosity() > 2))
		printf("Chromosome\tStart\tEnd\tLength\tState\tL2_ratio\tl1_threshold\tl2_threshold\n");
	if (pf_seg!=0)
		fprintf(pf_seg, "ID\tchrom\tloc.start\tloc.end\tnum.mark\tL2norm_ratio\tLambda2_threshold\tseg.mean\n");
	while (s >= 0)
	{
		n_seg++;
		int n_left = pS[s].m_leftbp;

		double cost_0, cost_1;
		FitModel_state0(pA, pG, mCounts_norm, pS, s, s_prev, cost_0);
		FitModel_state1(pC, pE, pX, pY, pG, mCounts_norm, pS, s, s_prev, s_g, ss_g, cost_1);
		double L2norm_ratio = cost_0 / cost_1;
		double cost_ratio = cost_0 / (Lambda1_adjusted * cost_1 + 2 * Lambda2_adjusted);
		double Lambda1_threshold = (cost_0 - 2 * Lambda2_adjusted) / cost_1;
		double Lambda2_adjusted_threshold = (cost_0 - Lambda1_adjusted * cost_1) / 2;
		double Lambda2_threshold = lambda2 / Lambda2_adjusted * Lambda2_adjusted_threshold;
		CString sID = "caQTLseg"; //  ::Format("%s l1=%g l2=%g", m_sRSID, m_Lambda1, m_Lambda2);

		if (state==1)
		{
			nBPstate1 += (n_right - n_left);
			if (pf_seg != 0)
			{
				fprintf(pf_seg, "%s\t%s\t%d\t%d\t%d\t%g\t%g\t%d\n", (const char*)sID, (const char*)m_sCHR, n_left, n_right, (n_right - n_left) / m_nBinSpacing, L2norm_ratio, Lambda2_threshold, 10 * state);
				if (GetVerbosity()>2)
					printf("%s\t%d\t%d\t%d\t%d\t%g\t%g\t%g\n", (const char*)m_sCHR, n_left, n_right, n_right - n_left, state, L2norm_ratio, Lambda1_threshold, Lambda2_threshold);
			}
		}
		n_right = n_left;
		state = pS[s].m_state;
		s_prev = s;
		s = pS[s].m_PreviousIndex;
	}

	if (pf_seg != 0)
		CloseFile(pf_seg);

	bool bOk = true;
	if (m_bExportIGV && !sBasename.IsEmpty())
		bOk= ExportIGV(sBasename, mCounts_norm, pG, pX, pY);
	
	return bOk;
}

/******************************************************************************/
// Main application

bool CApp::ParseGenotypes(const char *aFilename, CVector<double> &vG, CString &sA1, CString &sA2, CVector<CString> &vSampleID, int& n00, int& n01, int& n11, int& nNA)
{
	if (!m_idx_genotypes.IsRectangular())
	{
		ReportError("Genotype table must be rectangular", aFilename);
		return false;
	}
	if (m_idx_genotypes.GetMinColCount() < 3)
	{
		ReportError("Genotype table must have at least three columns", 0);
		return false;
	}

	int cRSID = Table_FindCol(m_idx_genotypes, "RSID");
	if (cRSID!=0)
	{
		ReportError("First column must be named 'RSID'", aFilename);
		return false;
	}
	int cA1= Table_FindCol(m_idx_genotypes, "A1");
	if (cA1!=1)
	{
		ReportError("Second column must be named 'A1'", aFilename);
		return false;
	}
	int cA2 = Table_FindCol(m_idx_genotypes, "A2");
	if (cA2!=2)
	{
		ReportError("Third column must be named 'A2'", aFilename);
		return false;
	}

	// find variant row
	int rSNP = Table_FindRow(m_idx_genotypes, m_sRSID);
	if (rSNP < 0)
	{
		ReportError(::Format("Row '%s' not found", (const char *)m_sRSID), aFilename);
		return false;
	}

	sA1 = Table_GetAt(m_idx_genotypes, rSNP, cA1);
	sA2 = Table_GetAt(m_idx_genotypes, rSNP, cA2);

	n00 = 0;
	n01 = 0;
	n11 = 0;
	nNA = 0;

	double g = 0;
	const int C = m_idx_genotypes.GetMinColCount();
	for (int c = 3; c < C; c++)
	{
		CString s_g_orig = Table_GetAt(m_idx_genotypes, rSNP, c);

		CString s_g;
		for (int i = 0; i < s_g_orig.GetLength(); i++)
			if (s_g_orig != '/') // remove any between-allele slashes
				s_g += s_g_orig[i];

		bool bNA = false;
		if (s_g == sA1 + sA1)
		{
			g = 0.0;
			n00++;
		}
		else if (s_g == sA1 + sA2 || s_g == sA2 + sA1)
		{
			g = 1.0;
			n01++;
		}
		else if (s_g == sA2 + sA2)
		{
			g = 2.0;
			n11++;
		}
		else if (s_g == "NA" || s_g == "na")
		{
			bool bNA = true;
			nNA++;
		}
		else
		{
			ReportError("Illegal genotypes", s_g);
			return false;
		}

		if (!bNA)
		{
			CString sSampleID = Table_GetAt(m_idx_genotypes, 0, c);
			if (vSampleID.Find(sSampleID) >= 0)
			{
				ReportError(::Format("Sample name '%s' found more than once", (const char*)sSampleID), aFilename);
				return false;
			}
			vG.Add(g);
			vSampleID.Add(sSampleID);
		}
	}

	if (vG.IsEmpty())
	{
		ReportError("No samples found in", aFilename);
		return false;
	}

	return true;
}

bool CApp::ParseGenotypes_multi(const char* aFilename, CVector<CString>& vSampleID, CMatrix<double>& mG, CVector<CString> &vRSID, CVector<CString>& vCHR, CVector<int>& vBP, CVector<CString>& vA1, CVector<CString>& vA2, CVector<int>& vn00, CVector<int>& vn01, CVector<int>& vn11, CVector<int>& vnNA)
{
	const int nAnnotRows = 1;
	const int nAnnotCols = 5;
	if (!m_idx_genotypes.IsRectangular())
	{
		ReportError("Genotype table must be rectangular", aFilename);
		return false;
	}
	if (m_idx_genotypes.GetMinColCount() < nAnnotCols)
	{
		ReportError(::Format("Genotype table must have at least %d columns", nAnnotCols), 0);
		return false;
	}

	//
	int cRSID = Table_FindCol(m_idx_genotypes, "RSID");
	if (cRSID != 0)
	{
		ReportError("First column must be named 'RSID'", aFilename);
		return false;
	}
	int cCHR = Table_FindCol(m_idx_genotypes, "CHR");
	if (cCHR != 1)
	{
		ReportError("Second column must be named 'CHR'", aFilename);
		return false;
	}
	int cBP = Table_FindCol(m_idx_genotypes, "BP");
	if (cBP != 2)
	{
		ReportError("Third column must be named 'BP'", aFilename);
		return false;
	}
	int cA1 = Table_FindCol(m_idx_genotypes, "A1");
	if (cA1 != 3)
	{
		ReportError("Fourth column must be named 'A1'", aFilename);
		return false;
	}
	int cA2 = Table_FindCol(m_idx_genotypes, "A2");
	if (cA2 != 4)
	{
		ReportError("Fifth column must be named 'A2'", aFilename);
		return false;
	}

	const int C = m_idx_genotypes.GetMinColCount();
	for (int c = nAnnotCols; c < C; c++)
	{
		CString sSampleID = Table_GetAt(m_idx_genotypes, 0, c);
		if (sSampleID.IsEmpty())
		{
			ReportError("Sample ID cannot be empty", aFilename);
			return false;
		}
		if (vSampleID.Find(sSampleID) >= 0)
		{
			ReportError(::Format("Sample ID '%s' found more than once", (const char*)sSampleID), aFilename);
			return false;
		}
		vSampleID.Add(sSampleID);
	}

	mG.ReInit(m_idx_genotypes.GetRowCount() - nAnnotRows, m_idx_genotypes.GetMinColCount() - nAnnotCols);

	//
	printf("%d\t%d\n", m_idx_genotypes.GetRowCount(), m_idx_genotypes.GetMinColCount());

	int nr = m_idx_genotypes.GetRowCount() - nAnnotRows;
	if (nr <= 0)
	{
		ReportError("No genotypes to parse", 0);
		return false;
	}

	vRSID.SetSize(nr);
	vCHR.SetSize(nr);
	vBP.SetSize(nr);
	vA1.SetSize(nr);
	vA2.SetSize(nr);
	vn00.SetSize(nr);
	vn01.SetSize(nr);
	vn11.SetSize(nr);
	vnNA.SetSize(nr);

	int ir = 0;
	for (int r = nAnnotRows; r < m_idx_genotypes.GetRowCount(); r++)
	{
		CString sRSID = Table_GetAt(m_idx_genotypes, r, cRSID);
		vRSID.SetAt(ir, sRSID);
		CString sCHR = Table_GetAt(m_idx_genotypes, r, cCHR);
		vCHR.SetAt(ir, sCHR);
		CString sBP = Table_GetAt(m_idx_genotypes, r, cBP);
		if (Table_CheckFloatFormat_complete(sBP) == 0)
		{
			ReportError("Illegal basepair format", sBP);
			return false;
		}
		int nBP = atoi(sBP);
		vBP.SetAt(ir, nBP);
		CString sA1 = Table_GetAt(m_idx_genotypes, r, cA1);
		vA1.SetAt(ir, sA1);
		CString sA2 = Table_GetAt(m_idx_genotypes, r, cA2);
		vA2.SetAt(ir, sA2);

		int n00 = 0;
		int n01 = 0;
		int n11 = 0;
		int nNA = 0;
		for (int c = nAnnotCols; c < C; c++)
		{
			const char* p = m_idx_genotypes.GetFirstPointer(r, c);
			double g = 0;
			bool bNA = false;

			if ((sA1.GetLength()==1 && sA2.GetLength()==1))
			{
				if (p[0] == sA1[0] && p[1] == sA1[0])
				{
					// A1A1
					g = 0.0;
					n00++;
				}
				else if ((p[0] == sA1[0] && p[1] == sA2[0]) || (p[0] == sA2[0] && p[1] == sA1[0]))
				{
					// A1A2 or A2A1
					g = 1.0;
					n01++;
				}
				else if (p[0] == sA2[0] && p[1] == sA2[0])
				{
					// A2A2
					g = 2.0;
					n11++;
				}
				else if ((p[0] == 'N' && p[1] == 'A') || (p[0] == 'n' && p[1] == 'a'))
				{
					// NA
					g = -1.0;
					nNA++;
				}
				else
				{
					// unrecognized
					g = -1.0;
					nNA++;
				}
			}
			else
			{
				CString s_g = Table_GetAt(m_idx_genotypes, r, c);
				//CString s_g;
				//for (int i = 0; i < s_g_orig.GetLength(); i++)
				//	if (s_g_orig != '/') // remove any between-allele slashes
				//		s_g += s_g_orig[i];

				if (s_g == sA1 + sA1)
				{
					g = 0.0;
					n00++;
				}
				else if (s_g == sA1 + sA2 || s_g == sA2 + sA1)
				{
					g = 1.0;
					n01++;
				}
				else if (s_g == sA2 + sA2)
				{
					g = 2.0;
					n11++;
				}
				else if (s_g == "NA" || s_g == "na")
				{
					g = -1.0;
					nNA++;
				}
				else
				{
					if (GetVerbosity() > 0)
						ReportWarning("Illegal genotype detected, using NA instead", ::Format("%s [%s '%s' '%s']", (const char*)s_g, (const char*)sRSID, (const char*)sA1, (const char*)sA2));
					g = -1.0;
					nNA++;
					// return false;
				}
			}
			mG.SetAt(r - nAnnotRows, c - nAnnotCols, g); 
		}
		vn00.SetAt(ir, n00);
		vn01.SetAt(ir, n01);
		vn11.SetAt(ir, n11);
		vnNA.SetAt(ir, nNA);
		ir++;
	}
	printf("No. of variants: %d\n", vRSID.GetSize());
	printf("No. of samples: %d\n", vSampleID.GetSize());

	if (mG.IsEmpty())
	{
		ReportError("No genotype data", aFilename);
		return false;
	}

	return true;
}

bool CApp::LoadCutSiteFiles(const char *aFolderName, const CVector<CString> &vSampleID, CMatrix<int> &mCounts, const CString &sCHR, int nMinBP, int nMaxBP, int nBinSpacing)
{
	if (GetVerbosity()>1)
		printf("Loading cut site files...\n");
	for (int i = 0; i < vSampleID.GetSize(); i++)
	{
		CVector<CString> vFiles;
		g_FileSystem.ListFiles(aFolderName, "*" + vSampleID[i] + "*", vFiles);
		if (vFiles.GetSize() == 0)
		{
			ReportError(::Format("No file found for sample '%s'", (const char*)vSampleID[i]), aFolderName);
			return false;
		}
		else if (vFiles.GetSize() > 1)
		{
			ReportError(::Format("More than one file found for sample '%s'", (const char*)vSampleID[i]), aFolderName);
			return false;
		}
		else
		{
			// unique file found, now load it
			// printf("Loading sample: %s\n", (const char *)vFiles[0]);
			int nCuts = 0;

			const int nBufSize = 100000;
			char linebuf[nBufSize];
			CLineReader lr;

			CVector<const char*> vItems(1000);
			int nItems;

			if (!lr.Open(vFiles[0]))
				return false;

			while (lr.GetLine(&linebuf[0], nBufSize, vItems, nItems))
			{
				if (nItems > 1)
				{
					if (vItems[0] == sCHR)
					{
						for (int c = 1; c < nItems; c++)
						{
							if (Table_CheckIntFormat_complete(vItems[c]) != 0)
							{
								int bp = atoi(vItems[c]);
								if (bp >= nMinBP && bp <= nMaxBP)
								{
									int col = (bp - nMinBP) / nBinSpacing;
									mCounts.SetAt(i, col, mCounts.GetAt(i, col)+1);
									// printf("%s\t%d\t%d\n", (const char *)sCHR, bp, col);
									nCuts++;
								}
							}
						}
					}
				}
			}
			lr.Close();
		}
	}
	return true;
}

bool CApp::BED2BIN(const char *aFilenameBED, const char *aFilenameBIN, const CVector<CString> &vsCHRdictionary, const bool bOverwriteCache)
{

	if (!g_FileSystem.FileExists(aFilenameBED))
	{
		ReportError("File does not exist", aFilenameBED);
		return false;
	}

#ifdef _WIN32
	FILETIME tCreateBED, tWriteBED;
	if (!g_FileSystem.GetFileTime(aFilenameBED, &tCreateBED, NULL, &tWriteBED))
#else
	time_t atime, mtime;
	if (!g_FileSystem.GetFileTime(aFilenameBED, atime, mtime))
#endif		
	{
		ReportError("Unable to get file time", aFilenameBED);
		return false;
	}

	// generate cached .bin files if needed
	bool bOk = false;
	if (!bOverwriteCache &&
#ifdef _WIN32		
		g_FileSystem.HasFileTime(aFilenameBIN, &tCreateBED, NULL, &tWriteBED))
#else
		g_FileSystem.HasFileTime(aFilenameBIN, NULL, &mtime))
#endif		
	{
		// cached .bin file exists and is up to date.
		bOk = true;
	}
	else
	{
		// convert
		//if (GetVerbosity()>1)
		//	printf("Creating binarized cache: %s\n", aFilenameBIN);

		const int nBufSize = 100000;
		char linebuf[nBufSize];
		CLineReader lr;
		CVector<const char*> vItems(1000);
		int nItems;
		if (!lr.Open(aFilenameBED))
			return false;

		FILE* pf_bin = OpenOutputFile(aFilenameBIN);
		if (!pf_bin)
		{
			ReportError("Cannot open output file", aFilenameBIN);
			return false;
		}

		CString sCHR_last;
		int line = 0;
		while (lr.GetLine(&linebuf[0], nBufSize, vItems, nItems))
		{
			if (line == 0)
			{
				if (nItems < 1)
				{
					ReportError("Unrecognized .bed format, first line is empty", aFilenameBED);
					return false;
				}
				if (strcmp(vItems[0], "CHR")!=0)
				{
					ReportError("Unrecognized .bed format, 'CHR' expected as first element", aFilenameBED);
					return false;
				}
			}
			else
			{
				if (nItems > 1)
				{
					if (vItems[0] != sCHR_last)
					{
						sCHR_last = vItems[0];
						int32_t i = 0;
						for (; i < vsCHRdictionary.GetSize(); i++)
							if (vItems[0] == vsCHRdictionary[i])
								break;
						if (i == vsCHRdictionary.GetSize())
						{
							ReportError(::Format("Unrecognized chromosome name '%s'", (const char*)vItems[0]), aFilenameBED);
							return false;
						}
						i = -i - 1;
						if (fwrite(&i, sizeof(i), 1, pf_bin) != 1)
						{
							ReportError("Cannot write output file", aFilenameBIN);
							return false;
						}
						// printf("Converting chromosome %s...\n", (const char *)sCHR_last);
					}

					for (int c = 1; c < nItems; c++)
					{
						if (Table_CheckIntFormat_complete(vItems[c]) != 0)
						{
							int32_t bp = atoi(vItems[c]);
							if (bp < 0)
							{
								ReportError(::Format("Base coordinate cannot be negative '%s'", (const char*)vItems[c]), aFilenameBED);
								return false;
							}
							if (fwrite(&bp, sizeof(bp), 1, pf_bin) != 1)
							{
								ReportError("Cannot write output file", aFilenameBIN);
								return false;
							}
						}
						else
						{
							ReportError(::Format("Unrecognized integer format '%s'", (const char*)vItems[c]), aFilenameBED);
							return false;
						}
					}
				}
			}
			line++;
		}
		CloseFile(pf_bin);

		// make sure BIN has same time stamp as BED
#ifdef _WIN32
		if (!g_FileSystem.SetFileTime(aFilenameBIN, &tCreateBED, NULL, &tWriteBED))
#else
		if (!g_FileSystem.SetFileTime(aFilenameBIN, atime, mtime))
#endif
		{
			ReportError("Cannot set file time", aFilenameBIN);
			return false;
		}

		bOk = true;
	}

	return bOk;
}

bool CheckOrder(const Table_SortItem_int *pSI, const int nSI)
{
	for (int i = 1; i < nSI; i++)
	{
		if (pSI[i - 1].m_sKey > pSI[i].m_sKey)
		{
			ReportError("Internal error: chromosomes not sorted correctly", 0);
			return false;
		}
	}
	return true;
}

bool CApp::LoadCutSiteFiles_multi(const char* aFolderName, const CVector<CString>& vSampleID, const CVector<CString>& vSNP_CHR, CVector<int>& vSNP_BP, const int& ww0, const int& ww1, CMatrix<double> &mData_out, CMatrix<double> &mG)
{
	if (vSNP_CHR.GetSize() != vSNP_BP.GetSize())
	{
		ReportError("Internal error", "vSNP_CHR.GetSize() != vSNP_BP.GetSize()");
		return false;
	}
	if (ww1 < ww0)
	{
		ReportError("Normalization window must be wider than effect windows, check -ww0 and -ww1 switches", 0);
		return false;
	}

	const int nSNP = vSNP_CHR.GetSize();

	// build SNP lookup index

	CVector<CString> vCHR_name;
	CVector<int> vCHR_i0;
	CVector<int> vCHR_i1;

	CVector<Table_SortItem_int> vSI;
	vSI.SetSize(nSNP);
	Table_SortItem_int* pSI = vSI.GetBuffer();
	for (int i = 0; i < nSNP; i++)
	{
		Table_SortItem_int si;
		si.m_nIndex = i;
		si.m_sKey = vSNP_CHR[i];
		si.m_Value = vSNP_BP[i];
		pSI[i] = si;
	}
	Table_SortItems_int(pSI, nSNP, true, false); // sort by chromosome
	if (!CheckOrder(pSI, nSNP))
		return false;

	int i0 = 0;
	while (i0 < nSNP)
	{
		vCHR_name.Add(pSI[i0].m_sKey);
		vCHR_i0.Add(i0);
		// printf("Adding chr '%s'\n", (const char*)pSI[i0].m_sKey);
		// printf("i0= %d\n", i0);
		int i1 = i0+1;
		while (i1 < nSNP && pSI[i0].m_sKey == pSI[i1].m_sKey)
			i1++;
		// if (i1 < nSNP)
		//	printf("Breaking because of pSI[%d]= '%s'\n", i1, pSI[i1].m_sKey);
		vCHR_i1.Add(i1);
		// printf("i1= %d\n", i1);
		Table_SortItems_int(&pSI[i0], i1 - i0, true, true); // sort by position within chromosome
		if (!CheckOrder(pSI, nSNP))
			return false;
		i0 = i1;
	}

	// count buffers
	CVector<int> vCuts_ww0(nSNP);
	CVector<int> vCuts_ww1(nSNP);
	int *pCuts_ww0 = vCuts_ww0.GetBuffer();
	int *pCuts_ww1 = vCuts_ww1.GetBuffer();
	
	mData_out.ReInit(vSNP_CHR.GetSize(), vSampleID.GetSize());
	const int* pSNP_BP = vSNP_BP.GetBuffer();

	CVector<CString> vsCHRdictionary; // NOTE: if this dictionary is ever changed, need to overwrite any cached .bin files
	vsCHRdictionary.Add("1");
	vsCHRdictionary.Add("2");
	vsCHRdictionary.Add("3");
	vsCHRdictionary.Add("4");
	vsCHRdictionary.Add("5");
	vsCHRdictionary.Add("6");
	vsCHRdictionary.Add("7");
	vsCHRdictionary.Add("8");
	vsCHRdictionary.Add("9");
	vsCHRdictionary.Add("10");
	vsCHRdictionary.Add("11");
	vsCHRdictionary.Add("12");
	vsCHRdictionary.Add("13");
	vsCHRdictionary.Add("14");
	vsCHRdictionary.Add("15");
	vsCHRdictionary.Add("16");
	vsCHRdictionary.Add("17");
	vsCHRdictionary.Add("18");
	vsCHRdictionary.Add("19");
	vsCHRdictionary.Add("20");
	vsCHRdictionary.Add("21");
	vsCHRdictionary.Add("22");
	vsCHRdictionary.Add("X");
	vsCHRdictionary.Add("Y");

	printf("Loading cut site files...\n");
	for (int s = 0;s < vSampleID.GetSize(); s++)
	{
		CVector<CString> vFiles;
		g_FileSystem.ListFiles(aFolderName, "*" + vSampleID[s] + "*.bed", vFiles);
		if (vFiles.GetSize() == 0)
		{
			ReportError(::Format("No file found for sample '%s'", (const char*)vSampleID[s]), aFolderName);
			return false;
		}
		else if (vFiles.GetSize() > 1)
		{
			ReportError(::Format("More than one file found for sample '%s'", (const char*)vSampleID[s]), aFolderName);
			return false;
		}
		else
		{
			// convert to .bin
			CString sFilenameBIN = g_FileSystem.GetNameWithoutExtension(vFiles[0]) + ".bin";
			if (!BED2BIN(vFiles[0], sFilenameBIN, vsCHRdictionary, m_bMultiOverwriteCache))
				return false;

			CLoadBuf *pLB= g_FileSystem.Load(sFilenameBIN);
			if (!pLB)
			{
				ReportError("Cannot open .bin file", sFilenameBIN);
				return false;
			}

			// reset density counters
			for (int j = 0; j < nSNP; j++)
			{
				pCuts_ww0[j] = 0;
				pCuts_ww1[j] = 0;
			}
			int nCHR_i0 = -1;
			int nCHR_i1 = -1;
			int last_i0 = -1;

			const int32_t* p = (const int32_t *)pLB->m_aLoadBuf;
			if (*p >= 0)
			{
				ReportError("First int32 must be a chromosome index", sFilenameBIN); // chromosome indices are stored as negative basepair int32 values
				return false;
			}
			int np = pLB->m_BufSize / sizeof(int32_t);

			printf("Processing %d of %d: %s [%g Mb]\n", s+1, vSampleID.GetSize(), (const char *)sFilenameBIN, double(pLB->m_BufSize)/1024/1024);
			while (np > 0)
			{
				int32_t bp = *p++;
				if (bp<0)
				{
					// negative value --> chromosome index
					int chr = -bp - 1;
					if (chr < 0 || chr >= vsCHRdictionary.GetSize())
					{
						ReportError(::Format("Illegal chromosome index (%d)", chr), sFilenameBIN);
						return false;
					}
					CString sCHR = vsCHRdictionary[-bp-1];

					int i = 0;
					for (;i < vCHR_name.GetSize(); i++)
						if (sCHR == vCHR_name[i])
							break;
					if (i<vCHR_name.GetSize())
					{
						nCHR_i0 = vCHR_i0[i];
						nCHR_i1 = vCHR_i1[i];
						last_i0 = -1;
					//	printf("Chr '%s'\n", (const char*)sCHR);
					}
					else
					{
						nCHR_i0 = -1;
						nCHR_i1 = -1;
						last_i0 = -1;
					}
				}
				else
				{
					// non-negative value -> basepair coordinate
					if (nCHR_i0>=0)
					{
						const int bp_start = bp - ww1;
						int i0;

						// int i0_cached = -1;
						if ((last_i0 > nCHR_i0) && (last_i0 < nCHR_i1) && (pSI[last_i0-1].m_Value<bp_start) && (pSI[last_i0].m_Value>=bp_start))
						{
							i0 = last_i0;
							// i0_cached = last_i0;
						}
						else
						{
							i0 = nCHR_i0;
							int i1 = nCHR_i1;
							while (i1 > i0)
							{
								int imid = (i0 + i1) >> 1;
								// printf("imid=%d, d_imid=%d\n", imid, bp_start - pSI[imid].m_Value);
								if (pSI[imid].m_Value < bp_start)
									i0 = imid + 1;
								else
									i1 = imid;
							}
							last_i0 = i0;
						}

						while (i0 < nCHR_i1)
						{
							int d_bp = pSI[i0].m_Value - bp;
							if (d_bp > ww1)
								break;
							if (d_bp >= -ww1)
							{
								const int snp = pSI[i0].m_nIndex;
								pCuts_ww1[snp]++; // region for normalization
								// printf("%d\t%d\n", bp, bp - pSI[i0].m_Value);
								if (__abs(bp - pSI[i0].m_Value) < ww0)
								{
									pCuts_ww0[snp]++; // region for signal
								}
							}
							i0++;
						}
					}
				}
				np--;
			}
			g_FileSystem.FreeBuf(pLB);

			for (int r = 0; r < nSNP; r++)
			{
				// printf("%d\t%d\t%g\t%g\t%g\t%g\n", pCuts_ww0[r], pCuts_ww1[r], double(pCuts_ww0[r])/ww0, double(pCuts_ww1[r])/ww1, (double(pCuts_ww0[r]) / ww0) / (double(pCuts_ww1[r]) / ww1), mG.GetAt(r, s));
				double density= (double(pCuts_ww0[r])/ ww0) / (double(pCuts_ww1[r])/ ww1);
				mData_out.SetAt(r, s, density);
			}
		}
	}
	return true;
}

double GetPi0(const CVector<double> &vP, double pcut)
{
	CVector<double> vPtmp = vP;

	// note: will change order in vP
	if (vPtmp.GetSize()<3)
		return Stat_GetNaN_double();

	double* pP = vPtmp.GetBuffer();
	Stat_Quicksort_double(pP, vP.GetSize(), true);

	int n = 0;
	for (int i = 0; i < vP.GetSize(); i++)
		if (vP[i] > pcut)
			n++;
	double pi0 = (double(n) / vP.GetSize())/(1.0-pcut);
	if (pi0 > 1.0)
		pi0 = 1.0;
	return pi0;

	CVector<double> vX, vY;
	for (double p = 0; p < 0.99; p += 0.01)
	{
		double x = Stat_GetQuantile_double(pP, vP.GetSize(), p);
		if (x > 0.5)
		{
			vX.Add(x);
			vY.Add(p);
			// printf("pi0_data\t%g\t%g\n", x, p);
		}
	}

	double k, m;
	Stat_GetLinreg_double(vX, vY, vX.GetSize(), k, m);
	// printf("k= %g\n", k);
	// printf("m= %g\n", m);

	return k;
}

CString GetSafeName(const CString& s)
{
	CString s_out;
	for (int i = 0; i < s.GetLength(); i++)
		if ((s[i] >= 'a' && s[i] < 'z') || (s[i] >= 'A' && s[i] < 'Z') || (s[i] >= '0' && s[i] < '9'))
			s_out += s[i];
		else
			s_out += '-';
	return s_out;
}

/*
bool CreateApeOutput(CString &sOutputBasename, CTableIndex &idx_ape, CVector<CString> &vRSID, 
	CVector <double> &vmuA1A1,
	CVector <double> &vmuA1A2,
	CVector <double> &vmuA2A2,
	CVector <double> &vK,
	CVector <double> &vR,
	CVector <double> &vP,
	CVector <double> &vP_lower, 
	CVector <double> &vP_upper)
{
	if (idx_ape.IsEmpty())
	{
		ReportError("Perfectos-Ape file is empty", 0);
		return false;
	}
	int cTF = Table_FindCol(idx_ape, "TF_Name");
	if (cTF < 0)
	{
		ReportError("Perfectos-Ape file must contain a column named 'TF_Name'", 0);
		return false;
	}
	int cRSID = Table_FindCol(idx_ape, "RSID");
	if (cRSID < 0)
	{
		ReportError("Perfectos-Ape file must contain a column named 'RSID'", 0);
		return false;
	}
	int cpA1 = Table_FindCol(idx_ape, "ref_P");
	if (cpA1 < 0)
	{
		ReportError("Perfectos - Ape file must contain a column named 'ref_P'", 0);
		return false;
	}
	int cpA2 = Table_FindCol(idx_ape, "alt_P");
	if (cpA2 < 0)
	{
		ReportError("Perfectos - Ape file must contain a column named 'alt_P'", 0);
		return false;
	}

	CVector<Table_SortItem_double> vSI;
	for (int i = 1; i < idx_ape.GetRowCount(); i++)
	{
		Table_SortItem_double si;
		si.m_nIndex = i;
		si.m_sKey = Table_GetAt(idx_ape, i, cTF);
		si.m_Value = 0;
		vSI.Add(si);
	}
	Table_SortItems_double(vSI, true, false);

	Table_SortItem_double* pSI = vSI.GetBuffer();
	for (int i0 = 0; i0 < vSI.GetSize();)
	{
		int i1 = i0 + 1;
		while (i1 < vSI.GetSize() && pSI[i0].m_sKey == pSI[i1].m_sKey)
			i1++;

		if (pSI[i0].m_sKey != "")
		{
			CString sFilename_ape = sOutputBasename + "." + GetSafeName(pSI[i0].m_sKey) + ".ape.txt";
			FILE* pf_ape = OpenOutputFile(sFilename_ape);
			if (!pf_ape)
				return false;

			fprintf(pf_ape, "RSID\tdensity_A1A1\tdensity_A1A2\tdensity_A2A2\tslope\tr\tp1sided\tp2sided\tlog10_fc\n");
			CVector<double> vP_ape;
			for (int i = i0; i < i1; i++)
			{
				CString sRSID = Table_GetAt(idx_ape, pSI[i].m_nIndex, cRSID);
				int j = 0;
				for (; j < vRSID.GetSize(); j++)
					if (vRSID[j] == sRSID)
						break;
				if (j<vRSID.GetSize())
				{
					fprintf(pf_ape, "%s\t%g\t%g\t%g\t%g\t%g", (const char*)vRSID[j], vmuA1A1[j], vmuA1A2[j], vmuA2A2[j], vK[j], vR[j]);

					double pA1 = atof(Table_GetAt(idx_ape, pSI[i].m_nIndex, cpA1));
					double pA2 = atof(Table_GetAt(idx_ape, pSI[i].m_nIndex, cpA2));
					double p1sided;
					if (pA1 < pA2)
					{
						// A1 is gain
						p1sided = vP_lower[j];
					}
					else
					{
						// A2 is gain
						p1sided = vP_upper[j];
					}
					double p2sided = vP[j];
					vP_ape.Add(p2sided);
					fprintf(pf_ape, "\t%g\t%g\t%g", p1sided, p2sided, log10(pA1/pA2));
					fprintf(pf_ape, "\n");
				}
			}
			CloseFile(pf_ape);
			printf("%s\t%d\t%g\n", (const char *)pSI[i0].m_sKey, vP_ape.GetSize(), GetPi0(vP_ape, 0.5));
			i0 = i1;
		}
	}

	return true;
}
*/


int CApp::Main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	bool bOk = false;
	SetDefaultParameters();

	if (ParseCommandLine(argc, argv, envp) && !m_bHelpShown)
	{
		if (m_vNonSwitchArgs.GetSize()<3)
		{
			ReportError("Not enough arguments", 0);
			return 100;
		}

		if (CheckParameters())
		{ 
			if (!m_bMulti)
			{
				CString sErr;
				if (GetVerbosity()>1)
					printf("Loading genotype file...\n");
				CLoadBuf* pLB0 = Table_LoadTabDataset(m_vNonSwitchArgs[0], m_idx_genotypes, sErr);
				if (pLB0)
				{
					CVector<CString> vSampleID;
					CVector<double> vG;
					CString sA1, sA2;
					int nA1A1, nA1A2, nA2A2, nNA;

					bOk = ParseGenotypes(m_vNonSwitchArgs[0], vG, sA1, sA2, vSampleID, nA1A1, nA1A2, nA2A2, nNA);
					if (bOk)
					{
						if (GetVerbosity() > 1)
						{
							printf("rsID= %s\n", (const char*)m_sRSID);
							printf("No. of %s/%s samples= %d\n", (const char*)sA1, (const char*)sA1, nA1A1);
							printf("No. of %s/%s samples= %d\n", (const char*)sA1, (const char*)sA2, nA1A2);
							printf("No. of %s/%s samples= %d\n", (const char*)sA2, (const char*)sA2, nA2A2);
							printf("No. of n/a samples= %d\n", nNA);
						}
						if ((nA1A2 == 0 && nA1A1 == 0) || (nA1A2 == 0 && nA2A2 == 0) || (nA1A1 == 0 && nA2A2 == 0))
						{
							ReportError("Only one genotype represented", 0);
							bOk = false;
						}
					}

					if (bOk)
					{
						int dBP = m_nMaxBP - m_nMinBP;
						int nDPcells = (dBP + m_nBinSpacing) / m_nBinSpacing;
						if (GetVerbosity() > 1)
						{
							printf("-chr: %s\n", (const char*)m_sCHR);
							printf("-minbp: %d\n", m_nMinBP);
							printf("-maxbp: %d\n", m_nMaxBP);
							printf("-binspacing: %d\n", m_nBinSpacing);
							printf("-binsize: %d\n", m_nBinSize);
							printf("-binsize/-binspacing: %d\n", m_nBinSize / m_nBinSpacing);
							printf("No. of bp in region: %d\n", dBP);
							printf("No. of subsampling cells: %d\n", nDPcells);
							printf("-p: %d\n", m_nPermutations);
							printf("-l1: ");
							for (int i = 0; i < m_vLambda1.GetSize(); i++)
							{
								if (i > 0)
									printf(", ");
								printf("%g", m_vLambda1[i]);
							}
							printf("\n");
							printf("-l2: ");
							for (int i = 0; i < m_vLambda2.GetSize(); i++)
							{
								if (i > 0)
									printf(", ");
								printf("%g", m_vLambda2[i]);
							}
							printf("\n");
						}

						CMatrix<int> mCounts;
						mCounts.ReInit(vG.GetSize(), nDPcells, 0);

						bOk = LoadCutSiteFiles(m_vNonSwitchArgs[1], vSampleID, mCounts, m_sCHR, m_nMinBP, m_nMaxBP, m_nBinSpacing);
						if (bOk)
						{
							CString sOutputBasename = g_FileSystem.GetNameWithoutExtension(m_vNonSwitchArgs[2]);

							CMatrix<double> mPI0, mH1proportion;
							mPI0.ReInit(m_vLambda1.GetSize(), m_vLambda2.GetSize(),1.0);
							mH1proportion.ReInit(m_vLambda1.GetSize(), m_vLambda2.GetSize(), 0.0);
							for (int i1 = 0; i1 < m_vLambda1.GetSize(); i1++)
							{
								for (int i2 = 0; i2 < m_vLambda2.GetSize(); i2++)
								{
									const double lambda1 = m_vLambda1[i1];
									const double lambda2 = m_vLambda2[i2];
									if (GetVerbosity() > 1)
										printf("Calculating solution... [%d/%d, l1=%g, l2=%g]\n", i1*m_vLambda2.GetSize()+i2+1, m_vLambda1.GetSize()*m_vLambda2.GetSize(), lambda1, lambda2);

									// calculations for this l1, l2 pair
									CVector<double> vBPstate1_H0;
									if (m_nPermutations > 0)
									{
										CVector<double> vG_perm = vG;
										double* pG_perm = vG_perm.GetBuffer();
										for (int p = 0; p < m_nPermutations; p++)
										{
											const int n = vG_perm.GetSize();
											for (int i = 0; i < n / 2; i++)
											{
												int j0 = rand_int31() % n;
												int j1 = rand_int31() % n;
												double tmp = pG_perm[j0];
												pG_perm[j0] = pG_perm[j1];
												pG_perm[j1] = tmp;
											}
											int nBPstate1_H0;
											bOk = Optimize(lambda1, lambda2, vG_perm, mCounts, m_nMinBP, m_nMaxBP, m_nBinSpacing, CString(""), nBPstate1_H0);
											vBPstate1_H0.Add(nBPstate1_H0);
											// printf("H0\t%d\n", nBPstate1_H0);
										}
									}

									int nBPstate1_H1;
									CString sOutputBasename_l1l2;
									if (m_vLambda1.GetSize() * m_vLambda2.GetSize() == 1)
										sOutputBasename_l1l2 = sOutputBasename;
									else
										sOutputBasename_l1l2 = sOutputBasename + ".l1-" + GetSafeName(m_vsLambda1[i1]) + "-l2-" + GetSafeName(m_vsLambda2[i2]);
									bOk = Optimize(lambda1, lambda2, vG, mCounts, m_nMinBP, m_nMaxBP, m_nBinSpacing, sOutputBasename_l1l2, nBPstate1_H1);

									// calculate pi0
									if (bOk && m_nPermutations > 0)
									{
										double h0_median = Stat_GetMedian_double(vBPstate1_H0);
										double h0_mean = Stat_GetMean_double(vBPstate1_H0);
										double pi0 = Stat_GetNaN_double();

										double H1proportion = double(nBPstate1_H1) / (m_nMaxBP - m_nMinBP);
										if (nBPstate1_H1 > 0)
										{
											pi0 = __min(h0_mean / (nBPstate1_H1), 1.0);
										}
										if (GetVerbosity() > 1)
										{
											printf("pi0= %g\n", pi0);
											// printf("H1%%= %g\n", H1proportion);
										}
										mPI0.SetAt(i1, i2, pi0);
										mH1proportion.SetAt(i1, i2, H1proportion);
									}
								}
							}

							if (m_bExportPI0 && m_nPermutations > 0)
								bOk= ExportPI0(sOutputBasename + ".pi0.txt", m_vLambda1, m_vLambda2, mPI0, mH1proportion);
						}
					}
					g_FileSystem.FreeBuf(pLB0);
				}
			}
			else
			{
				// multi-mode to get basic caQTL stats for multiple variants
				CString sErr;
				printf("Loading genotype file...\n");
				CLoadBuf* pLB0 = Table_LoadTabDataset(m_vNonSwitchArgs[0], m_idx_genotypes, sErr);
				if (pLB0)
				{
					CVector<CString> vSampleID;
					CMatrix<double> mG;
					CVector<CString> vRSID, vCHR, vA1, vA2;
					CVector<int> vBP, vnA1A1, vnA1A2, vnA2A2, vnNA;

					bOk = ParseGenotypes_multi(m_vNonSwitchArgs[0], vSampleID, mG, vRSID, vCHR, vBP, vA1, vA2, vnA1A1, vnA1A2, vnA2A2, vnNA);
					if (bOk && vSampleID.GetSize() > 0 && vCHR.GetSize() > 0)
					{
						CMatrix<double> mCutSiteDensity;
						bOk = LoadCutSiteFiles_multi(m_vNonSwitchArgs[1], vSampleID, vCHR, vBP, m_nMultiWW0, m_nMultiWW1, mCutSiteDensity, mG);
						if (bOk)
						{
							CString sOutputBasename = g_FileSystem.GetNameWithoutExtension(m_vNonSwitchArgs[2]);
							FILE* pf = OpenOutputFile(sOutputBasename + ".txt");
							if (pf)
							{
								const int R = mG.GetHeight();
								const int C = mG.GetWidth();

								CVector<double> vX(C); // valid genotypes
								CVector<double> vY(C); // cut site densities for samples with valid genotypes
								double* pX = vX.GetBuffer();
								double* pY = vY.GetBuffer();

								CVector<double> vP, vR, vP_upper, vP_lower, vP_valid;
								CVector<double> vmuA1A1, vmuA1A2, vmuA2A2, vK;
								fprintf(pf, "RSID\tCHR\tBP\tA1\tA2\tnA1A1\tnA1A2\tnA2A2\tnNA\tmuA1A1\tmuA1A2\tmuA2A2\tslope\tr\tp\n");
								for (int r = 0;r < R; r++)
								{
									fprintf(pf, "%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d", (const char *)vRSID[r], (const char*)vCHR[r], vBP[r], (const char*)vA1[r], (const char*)vA2[r], vnA1A1[r], vnA1A2[r], vnA2A2[r], vnNA[r]);

									double mu_00 = 0;
									double mu_01 = 0;
									double mu_11 = 0;
									int n_00 = 0;
									int n_01 = 0;
									int n_11 = 0;

									int n = 0;
									for (int c = 0; c < C; c++)
									{
										double x = mG.GetAt(r, c);
										if (x >= 0)
										{
											pX[n] = x;
											double d= mCutSiteDensity.GetAt(r, c);
											pY[n] = d;
											if (x < 0.5)
											{
												mu_00 += d;
												n_00++;
											}
											else if (x < 1.5)
											{
												mu_01 += d;
												n_01++;
											}
											else
											{
												mu_11 += d;
												n_11++;
											}
											n++;
										}
									}
									mu_00 /= n_00;
									mu_01 /= n_01;
									mu_11 /= n_11;
									vmuA1A1.Add(mu_00);
									vmuA1A2.Add(mu_01);
									vmuA2A2.Add(mu_11);

									if (n != vnA1A1[r] + vnA1A2[r] + vnA2A2[r])
									{
										ReportError("Internal error: allele counts do not add up correctly", 0);
										return 100;
									}

									double corr, p2sided, p1sided_uppertail, p1sided_lowertail;
									double k;
									if (n > 2)
									{
										if (false)
										{
											// randomly permute genotypes
											for (int i = 0; i < n; i++)
											{
												int j0 = rand_int31() % n;
												int j1 = rand_int31() % n;
												double tmp = pX[j0];
												pX[j0] = pX[j1];
												pX[j1] = tmp;
											}
										}

										double m;
										Stat_GetLinreg_double(pX, pY, n, k, m);

										if (m_bMultiQuantileNormalize)
											Stat_QuantileNormalize_Normal_double(pY, n);
										corr = Stat_GetCorrelation_double(pX, pY, n);
										double t_corr = corr * sqrt(n - 2) / sqrt(1.0 - corr * corr);
										p2sided = 2 * (1.0 - Stat_Tcdf(__abs(t_corr), n - 2));
										p1sided_lowertail = Stat_Tcdf(t_corr, n - 2);
										p1sided_uppertail = 1.0 - Stat_Tcdf(t_corr, n - 2);
										vP_valid.Add(p2sided);
									}
									else
									{
										corr = Stat_GetNaN_double();
										p2sided = Stat_GetNaN_double();
										p1sided_lowertail = Stat_GetNaN_double();
										p1sided_uppertail = Stat_GetNaN_double();
										k = Stat_GetNaN_double();
									}
									vR.Add(corr);
									vP.Add(p2sided); // also adding non-valid P here, to maintain order w.r.t. vRSID
									vP_lower.Add(p1sided_lowertail);
									vP_upper.Add(p1sided_uppertail);
									vK.Add(k);

									fprintf(pf, "\t%g\t%g\t%g\t%g\t%g\t%g", mu_00, mu_01, mu_11, k, corr, p2sided);
									fprintf(pf, "\n");
								}
								CloseFile(pf);

								// export *.param.txt
								double pi0_pcut = 0.5;
								double pi0 = GetPi0(vP_valid, pi0_pcut); // will alter order in vP_valid
								printf("pi0_pcut= %g\n", pi0_pcut);
								printf("pi0= %g\n", pi0);
								FILE* pf_param = OpenOutputFile(sOutputBasename + ".param.txt");
								if (pf_param)
								{
									fprintf(pf_param, "ww0\t%d\n", m_nMultiWW0);
									fprintf(pf_param, "ww1\t%d\n", m_nMultiWW0);
									fprintf(pf_param, "normalization\t%d\n", (int)m_bMultiQuantileNormalize);
									fprintf(pf_param, "n_valid_p\t%d\n", vP.GetSize());
									fprintf(pf_param, "pi0_pcut\t%g\n", pi0_pcut);
									fprintf(pf_param, "pi0\t%g\n", pi0);
									fprintf(pf_param, "pi1\t%g\n", 1.0-pi0);
									CloseFile(pf_param);
								}

								/*
								// export analysis per Perfectos-ape output (if there is such a fourth file specified)
								if (bOk && m_vNonSwitchArgs.GetSize() > 3)
								{
									CTableIndex idx_ape;
									CLoadBuf* pLBape = Table_LoadTabDataset(m_vNonSwitchArgs[3], idx_ape, sErr);
									if (pLBape)
									{
										bOk= CreateApeOutput(sOutputBasename, idx_ape, vRSID, vmuA1A1, vmuA1A2, vmuA2A2, vK, vR, vP, vP_lower, vP_upper);
										g_FileSystem.FreeBuf(pLBape);
									}
									else
									{
										ReportError(sErr, m_vNonSwitchArgs[3]);
										bOk = false;
									}
								}
								*/
							}
						}
					}
					g_FileSystem.FreeBuf(pLB0);
				}
				else
					ReportError("Cannot open genotype file", m_vNonSwitchArgs[0]);
			}
		}
	}

	return bOk ? 0 : 100; 
}
/******************************************************************************/

int main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	CApp app;
	return app.Main(argc, argv, envp);
}
