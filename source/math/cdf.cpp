/*******************************************************************************
 *
 *   cdf.cpp -- 
 *
 *   Björn Nilsson, 2006
 */

// #include <afx.h>
#include "cdf.h"
#include <math.h>
#include "statfun.h"
#include "rand.h"
#include "system/filesystem.h"
#include "spline_cubic.h"
#include <stdio.h>

#define TOL_double (0.0000000000001)
#define EDF_INTERPOLATE

// #define TEST_CDF

void FindIndices(const double *pX, int nX, double x, int &imin, int &imax)
{
	int i0= 0;
	int i1= nX-1;
	while (i1>i0)
	{
		int imid= (i0+i1) >> 1;
		if (pX[imid]<x)
			i0= imid+1; // ==> pX[i0-1]<x
		else
			i1= imid;
	}
	if (i0>0 && pX[i0]>x)
		i0--;
	imin= i0;

	// i0= 0;
	i1= nX-1;
	while (i1>i0)
	{
		int imid= (i0+i1+1) >> 1;
		if (pX[imid]<=x)
			i0= imid; // ==> pX[i0-1]<=x
		else
			i1= imid-1;
	}
	if (i1<nX-1 && pX[i1]<x)
		i1++;
	imax= i1;
}

/*
void TestSearch()
{
	CVector<double> v;
	v.Add(1.0);
	v.Add(2.0);
	v.Add(2.0);
	v.Add(3.0);
	v.Add(4.0);
	v.Add(5.0);
	v.Add(6.0);
	v.Add(7.0);

	int imin, imax;
	FindIndices(v, v.GetSize(), 70.0, imin, imax);
	printf("%d, %d\n", imin, imax);
}
*/

void FindIndices_double(const double *pX, int nX, double x, int &imin, int &imax)
{
	int i0= 0;
	int i1= nX-1;
	while (i1>i0)
	{
		int imid= (i0+i1) >> 1;
		if (pX[imid]<x)
			i0= imid+1; // ==> pX[i0-1]<x
		else
			i1= imid;
	}
	if (i0>0 && pX[i0]>x)
		i0--;
	imin= i0;

	// i0= 0;
	i1= nX-1;
	while (i1>i0)
	{
		int imid= (i0+i1+1) >> 1;
		if (pX[imid]<=x)
			i0= imid; // ==> pX[i0-1]<=x
		else
			i1= imid-1;
	}
	if (i1<nX-1 && pX[i1]<x)
		i1++;
	imax= i1;
}

/******************************************************************************/
// Nilsson algorithm for compact representation
// of empirical distribution functions at controlled
// P-value error rate.

inline int EDF_GetNext(int ki, int N_in, double epsilon)
{
	int k_next= int(N_in - (N_in-ki)/(1.0+epsilon));
	if (k_next>ki)
		return k_next;
	else
		return ki+1;
}

int EDF_GetSize(int N_in, 
				double epsilon)
{
	int i= 0;
	int k= 0;
	while (k<N_in)
	{
		i++;
		k= EDF_GetNext(k, N_in, epsilon);
	}
	return i;
}

void EDF_Reduce(const double *T_in, 
				int N_in, 
				double *T_out, 
				double epsilon)
{
	int i= 0;
	int k= 0;
	while (k<N_in)
	{
		T_out[i++]= T_in[k];
		k= EDF_GetNext(k, N_in, epsilon);
	}
}

double EDF_GetFi(int i, int N_in, double epsilon)
{
	// Assumes zero-based i

	if (i<0)
		return 0.0;

	int j= 0;
	int k= 0;
	while (k<N_in)
	{
		if (j==i)
			return double(k)/(N_in-1);
		k= EDF_GetNext(k, N_in, epsilon);
		j++;
	}

	return 1.0;
}

double EDF_GetFstar(const double T, 
					const double *T_out, 
					int N_out,
					int N_in, 
					double epsilon,
					bool bInterpolate)
{
	int i0, i1; // left and right indices
	FindIndices_double(T_out, N_out, T, i0, i1);

	double Fi0= EDF_GetFi(i0, N_in, epsilon);
	double Fi1= EDF_GetFi(i1, N_in, epsilon);
	if (bInterpolate && Fi1>Fi0+1)
	{
		double dt= (T-T_out[i0])/(T_out[i1]-T_out[i0]);
		return (1.0-dt)*Fi0 + dt*Fi1;
	}
	else
	{
		return Fi0;
	}
}

// For testing
double EDF_GetF(const double *pT, int nT, double T, bool bInterpolate)
{
	if (T<pT[0])
		return 0.0;
	if (T>pT[nT-1])
		return 1.0;
	int i0, i1; // left and right indices
	FindIndices_double(pT, nT, T, i0, i1);

	if (bInterpolate)
	{
		double Fi0= double(i0)/(nT-1);
		double Fi1= double(i1)/(nT-1);
		double dt= (T-pT[i0])/(pT[i1]-pT[i0]);
		return (1.0-dt)*Fi0 + dt*Fi1;
	}
	else
	{
		return double(i0)/(nT-1);
	}
}

// For article
void EDF_GetSequence()
{
	FILE *pf= fopen("C:\\reduced_edf_seq.txt", "wb");
	if (pf)
	{
		int N= 1000000;
		double epsilon= 0.01;
		int k= 0;
		int i= 1;

		while (k<N)
		{
			fprintf(pf, "%d\t%d\t%d\n", i, k, k+1); // 0-based and 1-based indices
			k= EDF_GetNext(k, N, epsilon);
			i++;
		}
		fclose(pf);
	}
}

void EDF_GetTable()
{	
	FILE *pf= fopen("C:\\reduced_edf_sizes.txt", "wb");
	if (pf)
	{
		const double epsilon[]= { 0.01, 0.001, 0.0001 };
		int N_in= 10;
		for (int j=0;j<=8;j++)
		{
			printf("%d\n", N_in);
			int N_exp= int(log(double(N_in))/log(10.0)+0.5);
			if (j>0)
				fprintf(pf, "$10^{%d}$ &", N_exp);
			else
				fprintf(pf, "$N$ & ");
			for (int i=0;i<3;i++)
			{
				if (j==0)
				{
					fprintf(pf, " & $\\epsilon$=%g", epsilon[i]);
				}
				else
				{
					int N_out= EDF_GetSize(N_in, epsilon[i]);
					fprintf(pf, " & %d (%g\\%%)", N_out, 100*(1.0-double(N_out)/N_in));
				}
			}

			if (j>0)
				fprintf(pf, " \\\\\n");
			else
				fprintf(pf, " \\\\\n\\\\\n");
			N_in *= 10;
			N_exp++;
		}
		fclose(pf);
	}
}

/******************************************************************************/

double CReducedCDF::GetFstar(const double T, bool bInterpolate) const
{
	ASSERT(IsValid());
	if (!IsValid())
		return Stat_GetNaN_double();
	if (T<=m_vT[0])
		return 0.0; 
	if (T>=m_vT[m_vT.GetSize()-1])
		return 1.0;
	return EDF_GetFstar(T, m_vT, m_vT.GetSize(), m_nN, m_dEpsilon, bInterpolate);
}

double safelog(double p)
{
	const double tol= 0.0000000001;
	if (p<tol)
		p= tol;
	return log(p)/log(10.0);
}

bool CReducedCDF::CreateErrorFile(CVector<double> &vT)
{
	int N= m_vT.GetSize();
	const double tol= 0.000000001;

	// Entire p=[0,1] range
	N= 10000;
	FILE *pf= fopen("C:\\p_error.txt", "wb");
	if (pf)
	{
		for (int i=0;i<N;i++)
		{
			double p= rand_real1();
			if (p>tol && p<(1.0-tol))
			{
				double T= Stat_GetQuantile_double(vT, vT.GetSize(), 1.0 - p);

				// Monte carlo p
				double p0= 1.0 - EDF_GetF(vT, vT.GetSize(), T, false);
				double p1= 1.0 - EDF_GetF(vT, vT.GetSize(), T, true);

				// Reconstructed p with errors
				double p_star0= 1.0 - GetFstar(T, false);
				double p_star1= 1.0 - GetFstar(T, true);
				double e0= p0 > p_star0 ? p0-p_star0 : p_star0-p0;
				double e1= p1 > p_star1 ? p1-p_star1 : p_star1-p1;

				fprintf(pf, "%g\t%g\t%g\t%g\t%g\n", safelog(p1), safelog(e0), safelog(e0/p0), safelog(e1), safelog(e1/p1));
			}
		}
		fclose(pf);
	}

	// Upper tail up close
	pf= fopen("C:\\p_error_tail.txt", "wb");
	if (pf)
	{
		for (int i=0;i<N;i++)
		{
			double p= rand_real1()/1000.0; // Yields p inside [10E-4, tol].
			if (p>tol && p<(1.0-tol))
			{
				double T= Stat_GetQuantile_double(vT, vT.GetSize(), 1.0 - p);

				// Monte carlo p
				double p0= 1.0 - EDF_GetF(vT, vT.GetSize(), T, false);
				double p1= 1.0 - EDF_GetF(vT, vT.GetSize(), T, true);

				// Reconstructed p with errors
				double p_star0= 1.0 - GetFstar(T, false);
				double p_star1= 1.0 - GetFstar(T, true);
				double e0= p0 > p_star0 ? p0-p_star0 : p_star0-p0;
				double e1= p1 > p_star1 ? p1-p_star1 : p_star1-p1;

				fprintf(pf, "%g\t%g\t%g\t%g\t%g\n", safelog(p1), safelog(e0), safelog(e0/p0), safelog(e1), safelog(e1/p1));
			}
		}
		fclose(pf);
	}
	return true;
}

void CReducedCDF::Uniquify(CVector<double> &vT)
{
	double min= vT[0] + TOL_double;
	for (int i=1;i<vT.GetSize();i++)
	{
		if (vT[i]<min)
			vT.SetAt(i, min);
		min= vT[i] + TOL_double;
	}
}

bool CReducedCDF::Initialize(CVector<double> &vT, double epsilon, int nDF)
{
	// Implements the Nilsson-Fioretos algorithm, for computing 
	// a reduced representation of the edf, while bounding the 
	// relative P-value reconstruction error by epsilon.

	// vT = Sorted vector of test statistics
	// epsilon= relative error bound
	// nDF= Degrees of freedom used in the simulation

	// STEP 1 : Get size of vector

	int N_in= vT.GetSize();
	int N_out= EDF_GetSize(N_in, epsilon);
	// printf("N_in= %d, epsilon= %g, N_out= %d\n", N_in, epsilon, N_out);

	m_dEpsilon= epsilon;
	m_nN= N_in;
	m_nDF= nDF;
	
	try { m_vT.SetSize(N_out); }
	catch(...) { return false; }

	// STEP 2: Get reduced vector

	double *pOut= m_vT.GetBuffer(m_vT.GetSize());
	EDF_Reduce(vT, N_in, pOut, epsilon);

	// STEP 3: Testing
#ifdef TEST_CDF

	int i0, i1;
	double T= m_vT[N_out-1]*0.999;
	printf("Original vector\n");
	FindIndices_double(vT, N_in, T, i0, i1);
	printf("T= %g, i0= %d, i1= %d, T[i0]=%g, T[i1]=%g\n", T, i0, i1, vT[i0], vT[i1]);
	double F= EDF_GetF(vT, N_in, T, true);
	printf("F= %g, Fi0= %g, Fi1=%g\n", F, double(i0)/(N_in-1), double(i1)/(N_in-1));

	printf("Reduced vector\n");
	FindIndices_double(pOut, N_out, T, i0, i1);
	printf("T= %g, i0= %d, i1= %d, T[i0]=%g, T[i1]=%g\n", T, i0, i1, m_vT[i0], m_vT[i1]);
	double Fstar= EDF_GetFstar(T, m_vT, m_vT.GetSize(), m_nN, m_dEpsilon, true);
	printf("F*= %g, Fi0= %g, Fi1= %g\n", Fstar, EDF_GetFi(i0, m_nN, m_dEpsilon), EDF_GetFi(i1, m_nN, m_dEpsilon));

	CreateErrorFile(vT);
	EDF_GetTable();
	EDF_GetSequence();
#endif
	return true;
}

/******************************************************************************/
// Serialization routines

bool CReducedCDF::SaveToFile(const char *aFilename)
{
	if (!IsValid())
		return false;

	int size= 0;
	CVector<double> v;
	v.SetSize(m_vT.GetSize() + 3);
	double *pv0= v.GetBuffer(v.GetSize());
	double *pv= pv0;

	// Store N (size of original vector)
	int *pi= (int *)pv;
	*pi++= m_nN;
	*pi++= m_nDF;
	pv= (double *)pi;
	size += sizeof(int);

	// Store epsilon
	*pv++ = m_dEpsilon;
	size += sizeof(double);

	// Store sample
	double *p= m_vT.GetBuffer(m_vT.GetSize());
	for (int i=0;i<m_vT.GetSize();i++)
	{
		*pv++ = *p++;
		size += sizeof(double);
	}

	return g_FileSystem.CreateBinary(aFilename, (const char *)pv0, size);
}

bool CReducedCDF::LoadFromFile(const char *aFilename)
{
	CLoadBuf *pLB= g_FileSystem.Load(aFilename);
	if (!pLB)
		return false;
	double *p= (double *)pLB->m_aLoadBuf;

	// Get N (size of original vector)
	int *pi= (int *)p;
	m_nN= *pi++;
	m_nDF= *pi++;
	p= (double *)pi;
	
	// Get epsilon 
	m_dEpsilon= *p++;

	// Compute sample size
	int N= (pLB->m_aLoadBuf + pLB->m_BufSize - (char *)p)/sizeof(double);
	m_vT.SetSize(N);

	// Get sample
	for (int i=0;i<N;i++)
		m_vT.SetAt(i, p[i]);

	g_FileSystem.FreeBuf(pLB);
	return true;
}

bool CReducedCDF::IsEqual(const CReducedCDF &other)
{
	if (!other.IsValid() || !IsValid())
		return false;
	if (other.m_dEpsilon!=m_dEpsilon)
	{
		ASSERT(false);
		return false;
	}
	if (other.m_nN!=m_nN)
	{
		ASSERT(false);
		return false;
	}
	if (other.m_nDF!=m_nDF)
	{
		ASSERT(false);
		return false;
	}
	for (int i=0;i<m_vT.GetSize();i++)
	{
		if (m_vT[i]!=other.m_vT[i])
		{
			ASSERT(false);
			return false;
		}
	}
	return true;
}
