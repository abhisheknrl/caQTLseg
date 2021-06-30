/*******************************************************************************
 *
 *   permutations.cpp -- 
 *
 *   Björn Nilsson, Veta Mera Data HB, 2004
 */


/******************************************************************************/

#include "types/vector.h"

/******************************************************************************/
// Internal functions
// Warning: Code is not thread-safe

/*===================================================================*/
/* C program for distribution from the Combinatorial Object Server.  */
/* Generates the k-combinations of [n] by transpositions via a       */
/* direct algorithm (see the book for what "direct" means).          */
/* No input error checking.  Assumes 0 <= k <= n <= MAX.             */
/* Outputs both the bitstring and the transposition (x,y) (meaning   */
/* that x leaves the subset and y enters).                           */
/* Algorithm is CAT (Constant Amortized Time).                       */
/* The program can be modified, translated to other languages, etc., */
/* so long as proper acknowledgement is given (author and source).   */     
/* Programmer: Frank Ruskey, 1995.                                   */
/* The latest version of this program may be found at the site       */
/* http://sue.uvic.ca/~cos/inf/comb/CombinationsInfo.html            */
/*===================================================================*/

CVector<bool> g_vBits; // true => position assigned to class 0 (k)
CVector<int> g_vTrans;
CVector< CVector<bool> > g_vComb;
int g_nCombCount;
int g_nCombLimit;
int g_n, g_k;

void Transpose(int i, int j)
{
	if (!g_nCombCount)
		return;
	g_nCombCount--;

	// Convert i and j from base-1 to base-0 indices.
	i--;j--;
	
	// i:th item is moved from Class 0 to Class 1
	// j:th item is moved from Class 1 to Class 0
	ASSERT(i!=j && !g_vBits[i] && g_vBits[j]);
	g_vBits.SetAt(i, true);
	g_vBits.SetAt(j, false);
	g_vTrans.Add(i); 
	g_vTrans.Add(j);
	g_vComb.Add(g_vBits);
}

void NEG(int n, int k);

void GEN(int n, int k) 
{
	if (k>0 && k<n && g_nCombCount)
	{
		GEN(n-1, k);
		if (k == 1)
			Transpose(n, n-1);
		else
			Transpose(n, k-1);
		NEG(n-1, k-1);
	}
};

void NEG(int n, int k) 
{
	if (k>0 && k<n && g_nCombCount)
	{
		GEN(n-1, k-1);
		if (k==1)
			Transpose(n-1, n);
		else
			Transpose(k-1, n);
		NEG(n-1, k);
	}
};

/******************************************************************************/

bool Perm_GetGroups(const CVector<int> &vIndex0,
					const CVector<int> &vIndex1,
					const int nLimit,
					const bool bBalanced,
					CVector< CVector<int> > &vG0,
					CVector< CVector<int> > &vG1)
{
	const int n0= vIndex0.GetSize();
	const int n1= vIndex1.GetSize();
	if (n0==0 || n1==0)
	{
		ASSERT(false);
		return false;
	}

	g_nCombCount= nLimit;
	g_vTrans.SetSize(0);
	g_vComb.SetSize(0);
	g_vBits.SetSize(0);
	for (int i=0;i<n0;i++)
		g_vBits.Add(true);
	for (int i=0;i<n1;i++)
		g_vBits.Add(false);
	g_n= n0+n1;
	g_k= n0;

	// Generate permutations
	GEN(g_n, g_k);

	CVector<int> vIndex= vIndex0;
	for (int i=0;i<vIndex1.GetSize();i++)
		vIndex.Add(vIndex1[i]);

	CVector<int> vTmp0, vTmp1;
	vTmp0.SetSize(n0);
	vTmp1.SetSize(n1);
	for (int i=0;i<g_vComb.GetSize();i++)
	{
		bool bKeep= true;
		if (bBalanced)
		{
			// Assert balanced permutation (of the smallest group)

			int n00= 0; // no of G0 members in this permutation of G0
			for (int j=0;j<n0;j++)
				if (g_vComb[i][j])
					n00++; 

			if (n0<n1)
				bKeep= (n00==(n0>>1) || n00==((n0+1)>>1) );
			else
			{
				int n01= n0-n00; // no of G0 members in permuted G1
				bKeep= (n01==(n1>>1) || n01==((n1+1)>>1) );				
			}
		}

		if (bKeep)
		{
			int g0= 0;
			int g1= 0;
			for (int j=0;j<g_vComb[i].GetSize();j++)
			{
				if (g_vComb[i][j])
					vTmp0.SetAt(g0++, vIndex[j]);
				else
					vTmp1.SetAt(g1++, vIndex[j]);
			}
			ASSERT(g0==n0 && g1==n1);
			vG0.Add(vTmp0);
			vG1.Add(vTmp1);
		}
	}

	return true;
}

bool Perm_GetTranspositions(const CVector<int> &vIndex0,
							const CVector<int> &vIndex1,
							const int nLimit,
							CVector<int> &vTrans)
{
	// Input:
	// vIndex0 and vIndex1 contains the column indices of 
	// the members of class 0 and class 1 respectively.
	// nLimit specifies the maximum number of tranpositions 
	// that should be generated.
	//
	// Output: 
	// vTrans contains the transposition pairs.

	const int n0= vIndex0.GetSize();
	const int n1= vIndex1.GetSize();
	if (n0==0 || n1==0)
	{
		ASSERT(false);
		return false;
	}
	
	g_nCombCount= nLimit;
	g_vTrans.SetSize(0);
	g_vComb.SetSize(0);
	g_vBits.SetSize(0);
	for (int i=0;i<n0;i++)
		g_vBits.Add(true);
	for (int i=0;i<n1;i++)
		g_vBits.Add(false);
	g_n= n0 + n1;
	g_k= n0;

	// Generate
	GEN(g_n, g_k);

	// Copy and convert transpositions to data offsets
	CVector<int> vIndex= vIndex0;
	for (int i=0;i<vIndex1.GetSize();i++)
		vIndex.Add(vIndex1[i]);

	int n= g_vTrans.GetSize();
	int *p0= g_vTrans.GetBuffer(n);
	int *p1= vTrans.GetBuffer(n);
	for (int i=n;--i>=0;)
		*p1++= vIndex[*p0++];

	return true;
}
