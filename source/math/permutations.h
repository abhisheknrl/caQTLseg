/*******************************************************************************
 *
 *   permutations.h -- 
 *
 *   Björn Nilsson, Veta Mera Data HB, 2004
 */


bool Perm_GetGroups(const CVector<int> &vIndex0,
					const CVector<int> &vIndex1,
					const int nLimit,
					const bool bBalanced,
					CVector< CVector<int> > &vG0,
					CVector< CVector<int> > &vG1);

bool Perm_GetTranspositions(const CVector<int> &vIndex0,
							const CVector<int> &vIndex1,
							const int nLimit,
							CVector<int> &vTrans);