#include "cholesky.h"

#include <float.h>
#include <math.h>

bool GetCholesky(CMatrix<double> &A, bool lower)
{
	//
	// Performs in-place Cholesky decomposition
	//

	if (A.GetWidth()!=A.GetHeight())
		return false;
	if (A.IsEmpty())
		return false;
	int n= A.GetHeight();

	/* pseudocode
	for k = 1 to n
		akk = sqrt(akk)
		for i = k + 1 to n
			aik = aik/akk
		end
		for j = k + 1 to n
			for i = j to n
				aij = aij - aik * ajk
			end
		end
	end
	*/

	if (lower)
	{
		// compute lower half of L, leave upper half unchanged
		for (int k=0;k<n;k++)
		{
			double akk= sqrt(A.GetAt(k,k));
			A.SetAt(k,k, akk);
			for (int i=k+1;i<n;i++)
			{
				A.SetAt(i,k,A.GetAt(i,k)/akk);
			}
			for (int j=k+1;j<n;j++)
			{
				double A_jk= A.GetAt(j,k);
				if (A_jk!=0.0)
				{
					for (int i=j;i<n;i++)
					{
						double A_ik= A.GetAt(i,k);
						if (A_ik!=0.0)
						{
							A.SetAt(i,j,A.GetAt(i,j)-A_ik*A_jk);
						}
					}
				}
			}
		}
	}
	else
	{
		// compute upper half of L, leave lower half unchanged
		for (int k=0;k<n;k++)
		{
			double *p_k= A.GetPointer(k, 0);
			p_k[k]= sqrt(p_k[k]);

			double Akk_inv= 1.0/p_k[k];
			for (int i=k+1;i<n;i++)
				p_k[i] *= Akk_inv;

			for (int j=k+1;j<n;j++)
			{
				double *p_j= A.GetPointer(j, 0);
				for (int i=j;i<n;i++)
					p_j[i] -= p_k[i]*p_k[j];
			}
		}
	}

	return true;
}

