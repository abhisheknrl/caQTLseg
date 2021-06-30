#include "mklhelpers.h"
#include "mkl.h"
#include "types/stringfun.h"
#include "system/filehelpers.h"
#include "math/statfun.h"

void printmatrixd(const CMatrix<double> &m)
{
	for (int r=0;r<m.GetHeight();r++)
	{
		for (int c=0;c<m.GetWidth();c++)
		{
			if (c>0)
				printf("\t");
			printf("%f", m.GetAt(r,c));
		}
		printf("\n");
	}
	printf("\n");
}

void printmatrixs(const CMatrix<float> &m)
{
	for (int r=0;r<m.GetHeight();r++)
	{
		for (int c=0;c<m.GetWidth();c++)
		{
			if (c>0)
				printf("\t");
			printf("%f", m.GetAt(r,c));
		}
		printf("\n");
	}
	printf("\n");
}

bool mkl_GetAB(CMatrix<double> &A, CMatrix<double> &B, CMatrix<double> &m_out)
{
	const int m= A.GetHeight();
	const int k= A.GetWidth();
	const int n= B.GetWidth();
	try
	{
		m_out.ReInit(m, n);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A.GetPointer(0,0), k, B.GetPointer(0,0), n, 0.0, m_out.GetPointer(0,0), n);
	return true;
}

bool mkl_GetAAt(CMatrix<double> &A, CMatrix<double> &m_out)
{
	if (A.IsEmpty())
		return false;

	const int R= A.GetHeight();
	const int C= A.GetWidth();
	try
	{
		m_out.ReInit(R, R);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, R, R, C, 1.0, A.GetPointer(0,0), C, A.GetPointer(0,0), C, 0.0, m_out.GetPointer(0,0), R);
	return true;
}

bool mkl_GetAAt(CMatrix<float> &A, CMatrix<float> &m_out)
{
	if (A.IsEmpty())
		return false;

	const int R= A.GetHeight();
	const int C= A.GetWidth();
	try
	{
		m_out.ReInit(R, R);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, R, R, C, 1.0, A.GetPointer(0,0), C, A.GetPointer(0,0), C, 0.0, m_out.GetPointer(0,0), R);
	return true;
}

bool mkl_GetTriangularAAt(CMatrix<double> &A, const char *uplo, CMatrix<double> &m_out)
{
	// assumes A square, lower-triangular (typically, a Cholesky decomposition)
	const int R= A.GetHeight();
	const int C= A.GetWidth();
	if (R==0)
		return false;
	if (R!=C)
		return false;
	if (m_out.GetHeight()!=R || m_out.GetWidth()!=C)
	{
		try
		{
			m_out.ReInit(R, C);
		}
		catch(...)
		{
			ReportError("Out of memory", 0);
			return false;
		}
	}

	// copy matrix
	double *p0= A.GetPointer(0,0);
	double *p1= m_out.GetPointer(0,0);
	for (int i=R*C;--i>=0;)
		*p1++ = *p0++;

	CBLAS_UPLO u;
	if (uplo[0]=='U' || uplo[0]=='u')
		u= CblasLower;
	else if (uplo[0]=='L' || uplo[0]=='l')
		u= CblasUpper;

	cblas_dtrmm(CblasRowMajor, CblasRight, u, CblasTrans, CblasNonUnit, R, R, 1.0, A.GetPointer(0,0), R, m_out.GetPointer(0,0), R); 
	return true;
}

/*
// SLOW! Use mkl_GetTriangularAtA
bool mkl_GetForwardCholesky(CMatrix<double> &A, CMatrix<double> &m_out)
{
	const int R= A.GetHeight();
	const int C= A.GetWidth();
	if (R==0)
		return false;
	if (R!=C)
		return false;
	try
	{
		m_out.ReInit(R, C);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	const int incx= 1;
	for (int i=0;i<R;i++)
	{
		double *p_i= A.GetPointer(i, 0);
		for (int j=0;j<=i;j++)
		{
			double *p_j= A.GetPointer(j, 0);
			int len= j+1;
			double m_ij= ddot(&len, p_i, &incx, p_j, &incx);
			m_out.SetAt(i,j,m_ij);
		}
	}
}
*/

bool mkl_GetTriangularAtA(CMatrix<double> &A, const char *uplo, CMatrix<double> &m_out)
{
	// assumes A square, lower-triangular (typically, a Cholesky decomposition)
	const int R= A.GetHeight();
	const int C= A.GetWidth();
	if (R==0)
		return false;
	if (R!=C)
		return false;
	try
	{
		m_out.ReInit(R, C);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	// copy matrix
	double *p0= A.GetPointer(0,0);
	double *p1= m_out.GetPointer(0,0);
	for (int i=R*C;--i>=0;)
		*p1++ = *p0++;

	CBLAS_UPLO u;
	if (uplo[0]=='U' || uplo[0]=='u')
		u= CblasLower;
	else if (uplo[0]=='L' || uplo[0]=='l')
		u= CblasUpper;

	cblas_dtrmm(CblasRowMajor, CblasLeft, u, CblasTrans, CblasNonUnit, R, R, 1.0, A.GetPointer(0,0), R, m_out.GetPointer(0,0), R); 
	return true;
}

bool mkl_GetAtA(CMatrix<double> &A, CMatrix<double> &m_out)
{
	if (A.IsEmpty())
		return false;

	const int R= A.GetHeight();
	const int C= A.GetWidth();
	try
	{
		m_out.ReInit(C, C);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, C, R, 1.0, A.GetPointer(0,0), C, A.GetPointer(0,0), C, 0.0, m_out.GetPointer(0,0), C);
	return true;
}

bool mkl_GetZeroFilled(CMatrix<double> &A, const char *uplo)
{
	const int R= A.GetHeight();
	const int C= A.GetWidth();
	if (R!=C)
		return false;

	if (uplo[0]=='U' || uplo[0]=='u')
	{
		for (int r=0;r<R;r++)
		{
			double *p= A.GetPointer(r, 0);
			for (int c=r+1;c<C;c++)
				p[c]= 0.0;
		}
		return true;
	}
	
	if (uplo[0]=='L' || uplo[0]=='l')
	{
		for (int r=0;r<R;r++)
		{
			double *p= A.GetPointer(r, 0);
			for (int c=0;c<r;c++)
				p[c]= 0.0;
		}
		return true;
	}

	return false;
}

bool mkl_GetCholesky(CMatrix<double> &A, const char *uplo, bool bZeroFill)
{
	// symmetric and posdef assumed
	const int R= A.GetHeight();
	const int C= A.GetWidth();
	if (R!=C)
		return false;
	if (R==0)
		return false;

	int info;
	dpotrf(uplo, &R, A.GetPointer(0,0), &R, &info); // 'U' yields lower-triangular matrix in C (row-major; Fortran column-major)
	if (info!=0)
	{
//#ifdef _DEBUG
//		ReportError("Cholesky decomposition failed", ::Format("%d", info));
//#endif _DEBUG
		return false;
	}

	if (bZeroFill)
	{
		if (uplo[0]=='U' || uplo[0]=='u')
		{
			for (int r=0;r<R;r++)
			{
				double *p= A.GetPointer(r, 0);
				for (int c=r+1;c<C;c++)
					p[c]= 0.0;
			}
		}
		else if (uplo[0]=='L' || uplo[0]=='l')
		{
			for (int r=0;r<R;r++)
			{
				double *p= A.GetPointer(r, 0);
				for (int c=0;c<r;c++)
					p[c]= 0.0;
			}
		}
	}
	
	return true;
}

bool mkl_GetTriangularInverse(CMatrix<double> &A, const char *uplo, bool bZeroFill)
{
	// inverse of lower-triangular square matrix
	const int R= A.GetHeight();
	const int C= A.GetWidth();
	if (R!=C)
		return false;
	if (R==0)
		return false;

	int info;
	dtrtri(uplo, "N", &R, A.GetPointer(0,0), &R, &info);
	if (info!=0)
	{
		ReportError("Triangular inverse failed", ::Format("%d", info));
		return false;
	}

	if (bZeroFill)
		mkl_GetZeroFilled(A, uplo);

	return true;
}

bool mkl_GetEigenDecomposition(CMatrix<double> &A, CVector<double> &eA, const char *uplo)
{
    //
    //    compute eigenvalues for positive symmetric matrix
    //
    //    input:    A positive definite, symmetric
    //            with uplo = "U" --> lower triangle (C vs Fortran issue)
    //            with uplo = "L" --> upper triangle (C vs Fortran issue)
    //            non-used triangle need not be zeroed
    //
    //    output: eA contains eigenvalues
    //            A contains eigenvectors, one per row (same indices as eigenvalues)
    //

    const int n= A.GetWidth();
    eA.SetSize(n);
    CVector<double> e(n);
    CVector<double> tau(n);
    CVector<double> work(1);

    int lwork= -1;
    int info;

    //
    // workspace query
    dsytrd(uplo, &n, A.GetPointer(0,0), &n, eA.GetBuffer(), e.GetBuffer(), tau.GetBuffer(), work.GetBuffer(), &lwork, &info);
    lwork= int(work.GetAt(0)+0.5);
    work.SetSize(lwork);

    //
    // [d e] <-- tridiagonal form
    dsytrd(uplo, &n, A.GetPointer(0,0), &n, eA.GetBuffer(), e.GetBuffer(), tau.GetBuffer(), work.GetBuffer(), &lwork, &info);
    if (info!=0)
    {
        ReportError("dsytrd failed", ::Format("%d", info));
        return false;
    }

    //
    // A <-- Q matrix
    lwork= -1;
    dorgtr(uplo, &n, A.GetPointer(0,0), &n, tau.GetBuffer(), work.GetBuffer(), &lwork, &info);
    lwork= int(work.GetAt(0)+0.5);
    work.SetSize(lwork);

    dorgtr(uplo, &n, A.GetPointer(0,0), &n, tau.GetBuffer(), work.GetBuffer(), &lwork, &info);
    if (info!=0)
    {
        ReportError("dorgtr failed", ::Format("%d", info));
        return false;
    }

    //
    // [d A] <-- eigenvalues, eigenvectors; assuming pos def matrix
    if (work.GetSize()<4*n)
        work.SetSize(4*n);
    dsteqr("V", &n, eA.GetBuffer(), e.GetBuffer(), A.GetPointer(0,0), &n, work.GetBuffer(), &info);
    if (info!=0)
    {
        ReportError("dpteqr failed", ::Format("%d", info));
        return false;
    }

    return true;
}


double mkl_GetLogdet(const CMatrix<double> &L)
{
	// input: L = Cholesky decomposition of positive definite matrix
	// output: log(det(L)) = sum of log:ed diagonal elements

	const int n= L.GetWidth();
	if (n!=L.GetHeight())
	{
		ASSERT(false);
		return Stat_GetNaN_double();
	}
	if (n==0)
	{
		ASSERT(false);
		return Stat_GetNaN_double();
	}

	double s= 0.0;
	for (int r=0;r<n;r++)
	{
		//printf("%d\t%g\t%g\n", r, L.GetAt(r,r), log(L.GetAt(r,r)));
		s += log(L.GetAt(r,r));
	}
	//printf("2*sum= %g\n", 2*s);
	return 2*s;
}
