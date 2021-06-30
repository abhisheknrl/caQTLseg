#include "sics.h"
#include "mkl.h"
#include "mkl_vml.h"
#include "mklhelpers.h"
#include "math/rand.h"

#include "types/stringfun.h"
#include "system/filehelpers.h"
#include "types/minmax.h"
#include "math/statfun.h"
#include "system/environment.h"

double GetUpdate(const CMatrix<double> &X, const CMatrix<double> &Xi, CMatrix<double> &Xb, const double step, const double lambda, bool bCompute_dKKT)
{
	// update current solution (lower half of X including diagonal)
	// leave empirical covariance matrix unchanged (upper half of X excluding diagonal) 

	// computes:
	// 1. X(lower half) + b*Xi(lower half), excluding diagonal
	// 2. soft-shrinks towards X(upper half=S)
	// 3. stores in full Xb (lower AND upper half)
	const int R= X.GetHeight();
	const int C= X.GetWidth();

	double dKKT_max= 0;
	if (bCompute_dKKT)
	{
		for (int r=0;r<R;r++)
		{
			double *pX= X.GetPointer(r,0);
			double *pS= X.GetPointer(0,r);
			double *pXi= Xi.GetPointer(r,0);
			double *pXbL= Xb.GetPointer(r,0);
			double *pXbU= Xb.GetPointer(0,r);

			//const double bound= max(0, lambda-1.0e-10);
			for (int c=0;c<r;c++)
			{
				double x= pX[c] + step*pXi[c];
				double s= pS[0];
				if (x>=s+lambda) // non-strict ineq important because of dKKT_max computation
					x= s+lambda;
				else if (x<=s-lambda)  // non-strict ineq important
					x= s-lambda;
				else
					dKKT_max= __max(__abs(pXi[c]), dKKT_max);

				pXbL[c]= x;
				pXbU[0]= x;
			
				pS += R;
				pXbU += R;
			}

			pXbL[r]= pX[r]; // do not update diagonal, should be X+lambda from start
		}
	}
	else
	{
		for (int r=0;r<R;r++)
		{
			double *pX= X.GetPointer(r,0);
			double *pS= X.GetPointer(0,r);
			double *pXi= Xi.GetPointer(r,0);
			double *pXbL= Xb.GetPointer(r,0);
			double *pXbU= Xb.GetPointer(0,r);

			//const double bound= max(0, lambda-1.0e-10);
			for (int c=0;c<r;c++)
			{
				double x= pX[c] + step*pXi[c];
				double s= pS[0];
				if (x>=s+lambda)
					x= s+lambda;
				else if (x<=s-lambda)
					x= s-lambda;

				pXbL[c]= x;
				pXbU[0]= x;
			
				pS += R;
				pXbU += R;
			}

			pXbL[r]= pX[r]; // do not update diagonal, should be X+lambda from start
		}
	}

	return dKKT_max;
}

double GetMagicBeta_HB(const double lambda)
{
	// GSE12417, 1000 genes
	// lambda 0.90 --> 0.06-0.08
	// lambda 0.85 --> 0.05-0.09
	// lambda 0.80 --> 0.08-0.10
	// lambda 0.75 --> 0.31
	// lambda 0.70 --> 0.55
	// lambda 0.65 --> 0.69..0.70
	// lambda 0.60 --> 0.75..0.77
	// lambda 0.50 --> 0.80..0.82
	// lambda 0.40 --> 0.86..0.87
	// lambda 0.20 --> 0.86..0.87
	// lambda 0.10 --> 0.90
	// lambda 0.05 --> 0.85

	double beta;
	if (lambda<0.65)
		beta= (0.65-lambda)*0.35+0.75;
	if (lambda>=0.65 && lambda<=0.80)
		beta= 0.75 - (lambda-0.65)*4;
	if (lambda>0.80)
		beta= (1.0-lambda)/2;

	return beta;
}

double GetUpdate_HB(const CMatrix<double> &X, const CMatrix<double> &Xprev, const CMatrix<double> &Xi, CMatrix<double> &Xb, const double alpha, const double beta, const double lambda, bool bCompute_dKKT, double &angle)
{
	// update current solution (lower half of X including diagonal)
	// leave empirical covariance matrix unchanged (upper half of X excluding diagonal) 

	// computes:
	// 1. X(lower half) + b*Xi(lower half), excluding diagonal
	// 2. soft-shrinks towards X(upper half=S)
	// 3. stores in full Xb (lower AND upper half)
	const int R= X.GetHeight();
	const int C= X.GetWidth();

	angle= 0.0;
	double s00= 0;
	double s01= 0;
	double s11= 0;
	double s01_min= 10000000;

	double dKKT_max= 0;
	if (true)
	{
		for (int r=0;r<R;r++)
		{
			double *pX= X.GetPointer(r,0);
			double *pXp= Xprev.GetPointer(r,0);
			double *pS= X.GetPointer(0,r);
			double *pXi= Xi.GetPointer(r,0);
			double *pXbL= Xb.GetPointer(r,0);
			double *pXbU= Xb.GetPointer(0,r);

			//const double bound= max(0, lambda-1.0e-10);
			for (int c=0;c<r;c++)
			{
				double dx_prev= pX[c]-pXp[c];
				double x= pX[c] + alpha*pXi[c] + beta*dx_prev;
				double s= pS[0];
				if (x>=s+lambda) // non-strict ineq important because of dKKT_max computation
					x= s+lambda;
				else if (x<=s-lambda)  // non-strict ineq important
					x= s-lambda;
				else
					dKKT_max= __max(__abs(pXi[c]), dKKT_max);

				double dx= x-pX[c];
				s00 += dx_prev*dx_prev;
				s01 += dx_prev*dx;
				s11 += dx*dx;

				pXbL[c]= x;
				pXbU[0]= x;
			
				pS += R;
				pXbU += R;
			}

			pXbL[r]= pX[r]; // do not update diagonal, should be X+lambda from start
		}
	}
	else
	{
		for (int r=0;r<R;r++)
		{
			double *pX= X.GetPointer(r,0);
			double *pXp= Xprev.GetPointer(r,0);
			double *pS= X.GetPointer(0,r);
			double *pXi= Xi.GetPointer(r,0);
			double *pXbL= Xb.GetPointer(r,0);
			double *pXbU= Xb.GetPointer(0,r);

			//const double bound= max(0, lambda-1.0e-10);
			for (int c=0;c<r;c++)
			{
				double x= (1+beta)*pX[c] + alpha*pXi[c] - beta*pXp[c];
				double s= pS[0];
				if (x>s+lambda)
					x= s+lambda;
				else if (x<=s-lambda)
					x= s-lambda;

				// double dx= x-pX[c];
				// s00 += dx_prev*dx_prev;
				// s01 += dx_prev*dx;
				// s11 += dx*dx;				

				pXbL[c]= x;
				pXbU[0]= x;
			
				pS += R;
				pXbU += R;
			}

			pXbL[r]= pX[r]; // do not update diagonal, should be X+lambda from start
		}
	}

	if (__min(s00, s11)<1.0e-10)
		angle= 0.0;
	else
		angle= s01/sqrt(s00)/sqrt(s11);

	return dKKT_max;
}

double SICS_getpairedupdate(const CMatrix<double> &X0, const CMatrix<double> &Xi0, 
	const CMatrix<double> &X1, const CMatrix<double> &Xi1, 
	const CMatrix<double> &Xb0, const CMatrix<double> &Xb1, 
	const CMatrix<double> &U, const double step, const double lambda)
{
	const int n= X0.GetWidth();

	// update lower half of U
	double dU_max= 0;
	int n_unconstrained= 0;
	for (int r=0;r<n;r++)
	{
		double *pUL= U.GetPointer(r,0);
		double *pUU= U.GetPointer(0,r);
		double *pXi0= Xi0.GetPointer(r,0);
		double *pXi1= Xi1.GetPointer(r,0);
		for (int c=0;c<r;c++)
		{
			double dU= (pXi0[c] - pXi1[c]);
			double u= pUU[c] + step*dU; 
			if (u<-lambda)
				u= -lambda;
			else if (u>lambda)
				u= lambda;
			else
			{
				dU_max= __max(__abs(dU), dU_max);
				n_unconstrained++;
			}

			pUL[c]= u;
			pUU += n;
		}
	}

	// new solution to Xb0
	for (int r=0;r<n;r++)
	{
		double *pU= U.GetPointer(r,0);
		double *pX0L= X0.GetPointer(r,0);
		double *pXb0L= Xb0.GetPointer(r,0);
		double *pXb0U= Xb0.GetPointer(0,r);
		for (int c=0;c<r;c++)
		{
			double x= pX0L[c] + pU[c];
			pXb0L[c]= x;
			pXb0U[0]= x;
			pXb0U += n;
		}
		pXb0L[r]= pX0L[r];
	}

	// new solution to Xb1
	for (int r=0;r<n;r++)
	{
		double *pU= U.GetPointer(r,0);
		double *pX1L= X1.GetPointer(r,0);
		double *pXb1L= Xb1.GetPointer(r,0);
		double *pXb1U= Xb1.GetPointer(0,r);
		for (int c=0;c<r;c++)
		{
			double x= pX1L[c] - pU[c];
			pXb1L[c]= x;
			pXb1U[0]= x;
			pXb1U += n;
		}
		pXb1L[r]= pX1L[r];
	}

	return dU_max;
}

double GetUpdate_noproj(const CMatrix<double> &X, const CMatrix<double> &Xi, CMatrix<double> &Xb, const double step, const double lambda)
{
	// update current solution (lower half of X including diagonal)
	// leave empirical covariance matrix unchanged (upper half of X excluding diagonal) 

	// computes:
	// 1. X(lower half) + b*Xi(lower half), excluding diagonal
	// 2. soft-shrinks towards X(upper half=S)
	// 3. stores in full Xb (lower AND upper half)
	const int R= X.GetHeight();
	const int C= X.GetWidth();

	for (int r=0;r<R;r++)
	{
		double *pX= X.GetPointer(r,0);
		double *pXi= Xi.GetPointer(r,0);
		double *pXbL= Xb.GetPointer(r,0);
		double *pXbU= Xb.GetPointer(0,r);
		for (int c=0;c<r;c++)
		{
			double x= pX[c] + step*pXi[c];
			pXbL[c]= x;
			pXbU[0]= x;			
			pXbU += R;
		}

		pXbL[r]= pX[r]; // do not update diagonal, should be X+lambda from start
	}

	return 0;
}

double GetDuchiGradient(const CMatrix<double> &X, const CMatrix<double> &Xi, const double lambda)
{
	// update current solution (lower half of X including diagonal)
	// leave empirical covariance matrix unchanged (upper half of X excluding diagonal) 

	const double barrier= 0.1;
	const int n= X.GetWidth();
	for (int r=0;r<n;r++)
	{
		double *pX= X.GetPointer(r,0);
		double *pS= X.GetPointer(0,r);
		double *pXiL= Xi.GetPointer(r,0);
		double *pXiU= Xi.GetPointer(0,r);
		for (int c=0;c<r;c++)
		{
			double x= pX[c];
			double s= pS[0];
			double dx= pXiL[c];
			if (x>=s+lambda && dx>=0)
			{
				pXiL[c]= 0;
				pXiU[0]= 0;
			}
			else if (x<=s-lambda && dx<=0)
			{
				pXiL[c]= 0;
				pXiU[0]= 0;
			}
			pS += n;
			pXiU += n;
		}
		pXiL[r]= 0;
	}

	return 0;
}

double GetDualityGap(const CMatrix<double> &S, const CMatrix<double> &Xi, const double lambda)
{
	//
	// computes duality gap: tr(S*Xi) + lambda*|Xi|_1 - n
	//
	// this code assumes upper half of X, excluding diagonal contains original S, and that 
	// the diagonal of S can be recovered by subtracting lambda from the diagonal of X.
	//
	// uses upper half of Xi.
	//

	const int n= S.GetWidth();
	double s0= 0;
	double s1= 0;

	for (int r=0;r<n;r++)
	{
		double *pS= S.GetPointer(r,0);
		double *pXi= Xi.GetPointer(r,0);

		double s0r= 0;
		double s1r= 0;
		for (int c=r+1;c<n;c++)
		{
			s0r += pS[c]*pXi[c];
			if (pXi[c]>=0)
				s1r += pXi[c];
			else
				s1r -= pXi[c];
		}
		s0r *= 2;
		s1r *= 2;
		s0r += (pS[r]-lambda)*pXi[r];
		s1r += __abs(pXi[r]);
		s0 += s0r;
		s1 += s1r;
	}

	return __max(s0 + lambda*s1 - n, 0);
}

double GetDualityGap2(const CMatrix<double> &S, const CMatrix<double> &Xi, const double lambda)
{
	//
	// computes duality gap: tr(S*Xi) + lambda*|Xi|_1 - n
	//
	// this code assumes upper half of X, excluding diagonal contains original S, and that 
	// the diagonal of S can be recovered by subtracting lambda from the diagonal of X.
	//
	// uses upper half of Xi.
	//

	const int n= S.GetWidth();
	double s0= 0;
	double s1= 0;

	for (int r=0;r<n;r++)
	{
		double *pS= S.GetPointer(0,r);
		double *pXi= Xi.GetPointer(r,0);

		double s0r= 0;
		double s1r= 0;
		for (int c=0;c<r;c++)
		{
			s0r += pS[0]*pXi[c];
			if (pXi[c]>=0)
				s1r += pXi[c];
			else
				s1r -= pXi[c];
			pS += n;
		}
		s0r *= 2;
		s1r *= 2;
		s0r += (pS[0]-lambda)*pXi[r];
		s1r += __abs(pXi[r]);
		s0 += s0r;
		s1 += s1r;
	}

	return s0 + lambda*s1 - n;
}

double SICS_GetDualObj_print(const CMatrix<double> &L)
{
	// L = Cholesky decomposition
	double ld= mkl_GetLogdet(L);
	printf("logdet= %g\n", ld);
	return ld + L.GetWidth();
}

double SICS_GetDualObj(const CMatrix<double> &L)
{
	// L = Cholesky decomposition
	double ld= mkl_GetLogdet(L);
	return ld + L.GetWidth();
}

bool SICS_pgd1(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int it)
{
	//
	// input: 
	//		X empirical covariance matrix (positive semidefinite)
	// output: 
	//		X, lower half and diagonal: lower half of regularized covariance matrix.
	//		X, upper half: off-diagonal elements of original covariance matrix.
	//		Xi, inverse of regularized covariance matrix.
	//

	if (lambda<0)
	{
		ReportError("Regularization parameter must be non-negative", 0);
		return false;
	}

	const int n= X.GetWidth();
	if (n!=X.GetHeight())
	{
		ReportError("Input is non-square", 0);
		return false;
	}
	if (n==0)
	{
		ReportError("Input is empty", 0);
		return false;
	}
	if (n==1)
	{
		// 1 x 1 matrix --> solve directly
		X.SetAt(0,0,X.GetAt(0,0)+lambda);
		Xi.ReInit(1,1);
		Xi.SetAt(0,0,1.0/X.GetAt(0,0));

		obj_this= log(X.GetAt(0,0))+1;
		duality_gap= 0.0;
		dKKT_max= 0;
		return true;
	}

	if (!(tol_gap>=0.0 || tol_kkt>=0.0 || tol_it>0 || tol_obj>0))
	{
		ReportError("Convergence criterion not specified", 0);
		return false;
	}

	// allocate memory for matrices (may be large)
	CMatrix<double> Xb;
	try
	{
		Xi.ReInit(n,n);
		Xb.ReInit(n,n);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}
	
	// get first guess (placed in lower half of X)
	for (int r=0;r<n;r++)
	{
		double *p0= X.GetPointer(0,r);
		double *p1= X.GetPointer(r,0);
		for (int c=0;c<r;c++)
		{
			// shrink off-diagonal by an amount lambda towards zero
			if (*p0>lambda)
				*p1= *p0-lambda;
			else if (*p0<-lambda)
				*p1= *p0+lambda;
			else
				*p1= 0.0;
			p1++;
			p0 += n;
		}
		*p1 += lambda; // increase on-diagonal
	}

	bool bFirstGuessOk= false;
	if (!bSkipShrink)
	{
		// copy first initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// compute Cholesky decomposition
		bFirstGuessOk= mkl_GetCholesky(Xb, "U", false);
	}

	// verify that initial guess is positive definite
	if (bFirstGuessOk)
	{
		int i=0;
		for (;i<n && Xb.GetAt(i,i)>0;i++) {}
		bFirstGuessOk= (i==n);
	}

	// if (bSkipShrink || !bFirstGuessOk)
	if (!bFirstGuessOk)
	{
		// initial guess not positive definite, fall back to standard initiation (on-diagonal increased by lambda, off-diagonal untouched)
		for (int r=0;r<n;r++)
		{
			// off-diagonal
			double *p0= X.GetPointer(0,r);
			double *p1= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				*p1= *p0;
				p1++;
				p0 += n;
			}
		
			// on-diagonal already increased (this was done when the first initial guess was constructed above)
		}

		// copy second initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// Cholesky decomposition of second initial guess
		bool bSecondGuessOk= mkl_GetCholesky(Xb, "U", false);

		// Verify that second initial guess is pos def
		if (bSecondGuessOk)
		{
			int i=0;
			for (;i<n && Xb.GetAt(i,i)>0;i++) {}
			bSecondGuessOk= (i==n);
		}

		if (!bSecondGuessOk)
		{
			ReportError("Cannot find positive definite initial guess, make sure input is a correlation matrix and/or try increasing lambda", 0);
			return false;
		}
	}

	// compute (dual) objective
	obj_this= SICS_GetDualObj(Xb); 
//	double obj_prev= obj_this;

	double b0= 1.0;
	it=0;
	int missteps= 0;
	dKKT_max= -1.0;
	duality_gap= -1.0;

	double o_prev= 1.0e10;
	while (tol_it<=0 || it<tol_it)
	{
		//
		// at this point, the following is assumed:
		// (a) lower half of X (including diagonal) contains current solution.
		// (b) upper half of X (excluding diagonal) contains empirical covariance matrix.
		// (c) lower half of Xb contains Cholesky decomposition of current solution.
		//

		if (GetVerbosity()>1 && tol_gap>0 && (n!=2 || !bFirstGuessOk))
		{
			if (it==0)
				printf("n\tdual_value\tduality_gap\tmissteps\n"); // header
		}

		mkl_GetTriangularInverse(Xb, "U", false);
		mkl_GetTriangularAtA(Xb, "U", Xi);

		if (n==2 && bFirstGuessOk)
		{
			duality_gap= 0;
			dKKT_max= 0;
			return true;
		}

		if (tol_gap>0)
			duality_gap= GetDualityGap(X, Xi, lambda)/n;

		if (GetVerbosity()>1 && tol_gap>0)
		{
			//if (it==0)
			//	printf("n\tdual_value\tduality_gap\tmissteps\n"); // header
			printf("%d\t%e\t%e\t%d\n", it, obj_this, duality_gap, missteps);
		}
		if (tol_gap>0 && duality_gap<tol_gap)
			return true;

		if (GetVerbosity()>1 && tol_kkt>0)
		{
			if (it>0)
			{
				if (it==1)
					printf("n\tdual_value\tdual_maxgrad\tmissteps\n"); // header
				printf("%d\t%e\t%e\t%d\n", it, obj_this, dKKT_max, missteps);
			}
		}

		// find an acceptable update
		double obj_next= Stat_GetNegInfty_double();
		double b= b0;
		
		missteps= 0;
		dKKT_max= 0;
		while (obj_next<=obj_this)
		{
			// extremely small step?
			if (b<tol_steplen)
			{
				ReportError("Unable to improve solution, increase lambda", 0);
				if (tol_gap<=0)
					duality_gap= GetDualityGap(X, Xi, lambda)/n; // compute gap on exit
				return false;
			}
			
			// Xb <-- Update
			dKKT_max= GetUpdate(X, Xi, Xb, b, lambda, (tol_kkt>0));
			if (tol_kkt>0 && dKKT_max<tol_kkt)
			{
				duality_gap= GetDualityGap(X, Xi, lambda)/n;  // compute gap on exit
				return true;
			}

			obj_next= obj_this;
			if (mkl_GetCholesky(Xb, "U", false))
				obj_next= SICS_GetDualObj(Xb);

			if (obj_next<=obj_this)
			{
				b /= sqrt(500.0);
				missteps++;
			}
		}

		// solution accepted, copy upper half of Xb to lower half of X (excluding diagonal, which is optimal from the beginning).
		// set upper half of Xb to zero.
//		obj_prev= obj_this;
		obj_this= obj_next;
		for (int r=0;r<n;r++)
		{
			double *p0= Xb.GetPointer(0,r);
			double *p1= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				p1[c]= p0[0];
				p0[0]= 0;
				p0 += n;
			}
		}

		// adaptively adjust initial b (this increases robustness to bad choices for beta and seems to decrease the number of iterations)
		if (missteps==0)
			b0 *= sqrt(sqrt(10.0)); // try a larger step next time
		else
			b0= b; // try last step size next time

		it++;
	}

	if (tol_gap<=0)
		duality_gap= GetDualityGap(X, Xi, lambda)/n; // compute gap on exit
	return true;
}

bool SICS_adm(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int it)
{
	//
	// input: 
	//		X empirical covariance matrix (positive semidefinite)
	// output: 
	//		X, lower half and diagonal: lower half of regularized covariance matrix.
	//		X, upper half: off-diagonal elements of original covariance matrix.
	//		Xi, inverse of regularized covariance matrix.
	//

	if (lambda<0)
	{
		ReportError("Regularization parameter must be non-negative", 0);
		return false;
	}

	const int n= X.GetWidth();
	if (n!=X.GetHeight())
	{
		ReportError("Input is non-square", 0);
		return false;
	}
	if (n==0)
	{
		ReportError("Input is empty", 0);
		return false;
	}
	if (n==1)
	{
		// 1 x 1 matrix --> solve directly
		X.SetAt(0,0,X.GetAt(0,0)+lambda);
		Xi.ReInit(1,1);
		Xi.SetAt(0,0,1.0/X.GetAt(0,0));

		obj_this= log(X.GetAt(0,0))+1;
		duality_gap= 0.0;
		dKKT_max= 0;
		return true;
	}

	if (!(tol_gap>=0.0 || tol_it>0))
	{
		ReportError("-tolgap or -toliter not specified", 0);
		return false;
	}

	// allocate memory for matrices (may be large)
	CMatrix<double> Xb;
	try
	{
		Xi.ReInit(n,n);
		Xb.ReInit(n,n);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}
	
	// get first guess (placed in lower half of X)
	for (int r=0;r<n;r++)
	{
		double *p0= X.GetPointer(0,r);
		double *p1= X.GetPointer(r,0);
		for (int c=0;c<r;c++)
		{
			// shrink off-diagonal by an amount lambda towards zero
			if (*p0>lambda)
				*p1= *p0-lambda;
			else if (*p0<-lambda)
				*p1= *p0+lambda;
			else
				*p1= 0.0;
			p1++;
			p0 += n;
		}
		*p1 += lambda; // increase on-diagonal
	}

	bool bFirstGuessOk= false;
	if (!bSkipShrink)
	{
		// copy first initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// compute Cholesky decomposition
		bFirstGuessOk= mkl_GetCholesky(Xb, "U", false);
	}

	// verify that initial guess is positive definite
	if (bFirstGuessOk)
	{
		int i=0;
		for (;i<n && Xb.GetAt(i,i)>0;i++) {}
		bFirstGuessOk= (i==n);
	}

	// if (bSkipShrink || !bFirstGuessOk)
	if (!bFirstGuessOk)
	{
		// initial guess not positive definite, fall back to standard initiation (on-diagonal increased by lambda, off-diagonal untouched)
		for (int r=0;r<n;r++)
		{
			// off-diagonal
			double *p0= X.GetPointer(0,r);
			double *p1= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				*p1= *p0;
				p1++;
				p0 += n;
			}
		
			// on-diagonal already increased (this was done when the first initial guess was constructed above)
		}

		// copy second initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// Cholesky decomposition of second initial guess
		bool bSecondGuessOk= mkl_GetCholesky(Xb, "U", false);

		// Verify that second initial guess is pos def
		if (bSecondGuessOk)
		{
			int i=0;
			for (;i<n && Xb.GetAt(i,i)>0;i++) {}
			bSecondGuessOk= (i==n);
		}

		if (!bSecondGuessOk)
		{
			ReportError("Cannot find positive definite initial guess, make sure input is a correlation matrix and/or try increasing lambda", 0);
			return false;
		}
	}

	// compute (dual) objective
	obj_this= SICS_GetDualObj(Xb); 
//	double obj_prev= obj_this;

	CVector<double> v_diag(n);
	double *p_diag= v_diag.GetBuffer();

	const double beta= __max(0.15*n, 50); // 10.0; // TODO: empirical rule to select beta

	double b0= beta;
	it=0;
	int missteps= 0;
	dKKT_max= -1.0;
	duality_gap= -1.0;

	// Initialize Z matrix, incl diag
	CVector<double> vZ_diag(n);
	for (int r=0;r<n;r++)
	{
		double *pZ= X.GetPointer(r,0);
		double *pS= X.GetPointer(0,r);
		for (int c=0;c<r;c++)
		{
			pZ[c]= pS[0];
			pS += n;
		}
		vZ_diag.SetAt(r, pS[0]-lambda);
	}

	// initialize X matrix
	mkl_GetTriangularInverse(Xb, "U", false);
	mkl_GetTriangularAtA(Xb, "U", Xi);

	// init Y
	for (int r=0;r<n;r++)
	{
		double *pX= Xi.GetPointer(r,0);
		double *pY= Xb.GetPointer(r,0);
		for (int c=0;c<=r;c++)
			*pY++ = *pX++;
	}

	double o_prev= 1.0e10;
	while (tol_it<=0 || it<tol_it)
	{
		//	Xi contains solution (X matrix), lower X Lagrange multipliers (Z matrix)

		//
		// compute Y
		//		Y0= X - Z/beta;
		//		Y= sign(Y0).*max(abs(Y0) - lambda/beta, 0);
		const double shrink= lambda/beta;
		const double beta_reciprocal= 1.0/beta;
		for (int r=0;r<n;r++)
		{	
			double *pX= Xi.GetPointer(r,0);
			double *pY= Xb.GetPointer(r,0);
			double *pZ= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				double d= pX[c] - pZ[c]*beta_reciprocal;
				if (d>shrink)
					pY[c]= d-shrink;
				else if (d<-shrink)
					pY[c]= d+shrink;
				else
					pY[c]= 0;
			}

			double d= pX[r] - vZ_diag[r]*beta_reciprocal;
			if (d>shrink)
				pY[r]= d-shrink;
			else if (d<-shrink)
				pY[r]= d+shrink;
			else
				pY[r]= 0;
		}

		//printf("Y= \n");
		//printmatrixd(Xb);

		//
		// compute Z
		//	    Z= Z - beta*(X-Y);
		for (int r=0;r<n;r++)
		{	
			double *pX= Xi.GetPointer(r,0);
			double *pY= Xb.GetPointer(r,0);
			double *pZ= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
				pZ[c] -= beta*(pX[c]-pY[c]);
			vZ_diag.SetAt(r, vZ_diag.GetAt(r) - beta*(pX[r]-pY[r]));
		}

		//
		// compute X
		//printf("X= \n");
		//printmatrixd(Xi);

		for (int r=0;r<n;r++)
		{	
			double *pX= Xi.GetPointer(r,0);
			double *pY= Xb.GetPointer(r,0);
			double *pZ= X.GetPointer(r,0);
			double *pS= X.GetPointer(0,r);
			for (int c=0;c<r;c++)
			{
				pX[c]= pY[c] - (pS[0] - pZ[c])/beta;
				pS += n;
			}
			pX[r]= pY[r] - (pS[0] - lambda - vZ_diag.GetAt(r))/beta;
		}

		//printf("A= \n");
		//printmatrixd(Xi);

		CVector<double> eA(n);
		if (!mkl_GetEigenDecomposition(Xi, eA, "U"))
			return false;
		//for (int i=0;i<eA.GetSize();i++)
		//	printf("L_A[%d]= %g\n", i, eA[i]);
		double s_new= 0;
		for (int i=0;i<n;i++)
		{
			double e_new= (eA[i] + sqrt(eA[i]*eA[i] + 4.0/beta))/2.0;
			eA.SetAt(i, e_new);
			s_new += 1.0/e_new;
		}

		// eigenvalues for dual (W)
		printf("e_max= %g, e_min= %g, cond= %g\n", 1.0/Stat_GetMin_double(eA, n), 1.0/Stat_GetMax_double(eA, n), Stat_GetMax_double(eA, n)/Stat_GetMin_double(eA, n));

		// printf("s_new= %g, theoretical= %g\n", s_new, n*(1.0+lambda));

		//for (int i=0;i<eA.GetSize();i++)
		//	printf("L_X[%d]= %g\n", i, eA[i]);

		// compute dual objective function
		obj_this= 0;
		for (int i=0;i<n;i++)
			obj_this += log(eA[i]);
		obj_this= -obj_this + n;
		// printf("o= %g\n", obj_this);
		// TODO: mkl
		for (int r=0;r<n;r++)
		{
			double e= sqrt(eA[r]);
			double *pX= Xi.GetPointer(r, 0);
			for (int i=0;i<n;i++)
				pX[i] *= e;
		}
		mkl_GetAtA(Xi, Xb);
		// TODO: avoid copying by pointer switching
		for (int r=0;r<n;r++)
		{	
			double *p0= Xb.GetPointer(r,0);
			double *p1= Xi.GetPointer(r,0);
			for (int c=0;c<=r;c++)
				*p1++ = *p0++;
		}
		//printmatrixd(Xi);

		if (GetVerbosity()>1 && tol_gap>0)
		{
			duality_gap= GetDualityGap2(X, Xb, lambda)/n;
			if (it==0)
				printf("n\tdual_value\tduality_gap\tmissteps\n"); // header
			printf("%d\t%e\t%e\t%d\n", it, obj_this, duality_gap, missteps);
			if (duality_gap>=0 && duality_gap<tol_gap)
				return true;
		}

		it++;
	}

	return true;
}

bool SICS_pgd2(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int it)
{
	//
	// input: 
	//		X empirical covariance matrix (positive semidefinite)
	// output: 
	//		X, lower half and diagonal: lower half of regularized covariance matrix.
	//		X, upper half: off-diagonal elements of original covariance matrix.
	//		Xi, inverse of regularized covariance matrix.
	//

	if (lambda<0)
	{
		ReportError("Regularization parameter must be non-negative", 0);
		return false;
	}

	const int n= X.GetWidth();
	if (n!=X.GetHeight())
	{
		ReportError("Input is non-square", 0);
		return false;
	}
	if (n==0)
	{
		ReportError("Input is empty", 0);
		return false;
	}
	if (n==1)
	{
		// 1 x 1 matrix --> solve directly
		X.SetAt(0,0,X.GetAt(0,0)+lambda);
		Xi.ReInit(1,1);
		Xi.SetAt(0,0,1.0/X.GetAt(0,0));

		obj_this= log(X.GetAt(0,0))+1;
		duality_gap= 0.0;
		dKKT_max= 0;
		return true;
	}

	if (!(tol_gap>=0.0 || tol_kkt>=0.0 || tol_it>0 || tol_obj>0))
	{
		ReportError("Convergence criterion not specified", 0);
		return false;
	}

	// allocate memory for matrices (may be large)
	CMatrix<double> Xb;
	try
	{
		Xi.ReInit(n,n);
		Xb.ReInit(n,n);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}
	
	// get first guess (placed in lower half of X)
	for (int r=0;r<n;r++)
	{
		double *p0= X.GetPointer(0,r);
		double *p1= X.GetPointer(r,0);
		for (int c=0;c<r;c++)
		{
			// shrink off-diagonal by an amount lambda towards zero
			if (*p0>lambda)
				*p1= *p0-lambda;
			else if (*p0<-lambda)
				*p1= *p0+lambda;
			else
				*p1= 0.0;
			p1++;
			p0 += n;
		}
		*p1 += lambda; // increase on-diagonal
	}

	bool bFirstGuessOk= false;
	if (!bSkipShrink)
	{
		// copy first initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// compute Cholesky decomposition
		bFirstGuessOk= mkl_GetCholesky(Xb, "U", false);
	}

	// verify that initial guess is positive definite
	if (bFirstGuessOk)
	{
		int i=0;
		for (;i<n && Xb.GetAt(i,i)>0;i++) {}
		bFirstGuessOk= (i==n);
	}

	// if (bSkipShrink || !bFirstGuessOk)
	if (!bFirstGuessOk)
	{
		// initial guess not positive definite, fall back to standard initiation (on-diagonal increased by lambda, off-diagonal untouched)
		for (int r=0;r<n;r++)
		{
			// off-diagonal
			double *p0= X.GetPointer(0,r);
			double *p1= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				*p1= *p0;
				p1++;
				p0 += n;
			}
		
			// on-diagonal already increased (this was done when the first initial guess was constructed above)
		}

		// copy second initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// Cholesky decomposition of second initial guess
		bool bSecondGuessOk= mkl_GetCholesky(Xb, "U", false);

		// Verify that second initial guess is pos def
		if (bSecondGuessOk)
		{
			int i=0;
			for (;i<n && Xb.GetAt(i,i)>0;i++) {}
			bSecondGuessOk= (i==n);
		}

		if (!bSecondGuessOk)
		{
			ReportError("Cannot find positive definite initial guess, make sure input is a correlation matrix and/or try increasing lambda", 0);
			return false;
		}
	}

	// compute (dual) objective
	obj_this= SICS_GetDualObj(Xb); 
//	double obj_prev= obj_this;

	const double beta= 1.0; // TODO: empirical rule to select beta

	double b0= beta;
	it=0;
	int missteps= 0;
	dKKT_max= -1.0;
	duality_gap= -1.0;

	double o_prev= 1.0e10;
	while (tol_it<=0 || it<tol_it)
	{
		//
		// at this point, the following is assumed:
		// (a) lower half of X (including diagonal) contains current solution.
		// (b) upper half of X (excluding diagonal) contains empirical covariance matrix.
		// (c) lower half of Xb contains Cholesky decomposition of current solution.
		//

		// inverse of Cholesky decomposition
		mkl_GetTriangularInverse(Xb, "U", false);

		// inverse of current solution
		mkl_GetTriangularAtA(Xb, "U", Xi);

		if (n==2 && bFirstGuessOk)
		{
			duality_gap= 0;
			dKKT_max= 0;
			return true;
		}

		if (tol_gap>0)
			duality_gap= GetDualityGap(X, Xi, lambda)/n;

		if (GetVerbosity()>1 && tol_gap>0)
		{
			if (it==0)
				printf("n\tdual_value\tduality_gap\tmissteps\n"); // header
			printf("%d\t%e\t%e\t%d\n", it, obj_this, duality_gap, missteps);
		}
		if (tol_gap>0 && duality_gap<tol_gap)
			return true;

		if (GetVerbosity()>1 && tol_kkt>0)
		{
			if (it>0)
			{
				if (it==1)
					printf("n\tdual_value\tdual_maxgrad\tmissteps\n"); // header
				printf("%d\t%e\t%e\t%d\n", it, obj_this, dKKT_max, missteps);
			}
		}

		// find an acceptable update
		double obj_next= Stat_GetNegInfty_double();
		double b= b0;
		
		missteps= 0;
		dKKT_max= 0;

		//
		// Duchi line search
		//
		//printmatrixd(X);
		//printmatrixd(Xi);
		GetDuchiGradient(X, Xi, lambda);

		/*
		// form Xi * X
		printmatrixd(Xi);
		cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, n, n, 1.0, X.GetPointer(0,0), n, Xi.GetPointer(0,0), n, 0.0, Xb.GetPointer(0,0), n);
		printmatrixd(Xb);

		double tr_SG= 0;
		for (int i=0;i<n;i++)
			tr_SG += Xb.GetAt(i,i);

		double tr_SGSG= 0;
		for (int r=0;r<n;r++)
		{
			double tr_row= 0;
			double *pXb= Xb.GetPointer(r,0);
			for (int c=0;c<r;c++)
				tr_row += pXb[c]*pXb[c];
			tr_SGSG += tr_row;
		}

		// solve: obj_this + t*tr_SG - 0.5*t*t*tr_SGSG
		double qb= tr_SG;
		double qa= 0.5*tr_SGSG;
		double qc= obj_this;
		//b= (-qb + sqrt(qb*qb-4*qa*qc))/(2*qa);

		printf("qa=%g\n", qa);
		printf("qb=%g\n", qb);
		printf("qc=%g\n", qc);
		printf("b= %g\n", b);
		return false;

		double sX= 0;
		double sXi= 0;
		for (int r=0;r<n;r++)
		{
			double *pX= X.GetPointer(r,0);
			double *pXi= Xi.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				sX += pX[c];
				sXi += pXi[c];
			}
		}
		b= log(sX/sXi)/2;
		printf("b_auto= %g\n", b);
		*/

		/*
		// make b-curve 
		double t= 1.0;
		double t_best;
		double o_best= -1.0;
		for (int i=0;i<20;i++)
		{
			GetUpdate_noproj(X, Xi, Xb, t, lambda);
			double o_i= 0;
			if (mkl_GetCholesky(Xb, "U", false))
				o_i= SICS_GetDualObj(Xb);
			if (o_i>o_best)
			{
				o_best= o_i;
				t_best= t;
			}
			t *= sqrt(sqrt(sqrt(10.0)));
		}
		printf("b_best= %g\n", t_best);

		// find new solution with better objective
		// printf("b_first= %g\n", b);

		// b= t_best;
		*/

		while (obj_next<=obj_this)
		{
			// extremely small step?
			if (b<tol_steplen)
			{
				ReportError("Unable to improve solution, increase lambda", 0);
				if (tol_gap<=0)
					duality_gap= GetDualityGap(X, Xi, lambda)/n; // compute gap on exit
				return false;
			}
			
			GetUpdate_noproj(X, Xi, Xb, b, lambda);
			/*
			if (tol_kkt>0 && dKKT_max<tol_kkt)
			{
				duality_gap= GetDualityGap(X, Xi, lambda)/n;  // compute gap on exit
				return true;
			}
			*/

			obj_next= obj_this;
			if (mkl_GetCholesky(Xb, "U", false))
				obj_next= SICS_GetDualObj(Xb);

			if (obj_next<=obj_this)
			{
				b /= sqrt(500.0);
				missteps++;
			}
		}

		// printf("b_chosen= %g\n", b);

		// solution accepted, copy upper half of Xb to lower half of X (excluding diagonal, which is optimal from the beginning).
		// set upper half of Xb to zero.
//		obj_prev= obj_this;
		obj_this= obj_next;
		for (int r=0;r<n;r++)
		{
			double *p0= Xb.GetPointer(0,r);
			double *p1= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				p1[c]= p0[0];
				p0[0]= 0;
				p0 += n;
			}
		}

		// adaptively adjust initial b (this increases robustness to bad choices for beta and seems to decrease the number of iterations)
		if (missteps==0)
		{
			b0 *= sqrt(sqrt(10.0)); // try a larger step next time
			printf("b0= %g\n", b0);
		}
		else
			b0= b; // try last step size next time

		it++;
	}

	if (tol_gap<=0)
		duality_gap= GetDualityGap(X, Xi, lambda)/n; // compute gap on exit
	return true;
}

double SICS_getpairedgap(	CMatrix<double> &X0, CMatrix<double> &Xi0,
							CMatrix<double> &X1, CMatrix<double> &Xi1, 
							double lambda)
{
	return 1.0; // TODO: compute
}

double SICS_getpairedobjective(CMatrix<double> &Xb0, CMatrix<double> &Xb1)
{
	// assumes X0, X1 contains cholesky decompositions of X0+U and X1-U, respectively
	return mkl_GetLogdet(Xb0) + mkl_GetLogdet(Xb1);
}

bool SICS_copytriangle(CMatrix<double> &X0, CMatrix<double> &X1, bool bLower)
{
	const int n= X0.GetWidth();
	if (X0.GetHeight()!=n || X1.GetWidth()!=n || X1.GetHeight()!=n)
		return false;
	if (bLower)
	{
		for (int r=0;r<n;r++)
		{
			double *pX0= X0.GetPointer(r,0);
			double *pX1= X1.GetPointer(r,0);
			for (int c=0;c<=n;c++)
				*pX1++ = *pX0++;
		}
	}
	else
	{
		for (int r=0;r<n;r++)
		{
			double *pX0= X0.GetPointer(r,r);
			double *pX1= X1.GetPointer(r,r);
			for (int c=r;c<n;c++)
				*pX1++ = *pX0++;
		}
	}
	return true;
}

bool SICS_pgd1_paired(CMatrix<double> &X0, CMatrix<double> &X1, const double lambda, CMatrix<double> &Xi0, CMatrix<double> &Xi1, 
	const double tol_steplen, const double tol_gap, const double tol_kkt, const int tol_it, double &obj_this, double &duality_gap, double &dU_max, int &it)
{
	//
	//	solves
	//	
	//	max(X0, X1) log det X0i + log det X1i - tr S0*X0i - tr S1*X1i - lambda*|| X0i - X1i ||_1  
	//
	//	s.t. X0, X1 positive definite; (X0, X1)_ii = (S0, S1)_ii + lambda
	//
	//	using mkl-accelerated projected gradient-descent
	//

	if (lambda<0)
	{
		ReportError("Regularization parameter must be non-negative", 0);
		return false;
	}

	const int n= X0.GetWidth();
	if (n!=X1.GetWidth() || n!=X0.GetHeight() || n!=X1.GetHeight())
	{
		ReportError("Illegal matrix sizes", 0);
		return false;
	}
	if (n==0)
	{
		ReportError("Input is empty", 0);
		return false;
	}

	// increase diagonal = first guess
	for (int i=0;i<n;i++)
	{
		X0.SetAt(i,i,X0.GetAt(i,i)+lambda);
		X1.SetAt(i,i,X1.GetAt(i,i)+lambda);
	}

	// solve 1x1 matrix directly
	if (n==1)
	{
		X0.SetAt(0,0,X0.GetAt(0,0));
		X1.SetAt(0,0,X1.GetAt(0,0));
		Xi0.ReInit(1,1);
		Xi0.SetAt(0,0,1.0/X0.GetAt(0,0));
		Xi1.ReInit(1,1);
		Xi1.SetAt(0,0,1.0/X1.GetAt(0,0));

		obj_this= log(X0.GetAt(0,0))+log(X1.GetAt(0,0));
		duality_gap= 0.0;
		return true;
	}

	// solve nxn matrix iteratively
	if (!(tol_gap>=0.0 || tol_it>0 || tol_kkt>0))
	{
		ReportError("Convergence criterion not specified", 0);
		return false;
	}

	// allocate memory for matrices (may be large)
	CMatrix<double> Xb0, Xb1, U;
	try
	{
		Xi0.ReInit(n,n);
		Xi1.ReInit(n,n);
		Xb0.ReInit(n,n);
		Xb1.ReInit(n,n);
		U.ReInit(n,n);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	// get initial solution (diagonal already increased)
	SICS_copytriangle(X0, Xb0, true);
	// printf("Xb0=\n");
	// printmatrixd(Xb0);

	bool bFirstGuessOk= mkl_GetCholesky(Xb0, "U", true);
	if (!bFirstGuessOk)
	{
		ReportError("First guess not positive definite (Xb0)", 0);
		return false;
	}
	SICS_copytriangle(X1, Xb1, true);
	// printf("Xb1=\n");
	// printmatrixd(Xb1);

	bFirstGuessOk= mkl_GetCholesky(Xb1, "U", true);
	if (!bFirstGuessOk)
	{
		ReportError("First guess not positive definite (Xb1)", 0);
		return false;
	}
	U.Fill(0);

	// compute (dual) objective
	obj_this= SICS_getpairedobjective(Xb0, Xb1);
//	double obj_prev= obj_this;

	const double beta= 1.0; // TODO: empirical rule to select beta

	double b0= beta;
	it=0;
	int missteps= 0;
	duality_gap= -1.0;
	dU_max= 0.0;
	double o_prev= 1.0e10;
	while (tol_it<=0 || it<tol_it)
	{
		//
		//		at this point, the following is assumed:
		//
		//	(a) lower halves of X0, X1 (including diagonal) contain Cholesky decompositions 
		//		of current solution (X0_upper+U and X1_upper-U, respectively)
		//	(b) upper halves of X0, X1 (including diagonal; because the diagonal
		//		is fixed) contain input correlation matrices.
		//

		mkl_GetTriangularInverse(Xb0, "U", false);
		mkl_GetTriangularAtA(Xb0, "U", Xi0);
		mkl_GetTriangularInverse(Xb1, "U", false);
		mkl_GetTriangularAtA(Xb1, "U", Xi1);

		if (tol_gap>0)
			duality_gap= SICS_getpairedgap(X0, Xi0, X1, Xi1, lambda);

		if (GetVerbosity()>1 && tol_gap>0)
		{
			if (it==0)
				printf("n\tdual_value\tduality_gap\tmissteps\n"); // header
			printf("%d\t%e\t%e\t%d\n", it, obj_this, duality_gap, missteps);
		}
		if (tol_gap>0 && duality_gap<tol_gap)
			return true;

		// find an acceptable update
		double obj_next= Stat_GetNegInfty_double();
		double b= b0;
		
		/*
		printf("X0=\n");
		printmatrixd(X0);
		printf("X1=\n");
		printmatrixd(X1);
		printf("Xi0=\n");
		printmatrixd(Xi0);
		printf("Xi1=\n");
		printmatrixd(Xi1);
		*/

		missteps= 0;
		while (obj_next<=obj_this)
		{
			if (b<tol_steplen)
			{
				ReportError("Unable to improve solution, adjust parameters", 0);
				if (tol_gap<=0)
					duality_gap= SICS_getpairedgap(X0, Xi0, X1, Xi1, lambda);
				return false;
			}
			
			// Xb0, Xb1, U <-- projecteds update
			dU_max= SICS_getpairedupdate(X0, Xi0, X1, Xi1, Xb0, Xb1, U, b, lambda);
			if (tol_kkt>0 && dU_max<tol_kkt)
			{
				if (GetVerbosity()>1)
				{
					if (it==0)
						printf("n\tdual_value\tdKKT_max\tmissteps\n"); // header
					printf("%d\t%e\t%e\t%d\n", it, obj_this, dU_max, missteps);
				}
				duality_gap= SICS_getpairedgap(X0, Xi0, X1, Xi1, lambda); // compute gap on exit
				return true;
			}

			obj_next= obj_this;
			if (mkl_GetCholesky(Xb0, "U", false) && mkl_GetCholesky(Xb1, "U", false))
				obj_next= SICS_getpairedobjective(Xb0, Xb1);
			if (obj_next<=obj_this)
			{
				b /= sqrt(500.0);
				missteps++;
			}
		}

		// accept solution
//		obj_prev= obj_this;
		obj_this= obj_next;

		// copy upper halves of Xb0, Xb1 to lower halves of X0, X1 (excluding diagonals, which are
		// optimal from the beginning), set upper halves of Xb0, Xb1 to zero.
		for (int r=0;r<n;r++)
		{
			double *p0= Xb0.GetPointer(0,r);
			double *p1= X0.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				p1[c]= p0[0];
				p0[0]= 0;
				p0 += n;
			}
		}
		for (int r=0;r<n;r++)
		{
			double *p0= Xb1.GetPointer(0,r);
			double *p1= X1.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				p1[c]= p0[0];
				p0[0]= 0;
				p0 += n;
			}
		}

		if (GetVerbosity()>1 && tol_kkt>0 && missteps==0)
		{
			if (it==0)
				printf("n\tdual_value\tdKKT_max\tmissteps\n"); // header
			printf("%d\t%e\t%e\t%d\n", it, obj_this, dU_max, missteps);
		}

		// adaptively adjust initial b (this increases robustness to bad choices for beta and seems to decrease the number of iterations)
		if (missteps==0)
			b0 *= sqrt(sqrt(10.0)); // try a larger step next time
		else
			b0= b; // try last step size next time

		it++;
	}

	if (tol_gap<=0)
		duality_gap= SICS_getpairedgap(X0, Xi0, X1, Xi1, lambda); // compute gap on exit
	return true;
}

int SICS_GetFirstGuess(CMatrix<double> &X, CMatrix<double> &Xb, const double lambda, bool bSkipShrink)
{
	const int n= X.GetWidth();
	for (int r=0;r<n;r++)
	{
		double *p0= X.GetPointer(0,r);
		double *p1= X.GetPointer(r,0);
		for (int c=0;c<r;c++)
		{
			// shrink off-diagonal by an amount lambda towards zero
			if (*p0>lambda)
				*p1= *p0-lambda;
			else if (*p0<-lambda)
				*p1= *p0+lambda;
			else
				*p1= 0.0;
			p1++;
			p0 += n;
		}
		*p1= 1.0 + lambda; // increase on-diagonal
	}

	bool bFirstGuessOk= false;
	if (!bSkipShrink)
	{
		// copy first initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// compute Cholesky decomposition
		if (mkl_GetCholesky(Xb, "U", false))
			return 1;
	}

	// initial guess not pos def, fall back to standard initiation
	for (int r=0;r<n;r++)
	{
		// off-diagonal
		double *p0= X.GetPointer(0,r);
		double *p1= X.GetPointer(r,0);
		for (int c=0;c<r;c++)
		{
			*p1= *p0;
			p1++;
			p0 += n;
		}
		
		// on-diagonal already increased (this was done when the first initial guess was constructed above)
	}

	// copy second initial guess to Cholesky buffer
	for (int r=0;r<n;r++)
	{
		double *p0= X.GetPointer(r,0);
		double *p1= Xb.GetPointer(r,0);
		int c=0;
		for (;c<=r;c++)
			*p1++= *p0++;
		for (;c<n;c++)
			*p1++= 0;
	}

	// Cholesky decomposition of second initial guess
	if (mkl_GetCholesky(Xb, "U", false))
		return 2;

	ReportError("Cannot find positive definite initial guess, make sure input is a correlation matrix and/or try increasing lambda", 0);
	return -1;
}

bool SICS_OptimizeMagicBeta(const CMatrix<double> &X, double &magicbeta, const double lambda, bool bSkipShrink)
{
	if (magicbeta<=-2)
	{
		if (lambda<=0.60 && X.GetWidth()>=1000) // difficult to estimate beta for high lambda, rely on reference table
		{
			const int n= 10; // subsample size

			printf("Optimizing solver parameters...\n");
			CVector<double> vB;
			double beta_mid= GetMagicBeta_HB(lambda);
			double beta= __max(0, beta_mid-0.15);
			for (;beta<=__min(1.00001, beta_mid+0.05);beta += 0.005)
				vB.Add(beta);
			const int nB= vB.GetSize();

			CVector<int> vB_it(nB);
			vB_it.Fill(0);

			int v= GetVerbosity();
			SetVerbosity(0);
			const int T= 2000; // replicates 
			for (int t=0;t<T;t++)
			{
				// printf("%d/%d\n", t, T);
				CMatrix<double> Xsub(n,n);
				CMatrix<double> Xisub;
				CVector<int> vI;
				for (int i=0;i<n;i++)
					vI.Add(rand_int31() % X.GetWidth());

				for (int r=0;r<n;r++)
					for (int c=0;c<n;c++)
						Xsub.SetAt(r,c,X.GetAt(vI[r], vI[c]));

				for (int r=0;r<n;r++)
					for (int c=0;c<r;c++)
						if (Xsub.GetAt(r,c)!=Xsub.GetAt(c,r))
						{
							SetVerbosity(v);
							ReportError("internal error: subsample not symmetric", 0);
							return false;
						}

				double o, g, k;
				for (int i=0;i<nB;i++)
				{
					// printf("vB[i]= %g\n", vB[i]);
					int it;
					CMatrix<double> Xtmp= Xsub;

					double b= vB[i];
					if (!SICS_pgd3(Xtmp, lambda, Xisub, -1, 1.0e-5, -1, 1.0e-10, -1, o, g, k, bSkipShrink, it, b))
						return false;
					vB_it.SetAt(i, vB_it.GetAt(i) + it);
				}
			}
			SetVerbosity(v);

			int i_min= 0;
			int it_min= vB_it[0];
			int n_min= 1;
			double s_min= vB[0];
			for (int i=1;i<nB;i++)
			{
				if (vB_it[i]<it_min)
				{
					it_min= vB_it[i];
					i_min= i;
					s_min= vB[i];
					n_min= 1;
				}
				else if (vB_it[i]==it_min)
				{
					s_min += vB[i];
					n_min ++;
				}

				double it_avg= double(vB_it[i])/T;
			}

			if (GetVerbosity()>2)
			{
				printf("beta, estimated= %g\n", s_min/n_min);
				printf("beta, reference= %g\n", GetMagicBeta_HB(lambda));
			}

			magicbeta= s_min/n_min;

		}
		else
			magicbeta= -1; // no point estimating beta for small matrices
	}

	return true;
}

bool SICS_pgd3(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int &it, double &magicbeta)
{
	//
	// heavy-ball solver
	//
	// input: 
	//		X empirical covariance matrix (positive semidefinite)
	// output: 
	//		X, lower half and diagonal: lower half of regularized covariance matrix.
	//		X, upper half: off-diagonal elements of original covariance matrix.
	//		Xi, inverse of regularized covariance matrix.
	//

	if (lambda<0)
	{
		ReportError("Regularization parameter must be non-negative", 0);
		return false;
	}

	const int n= X.GetWidth();
	if (n!=X.GetHeight())
	{
		ReportError("Input is non-square", 0);
		return false;
	}
	if (n==0)
	{
		ReportError("Input is empty", 0);
		return false;
	}
	if (n==1)
	{
		// 1 x 1 matrix --> solve directly
		X.SetAt(0,0,X.GetAt(0,0)+lambda);
		Xi.ReInit(1,1);
		Xi.SetAt(0,0,1.0/X.GetAt(0,0));

		obj_this= log(X.GetAt(0,0))+1;
		duality_gap= 0.0;
		dKKT_max= 0;
		return true;
	}

	if (!(tol_gap>=0.0 || tol_kkt>=0.0 || tol_it>0 || tol_obj>0))
	{
		ReportError("Convergence criterion not specified", 0);
		return false;
	}

	// allocate memory for matrices (may be large)
	CMatrix<double> Xb, Xprev;
	try
	{
		Xi.ReInit(n,n);
		Xb.ReInit(n,n);
		Xprev.ReInit(n,n);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	if (!SICS_OptimizeMagicBeta(X, magicbeta, lambda, bSkipShrink))
		return false;

	// get first guess (placed in lower half of X)
	for (int r=0;r<n;r++)
	{
		double *p0= X.GetPointer(0,r);
		double *p1= X.GetPointer(r,0);
		for (int c=0;c<r;c++)
		{
			// shrink off-diagonal by an amount lambda towards zero
			if (*p0>lambda)
				*p1= *p0-lambda;
			else if (*p0<-lambda)
				*p1= *p0+lambda;
			else
				*p1= 0.0;
			p1++;
			p0 += n;
		}
		*p1 += lambda; // increase on-diagonal
	}

	int nFirstGuess= SICS_GetFirstGuess(X, Xb, lambda, bSkipShrink);
	if (nFirstGuess<0)
		return false;

	/*
	bool bFirstGuessOk= false;
	if (!bSkipShrink)
	{
		// copy first initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// compute Cholesky decomposition
		bFirstGuessOk= mkl_GetCholesky(Xb, "U", false);
	}

	// verify that initial guess is positive definite
	if (bFirstGuessOk)
	{
		int i=0;
		for (;i<n && Xb.GetAt(i,i)>0;i++) {}
		bFirstGuessOk= (i==n);
	}

	// if (bSkipShrink || !bFirstGuessOk)
	if (!bFirstGuessOk)
	{
		// initial guess not positive definite, fall back to standard initiation (on-diagonal increased by lambda, off-diagonal untouched)
		for (int r=0;r<n;r++)
		{
			// off-diagonal
			double *p0= X.GetPointer(0,r);
			double *p1= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				*p1= *p0;
				p1++;
				p0 += n;
			}
		
			// on-diagonal already increased (this was done when the first initial guess was constructed above)
		}

		// copy second initial guess to Cholesky buffer
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xb.GetPointer(r,0);
			int c=0;
			for (;c<=r;c++)
				*p1++= *p0++;
			for (;c<n;c++)
				*p1++= 0;
		}

		// Cholesky decomposition of second initial guess
		bool bSecondGuessOk= mkl_GetCholesky(Xb, "U", false);

		// Verify that second initial guess is pos def
		if (bSecondGuessOk)
		{
			int i=0;
			for (;i<n && Xb.GetAt(i,i)>0;i++) {}
			bSecondGuessOk= (i==n);
		}

		if (!bSecondGuessOk)
		{
			ReportError("Cannot find positive definite initial guess, make sure input is a correlation matrix and/or try increasing lambda", 0);
			return false;
		}
	}
	*/

	// compute (dual) objective
	obj_this= SICS_GetDualObj(Xb); 
//	double obj_prev= obj_this;

	//  const double beta= 1.0;

	double hb_alpha0= 2.0;
	it=0;

	int missteps= 0;
//	int missteps_prev= 0;
	dKKT_max= -1.0;
	duality_gap= -1.0;

	Xprev.Fill(0);

	double o_prev= 1.0e10;
	double eta= 1.10; // 3.0;

	double hb_beta0= 0.05;
	double hb_beta= hb_beta0;

	double kappa_smooth= 0;
	double angle;

	while (tol_it<=0 || it<tol_it)
	{
		//
		// at this point, the following is assumed:
		// (a) lower half of X (including diagonal) contains current solution.
		// (b) upper half of X (excluding diagonal) contains empirical covariance matrix.
		// (c) lower half of Xb contains Cholesky decomposition of current solution.
		//

		mkl_GetTriangularInverse(Xb, "U", false);
		mkl_GetTriangularAtA(Xb, "U", Xi);

		if (n==2 && nFirstGuess==1)
		{
			duality_gap= 0;
			dKKT_max= 0;
			return true;
		}

		if (tol_gap>0)
			duality_gap= GetDualityGap(X, Xi, lambda)/n;

		if (GetVerbosity()>1 && tol_gap>0)
		{
			if (it==0)
				printf("n\tdual_value\tduality_gap\tmissteps\tangle\n"); // header
			printf("%d\t%e\t%e\t%d\t%g\n", it, obj_this, duality_gap, missteps, angle);

		//	if (it==0)
		//		printf("n\tdual_value\tduality_gap\tmissteps\n"); // header
		//	printf("%d\t%e\t%e\t%d\n", it, obj_this, duality_gap, missteps);
		}
		if (tol_gap>0 && duality_gap<tol_gap)
			return true;

		if (GetVerbosity()>1 && tol_kkt>0)
		{
			if (it>0)
			{
				if (it==1)
					printf("n\tdual_value\tdual_maxgrad\tmissteps\n"); // header
				printf("%d\t%e\t%e\t%d\n", it, obj_this, dKKT_max, missteps);
			}
		}

		// find an acceptable update
		double obj_next= Stat_GetNegInfty_double();
		double hb_alpha= hb_alpha0;

//		missteps_prev= missteps;
		missteps= 0;
		dKKT_max= 0;

		while (obj_next<=obj_this)
		{
			// extremely small step?
			if (hb_alpha<tol_steplen)
			{
				ReportError("Unable to improve solution, increase lambda", 0);
				if (tol_gap<=0)
					duality_gap= GetDualityGap(X, Xi, lambda)/n; // compute gap on exit
				return false;
			}
			
			// Xb <-- Update
			dKKT_max= GetUpdate_HB(X, Xprev, Xi, Xb, hb_alpha, hb_beta, lambda, (tol_kkt>0), angle);

			if (tol_kkt>0 && dKKT_max<tol_kkt)
			{
				duality_gap= GetDualityGap(X, Xi, lambda)/n;  // compute gap on exit
				return true;
			}

			obj_next= obj_this;
			if (mkl_GetCholesky(Xb, "U", false))
				obj_next= SICS_GetDualObj(Xb);

			if (obj_next<=obj_this)
			{
				hb_alpha /= eta*eta; // 10; // eta*eta*eta;  // 5..10 seems to be a good interval
				if (missteps>1)
					hb_beta= 0; //  /= eta; // = 0; // to guarantee convergence
				missteps++;
			}
		}

		// solution accepted

		// Xprev <-- X
		for (int r=0;r<n;r++)
		{
			double *p0= X.GetPointer(r,0);
			double *p1= Xprev.GetPointer(r,0);
			for (int c=0;c<=r;c++)
				*p1++ = *p0++;
		}
		
		// copy upper half of Xb to lower half of X (excluding diagonal, which is 
		// optimal from the beginning). Set upper half of Xb to zero.
//		obj_prev= obj_this;
		obj_this= obj_next;
		for (int r=0;r<n;r++)
		{
			double *p0= Xb.GetPointer(0,r);
			double *p1= X.GetPointer(r,0);
			for (int c=0;c<r;c++)
			{
				p1[c]= p0[0];
				p0[0]= 0;
				p0 += n;
			}
		}

		// adaptively adjust initial b (this increases robustness to bad choices for beta and seems to decrease the number of iterations)
		if (missteps==0)
		{
			hb_alpha0 *= eta; 
			if (magicbeta<0)
				hb_beta0= GetMagicBeta_HB(lambda);
			else
				hb_beta0= magicbeta;

			hb_beta= hb_beta0;
		}
		else
		{
			hb_alpha0= hb_alpha; // try last step size next time
		}

		it++;
	}

	if (tol_gap<=0)
		duality_gap= GetDualityGap(X, Xi, lambda)/n; // compute gap on exit
	return true;
}


bool SICS_pgd4(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int &it, double &magicbeta)
{
	//
	// heavy-ball continuation solver
	//
	// input: 
	//		X empirical covariance matrix (positive semidefinite)
	// output: 
	//		X, lower half and diagonal: lower half of regularized covariance matrix.
	//		X, upper half: off-diagonal elements of original covariance matrix.
	//		Xi, inverse of regularized covariance matrix.
	//

	if (lambda<0)
	{
		ReportError("Regularization parameter must be non-negative", 0);
		return false;
	}

	const int n= X.GetWidth();
	if (n!=X.GetHeight())
	{
		ReportError("Input is non-square", 0);
		return false;
	}
	if (n==0)
	{
		ReportError("Input is empty", 0);
		return false;
	}
	if (n==1)
	{
		// 1 x 1 matrix --> solve directly
		X.SetAt(0,0,X.GetAt(0,0)+lambda);
		Xi.ReInit(1,1);
		Xi.SetAt(0,0,1.0/X.GetAt(0,0));

		obj_this= log(X.GetAt(0,0))+1;
		duality_gap= 0.0;
		dKKT_max= 0;
		return true;
	}

	if (tol_gap<=0)
	{
		ReportError("continuation solver need -tolgap", 0);
		return false;
	}
	if (tol_kkt>=0.0)
	{
		ReportError("-tolkkt not implemented in continuation solver", 0);
		return false;
	}

	// allocate memory for matrices (may be large)
	CMatrix<double> Xb, Xprev;
	try
	{
		Xi.ReInit(n,n);
		Xb.ReInit(n,n);
		Xprev.ReInit(n,n);
	}
	catch(...)
	{
		ReportError("Out of memory", 0);
		return false;
	}

	if (!SICS_OptimizeMagicBeta(X, magicbeta, lambda, bSkipShrink))
		return false;

	// Check if lambda yields solution that satisifes coarse tol_gap (1.0)
	int nFirstGuess= SICS_GetFirstGuess(X, Xb, lambda, bSkipShrink);
	if (nFirstGuess<0)
		return false;

	obj_this= SICS_GetDualObj(Xb); 
	mkl_GetTriangularInverse(Xb, "U", false);
	mkl_GetTriangularAtA(Xb, "U", Xi);
	duality_gap= GetDualityGap(X, Xi, lambda)/n;
	if (tol_gap>0 && duality_gap>=0 && duality_gap<tol_gap)
	{
		it = 0;
		return true;
	}

	const double tol_gap_coarse= 1.0;
	double lambda_current;
	if (bSkipShrink || duality_gap<tol_gap_coarse || lambda<0.40) // continuation has an adverse effect for low lambda (and is pointless for high lambda as the initial shrink step serves the same purpose)
		lambda_current= lambda;
	else
	{
		lambda_current= __max(lambda, 0.7);
		if (!SICS_GetFirstGuess(X, Xb, lambda_current, bSkipShrink))
			return false;
		obj_this= SICS_GetDualObj(Xb); 
		mkl_GetTriangularInverse(Xb, "U", false);
		mkl_GetTriangularAtA(Xb, "U", Xi);
	}

	it=0;
//	double obj_prev= obj_this;
	double hb_alpha0= 1.0; // 10.0;
	int missteps= 0;
//	int missteps_prev= 0;
	dKKT_max= -1.0;
	duality_gap= -1.0;

	Xprev.Fill(0);

	double o_prev= 1.0e10;
	double eta= 1.10;

	double hb_beta0= 0.05;
	double hb_beta= hb_beta0;

	double angle= 0.0;

	double lambda_delta= 0.1;
	bool bFirstRun= true;
	while (true)
	{
		double tol_gap_current= lambda_current > lambda ? tol_gap_coarse : tol_gap; 

		bool bContinue= false;
		while (!bContinue)
		{
			// at this point, the following is assumed:
			// (a) lower half of X (including diagonal) contains current solution.
			// (b) upper half of X (excluding diagonal) contains empirical covariance matrix.
			// (c) lower half of Xb contains Cholesky decomposition of current solution.

			if (!bFirstRun) // initialized above on first run in pgd4
			{
				mkl_GetTriangularInverse(Xb, "U", false);
				mkl_GetTriangularAtA(Xb, "U", Xi);
			}
			bFirstRun= false;

			if (tol_gap_current>0)
				duality_gap= GetDualityGap(X, Xi, lambda_current)/n;

			if (GetVerbosity()>1 && tol_gap>0)
			{
				if (GetVerbosity()>2)
				{
					if (it==0)
						printf("n\tdual_value\tduality_gap\tn_miss\tlambda\tangle\n"); // header
					printf("%d\t%e\t%e\t%d\t%g\t%g\n", it, obj_this, duality_gap, missteps, lambda_current, angle);
				}
				else
				{
					if (it==0)
						printf("n\tdual_value\tduality_gap\tlambda\n"); // header
					printf("%d\t%e\t%e\t%g\n", it, obj_this, duality_gap, lambda_current);
				}
			}

			if (tol_gap_current>0 && duality_gap<tol_gap_current)
			{
				bContinue= true;
			}

			if (!bContinue)
			{
				if (GetVerbosity()>1 && tol_kkt>0)
				{
					if (it>0)
					{
						if (it==1)
							printf("n\tdual_value\tdual_maxgrad\tmissteps\tlambda\n"); // header
						printf("%d\t%e\t%e\t%d\t%g\n", it, obj_this, dKKT_max, missteps, lambda_current);
					}
				}
			}

			// find an acceptable update
			double obj_next= Stat_GetNegInfty_double();
			double hb_alpha= hb_alpha0;

//			missteps_prev= missteps;
			missteps= 0;
			dKKT_max= 0;

			while (!bContinue && obj_next<=obj_this)
			{
				// extremely small step?
				if (hb_alpha<tol_steplen)
				{
					ReportError("Unable to improve solution, increase lambda", 0);
					if (tol_gap<=0)
						duality_gap= GetDualityGap(X, Xi, lambda_current)/n; // compute gap on exit
					return false;
				}
			
				// Xb <-- Update

				dKKT_max= GetUpdate_HB(X, Xprev, Xi, Xb, hb_alpha, hb_beta, lambda_current, (tol_kkt>0), angle);
				if (tol_kkt>0 && dKKT_max<tol_kkt)
				{
					bContinue= true;
					duality_gap= GetDualityGap(X, Xi, lambda_current)/n;  // compute gap on exit
				}

				if (!bContinue)
				{
					obj_next= obj_this;
					if (mkl_GetCholesky(Xb, "U", false))
						obj_next= SICS_GetDualObj(Xb);

					if (obj_next<=obj_this)
					{
						hb_alpha /= eta*eta*eta; // 10; // eta*eta*eta;  // 5..10 seems to be a good interval
						if (missteps>1)
							hb_beta= 0; //  /= eta; // = 0; // to guarantee convergence
						missteps++;
					}
				}
			}

			if (!bContinue)
			{
				// solution accepted
				// Xprev <-- X
				for (int r=0;r<n;r++)
				{
					double *p0= X.GetPointer(r,0);
					double *p1= Xprev.GetPointer(r,0);
					for (int c=0;c<=r;c++)
						*p1++ = *p0++;
				}
		
				// copy upper half of Xb to lower half of X (excluding diagonal, which is 
				// optimal from the beginning). Set upper half of Xb to zero.
//				obj_prev= obj_this;
				obj_this= obj_next;
				for (int r=0;r<n;r++)
				{
					double *p0= Xb.GetPointer(0,r);
					double *p1= X.GetPointer(r,0);
					for (int c=0;c<r;c++)
					{
						p1[c]= p0[0];
						p0[0]= 0;
						p0 += n;
					}
				}

				// adaptively adjust initial b (this increases robustness to bad choices for beta and seems to decrease the number of iterations)
				if (missteps==0)
				{
					hb_alpha0 *= eta; 
					hb_beta0= GetMagicBeta_HB(lambda);
					hb_beta= hb_beta0;
				}
				else
				{
					hb_alpha0= hb_alpha; // try last step size next time
				}
			} // !bContinue

			it++;
		}

		// proceed to next lambda value
		if (lambda_current>lambda)
		{
			double lambda_prev= lambda_current;

			lambda_delta= __max(lambda_delta, 0.05);
			int i_tries= 10;
			bool bPosDef= false;
			while (!bPosDef && i_tries>0)
			{
				double lambda_next= __max(lambda, lambda_current - lambda_delta);
				if (lambda_delta<0.005)
				{
					lambda_next= lambda;
					i_tries= 1;
				}

				// printf("trying lambda %g\n", lambda_next);

				// first, try to continue by shrinkage (not guaranteed to be pos def)
				for (int r=0;r<n;r++)
				{
					// rescale current solution to make it feasible (is pos def already, will remain pos def), copy to Cholesky buffer
					double *pX= X.GetPointer(r,0);
					double *pS= X.GetPointer(0,r);
					double *pXbL= Xb.GetPointer(r,0);
					double *pXbU= Xb.GetPointer(0,r);

					double th= lambda_next;
					double sh= lambda_current - lambda_next;
					for (int c=0;c<r;c++) // lower triangle of X, Xb
					{
						double dx= pX[c]-pS[0];

						if (dx>th)
							dx= th;
						if (dx<-th)
							dx= -th;

						double x= pS[0] + dx;
						pXbL[c]= x;
						pXbU[0]= x;
						pXbU += n;
						pS += n;
					}
					pXbL[r]= 1.0 + lambda_next; // diagonal of X, Xb
				}

				if (mkl_GetCholesky(Xb, "U", false))
				{
					obj_this= SICS_GetDualObj(Xb);
					lambda_current= lambda_next;
					// if (GetVerbosity()>1)
					//	printf("Reconstraining... (%g)\n", lambda_current);
					bPosDef= true;
				}
				else
				{
					i_tries--;
					lambda_delta /= 2.0;
				}
			}

			if (!bPosDef)
			{
				// continuation failed, basically, jump to target lambda by interpolation
				lambda_current= lambda;

			//	if (GetVerbosity()>1)
			//		printf("Jumping to target constraints... (%g)\n", lambda_current);
				// second, try continuation by interpolation (inefficient)
				for (int r=0;r<n;r++)
				{
					// rescale current solution to make it feasible (is pos def already, will remain pos def), copy to Cholesky buffer
					double scale_factor= lambda_current/lambda_prev;
					double *pX= X.GetPointer(r,0);
					double *pS= X.GetPointer(0,r);
					double *pXbL= Xb.GetPointer(r,0);
					double *pXbU= Xb.GetPointer(0,r);

					int c=0;
					for (;c<r;c++) // lower triangle of X, Xb
					{
						double x= pS[0] + scale_factor*(pX[c]-pS[0]);
						pXbL[c]= x; 
						pXbU[0]= x;
						pS += n;
						pXbU += n;
					}
					pXbL[c]= 1.0 + lambda_current; // diagonal of X, Xb

					c++;
				}

				if (mkl_GetCholesky(Xb, "U", false))
					obj_this= SICS_GetDualObj(Xb);
				else
				{
					ReportError("internal error: rescaled solution not pos def", 0);
					return false;
				}
			}

			// X lower <-- Xb Upper <-- zero
			for (int r=0;r<n;r++)
			{
				double *p0= Xb.GetPointer(0,r);
				double *p1= X.GetPointer(r,0);
				for (int c=0;c<r;c++)
				{
					p1[c]= p0[0];
					p0[0]= 0;
					p0 += n;
				}
				p1[r]= 1.0 + lambda_current; // because choleksy has destroyed diagonal of Xb
			}

			// reset solver states
			hb_beta= 0;
//			missteps_prev= 0;
		}
		else
			break; // done
	}

	if (tol_gap<=0)
		duality_gap= GetDualityGap(X, Xi, lambda_current)/n; // compute gap on exit
	return true;
}

