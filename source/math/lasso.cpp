/*******************************************************************************
 *
 *   lasso.cpp -  Solvers for l1-regularized least squares
 *
 *   Björn Nilsson, 2011
 *
 */

#include "lasso.h"
#include "types/minmax.h"

// accessory functions
void printmatrix_double(const CMatrix<double> &m)
{
	for (int r=0;r<m.GetHeight();r++)
	{
		for (int c=0;c<m.GetWidth();c++)
		{
			if (c>0)
				printf("\t");
			printf("%.2g", m.GetAt(r,c));
		}
		printf("\n");
	}
	printf("\n");
}

void printvector_double(const CVector<double> &v)
{
	for (int i=0;i<v.GetSize();i++)
		printf("%g\t", v[i]);
	printf("\n");
}

// Hager-Zhang conjugate gradient descent
bool LASSO_HZ(const CMatrix<double> &At, const CVector<double> &b, const double lambda, CVector<double> &x, const int nVerbosity)
{
	return false;
}

// Friedman coordinate descent
bool LASSO_CCD_At(const CMatrix<double> &At, const CVector<double> &b, const double lambda, CVector<double> &x, const int nVerbosity)
{
	//
	// Solves ||Ax-b|| + lamdba|x|_1
	//
	// input:	At = data matrix transposed (!)
	//			b  = target vector
	//			x  = initial guess
	//			lambda = l1 penalty
	// output:  returns true if successful, in which case x = solution, returns false if failed
	//

	if (b.IsEmpty())
		return false;
	if (At.IsEmpty())
		return false;
	if (lambda<0.0)
		return false;
	if (b.GetSize()!=At.GetWidth())
		return false;
	if (x.GetSize()!=At.GetHeight())
		return false;

	// compute AtA
	if (nVerbosity>0)
		printf("Computing AtA...");
	const int R= At.GetWidth();
	const int C= At.GetHeight(); // R, C swapped because AtA is A tranposed 
	CMatrix<double> AtA;
	try
	{
		AtA.ReInit(R,R);
	}
	catch(...)
	{
		return false; // out of memory
	}
	for (int i=0;i<R;i++)
		for (int j=0;j<R;j++)
		{
			double *p_i= At.GetPointer(i, 0);
			double *p_j= At.GetPointer(j, 0);
			double s= 0;
			for (int k=0;k<C;k++)
				s += p_i[k]*p_j[k];
			AtA.SetAt(i,j,s);
		}
	if (nVerbosity>0)
		printf("ok\n");

	// compute Atb
	CVector<double> Atb(At.GetHeight());
	const double *p_b= b;
	for (int i=0;i<b.GetSize();i++)
	{
		double *p_i= At.GetPointer(i, 0);
		double s= 0;
		for (int k=b.GetSize();--k>=0;)
			s += p_i[k]*p_b[k];
		Atb.SetAt(i, s);
	}

	return LASSO_CCD_AtA(AtA, Atb, lambda, x, nVerbosity);
}

bool LASSO_CCD_At(const CMatrix<double> &At, const CVector<double> &b, const CVector<double> &vLambda, CVector<double> &x, const int nVerbosity)
{
	//
	// Solves ||Ax-b|| + lamdba|x|_1
	//
	// input:	At = data matrix transposed (!)
	//			b  = target vector
	//			x  = initial guess
	//			lambda = l1 penalty
	// output:  returns true if successful, in which case x = solution, returns false if failed
	//

	if (b.IsEmpty())
		return false;
	if (At.IsEmpty())
		return false;
	if (b.GetSize()!=At.GetWidth())
		return false;
	if (x.GetSize()!=At.GetHeight())
		return false;

	// compute AtA
	if (nVerbosity>0)
		printf("Computing AtA...");
	const int R= At.GetWidth();
	const int C= At.GetHeight(); // R, C swapped because AtA is A tranposed 
	CMatrix<double> AtA;
	try
	{
		AtA.ReInit(R,R);
	}
	catch(...)
	{
		return false; // out of memory
	}
	for (int i=0;i<R;i++)
		for (int j=0;j<R;j++)
		{
			double *p_i= At.GetPointer(i, 0);
			double *p_j= At.GetPointer(j, 0);
			double s= 0;
			for (int k=0;k<C;k++)
				s += p_i[k]*p_j[k];
			AtA.SetAt(i,j,s);
		}
	if (nVerbosity>0)
		printf("ok\n");

	// compute Atb
	CVector<double> Atb(At.GetHeight());
	const double *p_b= b;
	for (int i=0;i<b.GetSize();i++)
	{
		double *p_i= At.GetPointer(i, 0);
		double s= 0;
		for (int k=b.GetSize();--k>=0;)
			s += p_i[k]*p_b[k];
		Atb.SetAt(i, s);
	}

	return LASSO_CCD_AtA(AtA, Atb, vLambda, x, nVerbosity);
}


/*
	from Jerome Friedman's lasso

      subroutine lasso(rho,n,vv,s,thr,x,z,mm)                               
      real rho(n),vv(n,n),s(n),x(n),z(n)                                    
      integer mm(n)                                                         
      call fatmul(2,n,vv,x,s,z,mm)                                          
10320 continue                                                              
10321 continue                                                              
      dlx=0.0                                                               
10330 do 10331 j=1,n                                                        
      xj=x(j)                                                               
      x(j)=0.0                                                              
      t=s(j)+vv(j,j)*xj                                                     
      if (abs(t)-rho(j).gt.0.0) x(j)=sign(abs(t)-rho(j),t)/vv(j,j)          
      if(x(j).eq.xj)goto 10331                                              
      del=x(j)-xj                                                           
      dlx=max(dlx,abs(del))                                                 
      s=s-del*vv(:,j)                                                       
10331 continue                                                              
10332 continue                                                              
      if(dlx.lt.thr)goto 10322                                              
      goto 10321                                                            
10322 continue                                                              
      return                                                                
      end      
*/

bool LASSO_CCD_AtA(const CMatrix<double> &AtA, const CVector<double> &Atb, const double lambda, CVector<double> &x, const int nVerbosity)
{
	//
	// solves the KKT conditions AtA x - Atb + lambda*sign(x) = 0 corresponding to 
	// the minimization problem ||Ax-b||_F^2 + lambda * |x|_1 using coordinate descent
	//
	// input:	AtA= A^t A (symmetric)
	//			Atb= A^t b
	//			x= start vector
	//

	// sanity check
	const int n= AtA.GetWidth();
	if (n!=AtA.GetHeight())
		return false;
	if (n!=Atb.GetSize())
		return false;
	if (n!=x.GetSize())
		return false;
	if (lambda<0.0)
		return false;

	double *p_x= x.GetBuffer(n);

	CVector<double> v_s; 
	double *p_s= v_s.GetBuffer(n);
	for (int i=0;i<n;i++)
	{
		double *p_i= AtA.GetPointer(i,0);
		double s= 0;
		for (int j=n;j-->=0;)
			s += p_i[j]*p_x[j];
		p_s[i]= s - Atb[i];
	}

	const double thr= 1.0e-2; // TODO: argument
	const double eps= 1.0e-8;
	double dlx= thr;
	while (dlx>=thr)
	{
		dlx= 0.0;
		for (int j=n;--j>=0;)
		{
			double xj_old= p_x[j];
			double vv_jj= AtA.GetAt(j,j);

			double xj_new;
			double t= vv_jj*xj_old - p_s[j];
			if (t>lambda)
				xj_new= (t-lambda)/vv_jj;
			else if (t<-lambda)
				xj_new= (t+lambda)/vv_jj;
			else
				xj_new= 0.0;
			
			double dx= xj_new-xj_old;
			if (dx>eps || dx<-eps)
			{
				dlx= __max(dlx, __abs(dx));

				// trick: because AtA is symmetric, we can address rows (stored in linear memory) instead of columns
				double *p_i= AtA.GetPointer(j, 0);
				for (int i=n;--i>=0;)
					p_s[i] += p_i[i]*dx;
				p_x[j]= xj_new;
			}
		}

		if (nVerbosity>2)
			printf("dlx= %g\n", dlx);
	}

	if (nVerbosity>2)
	{
		printf("----- Final x -----\n");
		printvector_double(x);
	}

	return true;
}

bool LASSO_CCD_AtA(const CMatrix<double> &AtA, const CVector<double> &Atb, const CVector<double> &vLambda, CVector<double> &x, const int nVerbosity)
{
	//
	// solves the KKT conditions AtA x - Atb + lambda*sign(x) = 0 corresponding to 
	// the minimization problem ||Ax-b||_F^2 + lambda * |x|_1 using coordinate descent
	//
	// input:	AtA= A^t A (symmetric)
	//			Atb= A^t b
	//			x[]= start vector
	//			vLambda[]= lambda per x[]
	//

	// sanity check
	const int n= AtA.GetWidth();
	if (n!=AtA.GetHeight())
		return false;
	if (n!=Atb.GetSize())
		return false;
	if (n!=x.GetSize())
		return false;
	if (vLambda.GetSize()!=x.GetSize())
		return false;
	for (int i=0;i<vLambda.GetSize();i++)
		if (vLambda[i]<0.0)
			return false;

	double *p_x= x.GetBuffer(n);

	CVector<double> v_s; 
	double *p_s= v_s.GetBuffer(n);
	for (int i=0;i<n;i++)
	{
		double *p_i= AtA.GetPointer(i,0);
		double s= 0;
		for (int j=n;j-->=0;)
			s += p_i[j]*p_x[j];
		p_s[i]= s - Atb[i];
	}

	const double thr= 1.0e-2; // TODO: argument
	const double eps= 1.0e-8;
	double dlx= thr;
	while (dlx>=thr)
	{
		dlx= 0.0;
		for (int j=n;--j>=0;)
		{
			double xj_old= p_x[j];
			double vv_jj= AtA.GetAt(j,j);

			double t= vv_jj*xj_old - p_s[j];
			const double lambda= vLambda[j];
			double xj_new;
			if (t>lambda)
				xj_new= (t-lambda)/vv_jj;
			else if (t<-lambda)
				xj_new= (t+lambda)/vv_jj;
			else
				xj_new= 0.0;
			
			double dx= xj_new-xj_old;
			if (dx>eps || dx<-eps)
			{
				dlx= __max(dlx, __abs(dx));

				// trick: because AtA is symmetric, we can address rows (stored in linear memory) instead of columns
				double *p_i= AtA.GetPointer(j, 0);
				for (int i=n;--i>=0;)
					p_s[i] += p_i[i]*dx;
				p_x[j]= xj_new;
			}
		}

		if (nVerbosity>2)
			printf("dlx= %g\n", dlx);
	}

	if (nVerbosity>2)
	{
		printf("----- Final x -----\n");
		printvector_double(x);
	}

	return true;
}

/*

	// main loop
	const double kkt_tol= 1.0e-4;
	while (true)
	{
		double d_kkt_max= kkt_tol; 
		int i_kkt_max= -1;

		// compute deviations in kkt equations
		for (int i=0;i<n;i++)
		{
			double d_kkt;
			if (x[i]>0)
				d_kkt= p_kkt[i] + lambda;
			else if (x[i]<0)
				d_kkt= p_kkt[i] - lambda;
			else
				d_kkt= __max(0, __abs(p_kkt[i]) - lambda);
		
			if (d_kkt<0)
				d_kkt= -d_kkt;
			if (d_kkt>d_kkt_max)
			{
				d_kkt_max= d_kkt;
				i_kkt_max= i;
			}
		}
		if (i_kkt_max<0)
			return true; // all kkt conditions within tolerance --> done

		// solve most deviating kkt condition
		double x_old= p_x[i_kkt_max];
		double x_new;

		// update kkt conditions
		if (x_old!=x_new)
		{
			double x_delta= x_new-x_old;

			// trick: because AtA is symmetric, we can address rows (stored in linear memory) instead of columns
			double *p_i= AtA.GetPointer(i_kkt_max, 0);
			for (int j=n;--j>=0;)
				v_kkt[j] += x_delta*p_i[j];
		}
	}

	// 
}
*/

/*
bool LASSO_CoordinateDescent2(const CMatrix<double> &AtA, const CVector<double> &Atb, const double lambda, CVector<double> &x)
{
	//
	// solves the KKT conditions AtA x - Atb + lambda*sign(x) = 0 corresponding to 
	// the minimization problem ||Ax-b||_F^2 + lambda * |x|_1 using coordinate descent
	//
	// input:	AtA= A^t A (symmetric)
	//			Atb= A^t b
	//			x= start vector
	//

	printf("Not implemented yet\n");

	return false;

	// sanity check
	const int n= AtA.GetWidth();
	if (n!=AtA.GetHeight())
		return false;
	if (n!=Atb.GetSize())
		return false;
	if (n!=x.GetSize())
		return false;

	double *p_x= x.GetBuffer(n);

	// init KKTs, initiatize active/passive set
	CVector<double> v_kkt; // kkt equation AtAx-Atb without slack
	double *p_kkt= v_kkt.GetBuffer(n);
	for (int i=0;i<n;i++)
	{
		double *p_i= AtA.GetPointer(i,0);
		double s= 0;
		for (int j=n;j-->=0;)
			s += p_i[j]*p_x[j];
		p_kkt[i]= s - Atb[i];
	}

	// main loop
	const double kkt_tol= 1.0e-4;
	while (true)
	{
		double d_kkt_max= kkt_tol; 
		int i_kkt_max= -1;

		// compute deviations in kkt equations
		for (int i=0;i<n;i++)
		{
			double d_kkt;
			if (x[i]>0)
				d_kkt= p_kkt[i] + lambda;
			else if (x[i]<0)
				d_kkt= p_kkt[i] - lambda;
			else
				d_kkt= __max(0, __abs(p_kkt[i]) - lambda);
		
			if (d_kkt<0)
				d_kkt= -d_kkt;
			if (d_kkt>d_kkt_max)
			{
				d_kkt_max= d_kkt;
				i_kkt_max= i;
			}
		}
		if (i_kkt_max<0)
			return true; // all kkt conditions within tolerance --> done

		// solve most deviating kkt condition
		double x_old= p_x[i_kkt_max];
		double x_new;

		// update kkt conditions
		if (x_old!=x_new)
		{
			double x_delta= x_new-x_old;

			// trick: because AtA is symmetric, we can address rows (stored in linear memory) instead of columns
			double *p_i= AtA.GetPointer(i_kkt_max, 0);
			for (int j=n;--j>=0;)
				p_kkt[j] += x_delta*p_i[j];
		}
	}

	// 
}
*/