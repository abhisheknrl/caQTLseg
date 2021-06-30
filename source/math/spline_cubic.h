/*******************************************************************************
 *
 *   spline_cubic.h -- 
 *
 *   Björn Nilsson, Medical Mathematics, 2006
 *	 Excerpt from the code by John Burkhardt, 2004.
 *
 */

double *spline_cubic_set ( int n, const double *t, const double *y, int ibcbeg, double ybcbeg, 
  int ibcend, double ybcend );

double spline_cubic_val ( int n, const double *t, double tval, const double *y, 
  const double *ypp, double *ypval, double *yppval );

void spline_cubic_val2 ( int n, const double *t, double tval, int *left, const double *y, 
  const double *ypp, double *yval, double *ypval, double *yppval );

void spline_linear_val ( int ndata, const double *t, const double *ydata, 
  double tval, double *yval, double *ypval );
