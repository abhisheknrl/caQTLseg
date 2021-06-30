/*******************************************************************************
 *
 *   spline_cubic.cpp -- 
 *
 *   Björn Nilsson, Medical Mathematics, 2006
 *	 Excerpt from the code by John Burkhardt, 2004.
 *
 */

#include <cstdlib>
#include <cmath>


void dvec_bracket ( int n, const double *x, double xval, int *left, 
  int *right )

//********************************************************************
//
//  Purpose:
//
//    DVEC_BRACKET searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    nearest to or containing the given value.
//
//    It is always true that RIGHT = LEFT+1.
//
//    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
//      XVAL   < X[0] < X[1];
//    If X(1) <= XVAL < X[N-1], then
//      X[LEFT-1] <= XVAL < X[RIGHT-1]; 
//    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
//      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
//
//    For consistency, this routine computes indices RIGHT and LEFT
//    that are 1-based, although it would be more natural in C and
//    C++ to use 0-based values.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input, double X[N], an array that has been sorted into ascending order.
//
//    Input, double XVAL, a value to be bracketed.
//
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
  int i;

  for ( i = 2; i <= n - 1; i++ ) 
  {
    if ( xval < x[i-1] ) 
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
}

void spline_linear_val ( int ndata, const double *tdata, const double *ydata, 
  double tval, double *yval, double *ypval )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.
//
//  Discussion:
//
//    Because of the extremely simple form of the linear spline,
//    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
//    evaluate the spline at any point.  No processing of the data
//    is required.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points defining the spline.
//
//    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
//    and dependent variables at the data points.  The values of TDATA should
//    be distinct and increasing.
//
//    Input, double TVAL, the point at which the spline is to be evaluated.
//
//    Output, double *YVAL, *YPVAL, the value of the spline and its first
//    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
//    equal to TDATA(I) for some I.
//
{
  int left;
  int right;
//
//  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
//  nearest to, TVAL.
//
  dvec_bracket ( ndata, tdata, tval, &left, &right );
//
//  Now evaluate the piecewise linear function.
//
  *ypval = ( ydata[right-1] - ydata[left-1] ) 
         / ( tdata[right-1] - tdata[left-1] );

  *yval = ydata[left-1] +  ( tval - tdata[left-1] ) * (*ypval);

  return;
}

void dvec_bracket3 ( int n, const double *t, double tval, int *left )

//********************************************************************
//
//  Purpose:
//
//    DVEC_BRACKET3 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[N-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    For consistency with other versions of this routine, the
//    value of LEFT is assumed to be a 1-based index.  This is
//    contrary to the typical C and C++ convention of 0-based indices.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, double T[N], an array that has been sorted into ascending order.
//
//    Input, double TVAL, a value to be bracketed by entries of T.
//
//    Input/output, int *LEFT.
//
//    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
//    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
//    is searched first, followed by the appropriate interval to the left
//    or right.  After that, a binary search is used.
//
//    On output, LEFT is set so that the interval [ T[LEFT-1], T[LEFT] ]
//    is the closest to TVAL; it either contains TVAL, or else TVAL
//    lies outside the interval [ T[0], T[N-1] ].
//
{
  int high;
  int low;
  int mid;
//  
//  Check the input data.
//
  if ( n < 2 ) 
  {
	  /*
    cout << "\n";
    cout << "DVEC_BRACKET3 - Fatal error!\n";
    cout << "  N must be at least 2.\n";
	*/
    exit ( 1 );
  }
//
//  If *LEFT is not between 1 and N-1, set it to the middle value.
//
  if ( *left < 1 || n - 1 < *left ) 
  {
    *left = ( n + 1 ) / 2;
  }

//
//  CASE 1: TVAL < T[*LEFT]:
//  Search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-1.
//
  if ( tval < t[*left] ) 
  {

    if ( *left == 1 ) 
    {
      return;
    }
    else if ( *left == 2 ) 
    {
      *left = 1;
      return;
    }
    else if ( t[*left-2] <= tval )
    {
      *left = *left - 1;
      return;
    }
    else if ( tval <= t[1] ) 
    {
      *left = 1;
      return;
    }
// 
//  ...Binary search for TVAL in (T[I-1],T[I]), for I = 2 to *LEFT-2.
//
    low = 2;
    high = *left - 2;

    for (;;)
    {

      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid-1] <= tval ) 
      {
        low = mid;
      }
      else 
      {
        high = mid - 1;
      }

    }
  }
// 
//  CASE 2: T[*LEFT] < TVAL:
//  Search for TVAL in (T[I-1],T[I]) for intervals I = *LEFT+1 to N-1.
//
  else if ( t[*left] < tval ) 
  {

    if ( *left == n - 1 ) 
    {
      return;
    }
    else if ( *left == n - 2 ) 
    {
      *left = *left + 1;
      return;
    }
    else if ( tval <= t[*left+1] )
    {
      *left = *left + 1;
      return;
    }
    else if ( t[n-2] <= tval ) 
    {
      *left = n - 1;
      return;
    }
// 
//  ...Binary search for TVAL in (T[I-1],T[I]) for intervals I = *LEFT+2 to N-2.
//
    low = *left + 2;
    high = n - 2;

    for ( ; ; ) 
    {

      if ( low == high ) 
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid-1] <= tval ) 
      {
        low = mid;
      }
      else 
      {
        high = mid - 1;
      }
    }
  }
//
//  CASE 3: T[*LEFT-1] <= TVAL <= T[*LEFT]:
//  T is just where the user said it might be.
//
  else 
  {
  }

  return;
}
//******************************************************************************


double *d3_np_fs ( int n, double a[], double b[] )

//**********************************************************************
//
//  Purpose:
//
//    D3_NP_FS factors and solves a D3 system.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    This algorithm requires that each diagonal entry be nonzero.
//    It does not use pivoting, and so can fail on systems that
//    are actually nonsingular.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[3*N].
//    On input, the nonzero diagonals of the linear system.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input, double B[N], the right hand side.
//
//    Output, double D3_NP_FS[N], the solution of the linear system.
//    This is NULL if there was an error because one of the diagonal
//    entries was zero.
//
{
  int i;
  double *x;
  double xmult;
//
//  Check.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      return NULL;
    }
  }
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( i = 1; i < n; i++ )
  {
    xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
    a[1+i*3] = a[1+i*3] - xmult * a[0+i*3];
    x[i] = x[i] - xmult * x[i-1];
  }

  x[n-1] = x[n-1] / a[1+(n-1)*3];
  for ( i = n-2; 0 <= i; i-- )
  {
    x[i] = ( x[i] - a[0+(i+1)*3] * x[i+1] ) / a[1+i*3];
  }

  return x;
}
//******************************************************************************

double *spline_cubic_set ( int n, const double *t, const double *y, int ibcbeg, double ybcbeg, 
  int ibcend, double ybcend )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
//
//  Discussion:
//
//    For data interpolation, the user must call SPLINE_SET to determine
//    the second derivative data, passing in the data to be interpolated,
//    and the desired boundary conditions.
//
//    The data to be interpolated, plus the SPLINE_SET output, defines
//    the spline.  The user may then call SPLINE_VAL to evaluate the
//    spline at any point.
//
//    The cubic spline is a piecewise cubic polynomial.  The intervals
//    are determined by the "knots" or abscissas of the data to be
//    interpolated.  The cubic spline has continous first and second
//    derivatives over the entire interval of interpolation.
//
//    For any point T in the interval T(IVAL), T(IVAL+1), the form of
//    the spline is
//
//      SPL(T) = A(IVAL)
//             + B(IVAL) * ( T - T(IVAL) )
//             + C(IVAL) * ( T - T(IVAL) )**2
//             + D(IVAL) * ( T - T(IVAL) )**3
//
//    If we assume that we know the values Y(*) and YPP(*), which represent
//    the values and second derivatives of the spline at each knot, then
//    the coefficients can be computed as:
//
//      A(IVAL) = Y(IVAL)
//      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
//        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
//      C(IVAL) = YPP(IVAL) / 2
//      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
//
//    Since the first derivative of the spline is
//
//      SPL'(T) =     B(IVAL)
//              + 2 * C(IVAL) * ( T - T(IVAL) )
//              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
//
//    the requirement that the first derivative be continuous at interior
//    knot I results in a total of N-2 equations, of the form:
//
//      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
//      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
//
//    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
//
//      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
//      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
//      + YPP(IVAL-1) * H(IVAL-1)
//      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
//      =
//      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
//      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
//
//    or
//
//      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
//      + YPP(IVAL) * H(IVAL)
//      =
//      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
//      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
//
//    Boundary conditions must be applied at the first and last knots.  
//    The resulting tridiagonal system can be solved for the YPP values.
//
//  Modified:
//
//    06 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.  N must be at least 2.
//    In the special case where N = 2 and IBCBEG = IBCEND = 0, the
//    spline will actually be linear.
//
//    Input, double T[N], the knot values, that is, the points were data is
//    specified.  The knot values should be distinct, and increasing.
//
//    Input, double Y[N], the data values to be interpolated.
//
//    Input, int IBCBEG, left boundary condition flag:
//      0: the cubic spline should be a quadratic over the first interval;
//      1: the first derivative at the left endpoint should be YBCBEG;
//      2: the second derivative at the left endpoint should be YBCBEG.
//
//    Input, double YBCBEG, the values to be used in the boundary
//    conditions if IBCBEG is equal to 1 or 2.
//
//    Input, int IBCEND, right boundary condition flag:
//      0: the cubic spline should be a quadratic over the last interval;
//      1: the first derivative at the right endpoint should be YBCEND;
//      2: the second derivative at the right endpoint should be YBCEND.
//
//    Input, double YBCEND, the values to be used in the boundary
//    conditions if IBCEND is equal to 1 or 2.
//
//    Output, double SPLINE_CUBIC_SET[N], the second derivatives of the cubic spline.
//
{
  double *a;
  double *b;
  int i;
  double *ypp;
//
//  Check.
//
  if ( n <= 1 )
  {
	/*
    cout << "\n";
    cout << "SPLINE_CUBIC_SET - Fatal error!\n";
    cout << "  The number of data points N must be at least 2.\n";
    cout << "  The input value is " << n << ".\n";
	*/
    return NULL;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    if ( t[i+1] <= t[i] )
    {
	  /*
      cout << "\n";
      cout << "SPLINE_CUBIC_SET - Fatal error!\n";
      cout << "  The knots must be strictly increasing, but\n";
      cout << "  T(" << i   << ") = " << t[i]   << "\n";
      cout << "  T(" << i+1 << ") = " << t[i+1] << "\n";
	  */
      return NULL;
    }
  }
  a = new double[3*n];
  b = new double[n];
//
//  Set up the first equation.
//
  if ( ibcbeg == 0 )
  {
    b[0] = 0.0;
    a[1+0*3] = 1.0;
    a[0+1*3] = -1.0;
  }
  else if ( ibcbeg == 1 )
  {
    b[0] = ( y[1] - y[0] ) / ( t[1] - t[0] ) - ybcbeg;
    a[1+0*3] = ( t[1] - t[0] ) / 3.0;
    a[0+1*3] = ( t[1] - t[0] ) / 6.0;
  }
  else if ( ibcbeg == 2 )
  {
    b[0] = ybcbeg;
    a[1+0*3] = 1.0;
    a[0+1*3] = 0.0;
  }
  else
  {
	/*
    cout << "\n";
    cout << "SPLINE_CUBIC_SET - Fatal error!\n";
    cout << "  IBCBEG must be 0, 1 or 2.\n";
    cout << "  The input value is " << ibcbeg << ".\n";
	*/
    delete [] a;
    delete [] b;
    return NULL;
  }
//
//  Set up the intermediate equations.
//
  for ( i = 1; i < n-1; i++ )
  {
    b[i] = ( y[i+1] - y[i] ) / ( t[i+1] - t[i] )
      - ( y[i] - y[i-1] ) / ( t[i] - t[i-1] );
    a[2+(i-1)*3] = ( t[i] - t[i-1] ) / 6.0;
    a[1+ i   *3] = ( t[i+1] - t[i-1] ) / 3.0;
    a[0+(i+1)*3] = ( t[i+1] - t[i] ) / 6.0;
  }
//
//  Set up the last equation.
//
  if ( ibcend == 0 )
  {
    b[n-1] = 0.0;
    a[2+(n-2)*3] = -1.0;
    a[1+(n-1)*3] = 1.0;
  }
  else if ( ibcend == 1 )
  {
    b[n-1] = ybcend - ( y[n-1] - y[n-2] ) / ( t[n-1] - t[n-2] );
    a[2+(n-2)*3] = ( t[n-1] - t[n-2] ) / 6.0;
    a[1+(n-1)*3] = ( t[n-1] - t[n-2] ) / 3.0;
  }
  else if ( ibcend == 2 )
  {
    b[n-1] = ybcend;
    a[2+(n-2)*3] = 0.0;
    a[1+(n-1)*3] = 1.0;
  }
  else
  {
	/*
    cout << "\n";
    cout << "SPLINE_CUBIC_SET - Fatal error!\n";
    cout << "  IBCEND must be 0, 1 or 2.\n";
    cout << "  The input value is " << ibcend << ".\n";
	*/
    delete [] a;
    delete [] b;
    return NULL;
  }
//
//  Solve the linear system.
//
  if ( n == 2 && ibcbeg == 0 && ibcend == 0 )
  {
    ypp = new double[2];

    ypp[0] = 0.0;
    ypp[1] = 0.0;
  }
  else
  {
    ypp = d3_np_fs ( n, a, b );

    if ( !ypp )
    {
	  /*
      cout << "\n";
      cout << "SPLINE_CUBIC_SET - Fatal error!\n";
      cout << "  The linear system could not be solved.\n";
	  */
      delete [] a;
      delete [] b;
      return NULL;
    }

  }

  delete [] a;
  delete [] b;
  return ypp;
}
//**********************************************************************

double spline_cubic_val ( int n, const double *t, double tval, const double *y, 
  const double *ypp, double *ypval, double *yppval )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
//
//  Discussion:
//
//    SPLINE_CUBIC_SET must have already been called to define the values of YPP.
//
//    For any point T in the interval T(IVAL), T(IVAL+1), the form of
//    the spline is
//
//      SPL(T) = A
//             + B * ( T - T(IVAL) )
//             + C * ( T - T(IVAL) )**2
//             + D * ( T - T(IVAL) )**3
//
//    Here:
//      A = Y(IVAL)
//      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
//        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
//      C = YPP(IVAL) / 2
//      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
//
//  Modified:
//
//    04 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int n, the number of knots.
//
//    Input, double Y[N], the data values at the knots.
//
//    Input, double T[N], the knot values.
//
//    Input, double TVAL, a point, typically between T[0] and T[N-1], at
//    which the spline is to be evalulated.  If TVAL lies outside
//    this range, extrapolation is used.
//
//    Input, double Y[N], the data values at the knots.
//
//    Input, double YPP[N], the second derivatives of the spline at
//    the knots.
//
//    Output, double *YPVAL, the derivative of the spline at TVAL.
//
//    Output, double *YPPVAL, the second derivative of the spline at TVAL.
//
//    Output, double SPLINE_VAL, the value of the spline at TVAL.
//
{
  double dt;
  double h;
  int i;
  int ival;
  double yval;
//
//  Determine the interval [ T(I), T(I+1) ] that contains TVAL.
//  Values below T[0] or above T[N-1] use extrapolation.
//
  ival = n - 2;

  for ( i = 0; i < n-1; i++ )
  {
    if ( tval < t[i+1] )
    {
      ival = i;
      break;
    }
  }
//
//  In the interval I, the polynomial is in terms of a normalized
//  coordinate between 0 and 1.
//
  dt = tval - t[ival];
  h = t[ival+1] - t[ival];

  yval = y[ival]
    + dt * ( ( y[ival+1] - y[ival] ) / h
           - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
    + dt * ( 0.5 * ypp[ival]
    + dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * h ) ) ) );

  *ypval = ( y[ival+1] - y[ival] ) / h
    - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
    + dt * ( ypp[ival]
    + dt * ( 0.5 * ( ypp[ival+1] - ypp[ival] ) / h ) );

  *yppval = ypp[ival] + dt * ( ypp[ival+1] - ypp[ival] ) / h;

  return yval;
}
//**********************************************************************

void spline_cubic_val2 ( int n, const double *t, double tval, int *left, const double *y, 
  const double *ypp, double *yval, double *ypval, double *yppval )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_CUBIC_VAL2 evaluates a piecewise cubic spline at a point.
//
//  Discussion:
//
//    This routine is a modification of SPLINE_CUBIC_VAL; it allows the
//    user to speed up the code by suggesting the appropriate T interval
//    to search first.
//
//    SPLINE_CUBIC_SET must have already been called to define the
//    values of YPP.
//
//    In the LEFT interval, let RIGHT = LEFT+1.  The form of the spline is
//
//    SPL(T) =
//      A
//    + B * ( T - T[LEFT] )
//    + C * ( T - T[LEFT] )**2
//    + D * ( T - T[LEFT] )**3
//
//    Here:
//      A = Y[LEFT]
//      B = ( Y[RIGHT] - Y[LEFT] ) / ( T[RIGHT] - T[LEFT] )
//        - ( YPP[RIGHT] + 2 * YPP[LEFT] ) * ( T[RIGHT] - T[LEFT] ) / 6
//      C = YPP[LEFT] / 2
//      D = ( YPP[RIGHT] - YPP[LEFT] ) / ( 6 * ( T[RIGHT] - T[LEFT] ) )
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of knots.
//
//    Input, double T[N], the knot values.
//
//    Input, double TVAL, a point, typically between T[0] and T[N-1], at
//    which the spline is to be evalulated.  If TVAL lies outside
//    this range, extrapolation is used.
//
//    Input/output, int *LEFT, the suggested T interval to search.
//    LEFT should be between 1 and N-1.  If LEFT is not in this range,
//    then its value will be ignored.  On output, LEFT is set to the
//    actual interval in which TVAL lies.
//
//    Input, double Y[N], the data values at the knots.
//
//    Input, double YPP[N], the second derivatives of the spline at
//    the knots.
//
//    Output, double *YVAL, *YPVAL, *YPPVAL, the value of the spline, and
//    its first two derivatives at TVAL.
//
{
  double dt;
  double h;
  int right;
//
//  Determine the interval [T[LEFT], T[RIGHT]] that contains TVAL.  
//  
//  What you want from DVEC_BRACKET3 is that TVAL is to be computed
//  by the data in interval [T[LEFT-1], T[RIGHT-1]].  
//
  dvec_bracket3 ( n, t, tval, left );
//
// In the interval LEFT, the polynomial is in terms of a normalized
// coordinate  ( DT / H ) between 0 and 1.
//
  right = *left + 1;

  dt = tval - t[*left-1];
  h = t[right-1] - t[*left-1];

  *yval = y[*left-1]
     + dt * ( ( y[right-1] - y[*left-1] ) / h
        - ( ypp[right-1] / 6.0 + ypp[*left-1] / 3.0 ) * h
     + dt * ( 0.5 * ypp[*left-1]
     + dt * ( ( ypp[right-1] - ypp[*left-1] ) / ( 6.0 * h ) ) ) );

  *ypval = ( y[right-1] - y[*left-1] ) / h
     - ( ypp[right-1] / 6.0 + ypp[*left-1] / 3.0 ) * h
     + dt * ( ypp[*left-1]
     + dt * ( 0.5 * ( ypp[right-1] - ypp[*left-1] ) / h ) );

  *yppval = ypp[*left-1] + dt * ( ypp[right-1] - ypp[*left-1] ) / h;

  return;
}