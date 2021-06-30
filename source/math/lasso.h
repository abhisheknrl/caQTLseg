/*******************************************************************************
 *
 *   lasso.cpp -  Solvers for l1-regularized least squares
 *
 *   Björn Nilsson, 2011
 *
 */

#ifndef LASSO_H
#define LASSO_H

#include "types/vector.h"
#include "types/matrix.h"

void printmatrix_double(const CMatrix<double> &m);
void printvector_double(const CVector<double> &v);

// constant lambda
bool LASSO_CCD_At(const CMatrix<double> &At, const CVector<double> &b, const double lambda, CVector<double> &x, const int nVerbosity);
bool LASSO_CCD_AtA(const CMatrix<double> &AtA, const CVector<double> &Atb, const double lambda, CVector<double> &x, const int nVerbosity);

// vectorized lambda
bool LASSO_CCD_At(const CMatrix<double> &At, const CVector<double> &b, const CVector<double> &vLambda, CVector<double> &x, const int nVerbosity);
bool LASSO_CCD_AtA(const CMatrix<double> &AtA, const CVector<double> &Atb, const CVector<double> &vLambda, CVector<double> &x, const int nVerbosity);

#endif