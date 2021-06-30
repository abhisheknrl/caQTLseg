//
// matrix algebra helper routines (most of which wrap an mkl routine)
//
// Björn Nilsson, 2012
//

#ifndef MKLHELPERS_H
#define MKLHELPERS_H

#include "types/matrix.h"
#include "types/vector.h" 

bool mkl_GetAB(CMatrix<double> &A, CMatrix<double> &B, CMatrix<double> &m_out);
bool mkl_GetAAt(CMatrix<double> &A, CMatrix<double> &m_out);
bool mkl_GetAAt(CMatrix<float> &A, CMatrix<float> &m_out);

bool mkl_GetAtA(CMatrix<double> &A, CMatrix<double> &m_out);
bool mkl_GetTriangularAAt(CMatrix<double> &A, const char *uplo, CMatrix<double> &m_out);
bool mkl_GetTriangularAtA(CMatrix<double> &A, const char *uplo, CMatrix<double> &m_out);
bool mkl_GetCholesky(CMatrix<double> &A, const char *uplo, bool bZeroFill);
bool mkl_GetTriangularInverse(CMatrix<double> &A, const char *uplo, bool bZeroFill);
bool mkl_GetEigenDecomposition(CMatrix<double> &A, CVector<double> &eA, const char *uplo);
double mkl_GetLogdet(const CMatrix<double> &L);
void printmatrixd(const CMatrix<double> &m);
void printmatrixs(const CMatrix<float> &m);

#endif
