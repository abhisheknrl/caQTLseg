#ifndef SICS_H
#define SICS_H

#include "types/matrix.h"

//
// 1-sample solvers
//

// Projected-gradient solver
bool SICS_pgd1(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int it);

// Duchi solver
bool SICS_pgd2(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int it);

// Heavy ball solver
bool SICS_pgd3(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int &it, double &magicbeta);

// Heavy ball contination solver
bool SICS_pgd4(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int &it, double &magicbeta);

// ADM
bool SICS_adm(CMatrix<double> &X, const double lambda, CMatrix<double> &Xi, const double tol_obj, const double tol_gap, const double tol_kkt, const double tol_steplen, const int tol_it, double &obj_this, double &duality_gap, double &dKKT_max, const bool &bSkipShrink, int it);

//
// 2-sample solvers
//

// Projected-gradient 
bool SICS_pgd1_paired(CMatrix<double> &X0, CMatrix<double> &X1, const double lambda, CMatrix<double> &Xi0, CMatrix<double> &Xi1, 
	const double tol_steplen, const double tol_gap, const double tol_kkt, const int tol_it, double &obj_this, double &duality_gap, double &dU_max, int &it);

#endif