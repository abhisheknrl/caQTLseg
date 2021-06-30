/*******************************************************************************
 *
 *   cdf.h -- Implements the Nilsson-Fioretos algorithm
 *			  for computing reduced representations of 
 *			  an empirical distribution function of any
 *			  continuous-valued random variable.
 *
 *   Björn Nilsson, Medical Mathematics, 2006
 */

#ifndef CDF_H
#define CDF_H

#include "types/vector.h"

class CReducedCDF
{
protected:
	int	m_nN;
	int m_nDF;
	double m_dEpsilon;
	CVector<double> m_vT;

	void Uniquify(CVector<double> &vT);
	bool CreateErrorFile(CVector<double> &T);
public:
	// Initialization
	bool Initialize(CVector<double> &vT, double epsilon, int nDF);
	bool IsValid() const { return m_vT.GetSize()>0 && m_nN>0 && m_dEpsilon>0; }

	// Returns size of original sample
	int	GetN() const	{ ASSERT(IsValid()); return m_nN;  }
	int GetDF()	const	{ ASSERT(IsValid()); return m_nDF; }

	// Returns the relative error bound epsilon
	double GetEpsilon() { ASSERT(IsValid()); return m_dEpsilon; }

	// Returns reconstructed distribution function value 
	// The algorithm guarantees (1-\epsilon)F < Fstar < (1+\epsilon)F.
	double GetFstar(const double T, bool bInterpolate=true) const;

	// Serialization
	bool SaveToFile(const char *aFilename);
	bool LoadFromFile(const char *aFilename);
	bool IsEqual(const CReducedCDF &other);
};

#endif