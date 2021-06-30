/*******************************************************************************
 *
 *   tableexpr.h -- Computation of expressions related to tables
 *
 *   Björn Nilsson, 2008
 *
 */

#ifndef TABLEEXPR_H
#define TABLEEXPR_H

#include "types/string.h"
#include "types/vector.h"
#include "math/statfun.h"
#include "tableindex.h"

enum 
{
	// Operator tokens
	TOKEN_ERROR= -1,
	TOKEN_OPERATOR,
	TOKEN_PAROPEN,
	TOKEN_PARCLOSE,
	TOKEN_ADD, TOKEN_SUB,
	TOKEN_MUL, TOKEN_DIV,
	TOKEN_EQ, TOKEN_NEQ, TOKEN_GT, TOKEN_GET, TOKEN_LT, TOKEN_LET, TOKEN_LIKE, TOKEN_HAS, TOKEN_INSET,
	TOKEN_LOGICAND, TOKEN_LOGICOR, TOKEN_LOGICNOT,
	TOKEN_EXP, TOKEN_UNARYMINUS,

	// Operand tokens
	TOKEN_OPERAND,
	TOKEN_STR,
	TOKEN_INT,
	TOKEN_DBL,
	TOKEN_VAR,
	TOKEN_STR_AND_DBL
};

class expr_token
{
public:
	int t;
	int i; // integer value, string index, or variable index
	double x;// double value for constants
	expr_token() { t= TOKEN_ERROR; i= 0; }
	expr_token(int t_init, int i_init)
	{
		t= t_init;
		i= i_init;
		x= 0;
	}
	expr_token(int t_init, int i_init, double x_init)
	{
		t= t_init;
		i= i_init;
		x= x_init;
	}
};

class expr_operand
{
public:
	int t;
	CString val_str;
	double val_dbl;
	int val_int;

	// Type checking
	bool HasStr() { return t==TOKEN_STR || t==TOKEN_VAR || t==TOKEN_STR_AND_DBL; }
	bool HasDbl() { return t!=TOKEN_STR || t==TOKEN_STR_AND_DBL; }
	bool HasInt() { return t==TOKEN_INT; }
	bool HasNum() { return HasDbl() | HasInt(); }

	// Type casting
	double GetDbl()
	{
		if (HasDbl())
			return (double) val_dbl;
		else if (HasInt())
			return double(val_int);
		else
			return Stat_GetNaN_double();
	}
	bool GetBool()
	{
		if (HasInt())
			return val_int != 0;
		else
			return val_dbl != 0.0f;
	}
};

bool expr_Tokenize(const CTableIndex &ti, const char *ach, bool bRowId, CVector<expr_token> &vT, CVector<CString> &vTstr);
int expr_Evaluate(const CTableIndex &ti, const CVector<expr_token> &vT, const CVector<CString> &vTstr, bool bRowID, int col_or_row, expr_operand &res, expr_operand *pStack);
const char *expr_ScanIdentifier(const char *ach, CString &sID);
const char *expr_ScanQuotedIdentifier(const char *ach, CString &sID);
const char *expr_ScanPossiblyQuotedIdentifier(const char *ach, CString &sID);

#endif