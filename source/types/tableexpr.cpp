#include "tableexpr.h"
#include "system/filehelpers.h"
#include "types/stringfun.h"
#include "tablefun.h"
/*
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

class token
{
public:
	int t;
	int i; // integer value, string index, or variable index
	double x;// double value for constants
	token() { t= TOKEN_ERROR; i= 0; }
	token(int t_init, int i_init)
	{
		t= t_init;
		i= i_init;
		x= 0;
	}
	token(int t_init, int i_init, double x_init)
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
*/

const char *expr_ScanIdentifier(const char *ach, CString &sID)
{
	if (!ach)
		return 0;
	// Parses quoted strings and c-style identifier
	ach= Scan_SkipRichSpace(ach);

	const char *ach0= ach;
	if (*ach=='@' || *ach=='$') // Inital @ or $ allowed, @ --> tokenize as col/row index even if numeric, $ --> tokenize as col/row name even if numeric
		ach++;
	if (*ach=='!') // name itself can start with ! (motivated by GEO matrix file format, where row names start with !)
		ach++;
	while (*ach && ((*ach>='0' && *ach<='9') || 
					(*ach>='A' && *ach<='Z') ||
					(*ach>='a' && *ach<='z') ||
					(*ach=='_' || *ach=='.' || *ach=='#' || *ach==':')))
		ach++;
	sID= CString(ach0, ach-ach0);
	return ach;	
}

const char *expr_ScanQuotedIdentifier(const char *ach, CString &sID)
{
	if (ach && ach[0]=='\'')
	{
		ach++;
		const char *ach0= ach;
		while (*ach && *ach!='\'')
			ach++;
		if (*ach)
		{
			sID= CString(ach0, ach-ach0);
			return ach+1;
		}
	}
	return 0;
}

const char *expr_ScanPossiblyQuotedIdentifier(const char *ach, CString &sID)
{
	if (ach && ach[0]=='\'')
		return expr_ScanQuotedIdentifier(ach, sID);
	else
		return expr_ScanIdentifier(ach, sID);
}

bool IsRightAssociative(int t)
{
	// true if t is a right-associative operator, false otherwise
	return t==TOKEN_EXP;
	return false;
}

int GetOperatorPrecedence(int t)
{
	/*
	operator_Priority["="] = 0; // equal
	
	operator_Priority["&&"] = 1;// and
	operator_Priority["||"] = 1;// or
	
	operator_Priority["!"] = 2;// not

	operator_Priority["=="] = 3;// equality check
	operator_Priority["!="] = 3;// not equal
	operator_Priority["<" ] = 3;// lower
	operator_Priority[">" ] = 3;// higher
	operator_Priority["<="] = 3;// lower or equal
	operator_Priority[">="] = 3;// higher or equal
	
	operator_Priority["+"] = 4; // addition
	operator_Priority["-"] = 4; // subtraction

	operator_Priority["*"] = 5; // multiply
	operator_Priority["/"] = 5; // division
	operator_Priority["%"] = 5; // mod

	operator_Priority["~"] = 6; // unary -
	operator_Priority["^"] = 6; // exponantion
	
	operator_Priority["("] = -1; // parantheses
	operator_Priority[")"] = -2; // prantheses}
	*/

	switch (t)
	{
		case TOKEN_PAROPEN :
			return -1;
		case TOKEN_PARCLOSE : 
			return -2;
		case TOKEN_LOGICAND :
		case TOKEN_LOGICOR : 
			return 1;
		case TOKEN_LOGICNOT : 
			return 2;
		case TOKEN_EQ :
		case TOKEN_NEQ : 
		case TOKEN_LT :
		case TOKEN_LET :
		case TOKEN_GT :
		case TOKEN_GET :
		case TOKEN_LIKE :
		case TOKEN_HAS :
		case TOKEN_INSET : 
			return 3;
		case TOKEN_ADD :
		case TOKEN_SUB : 
			return 4;
		case TOKEN_MUL :
		case TOKEN_DIV : 
			return 5;
		case TOKEN_EXP : 
		case TOKEN_UNARYMINUS : 
			return 6;
		default: 
			ReportError("Unknown operator precedence", ::Format("%d", int(t)));
			return -1000;
	}
}

bool ToRPN(CVector<expr_token> &vT, CVector<expr_token> &vT_out)
{
	// Implements Dijkstraa's shunting yard algorithm -- from Wikipedia
	CVector<expr_token> vS;
	vS.SetSize(100000);
	int s= 0;

	for (int i=0;i<vT.GetSize();i++)
	{
		expr_token o1= vT[i];
		if (o1.t>TOKEN_OPERAND)
			vT_out.Add(o1); // operand type
		else if (o1.t==TOKEN_PAROPEN) // If the token is a left parenthesis, then push it onto the stack.
			vS.SetAt(s++, o1);
		else if (o1.t==TOKEN_PARCLOSE) // Right paranthesis
		{
			// Until the token at the top of the stack is a left parenthesis, pop operators off the stack onto the output queue.
			while (s>0 && vS[s-1].t!=TOKEN_PAROPEN)
				vT_out.Add(vS[--s]);

			// If the stack runs out without finding a left parenthesis, then there are mismatched parentheses.
			if (s==0)
			{
				ReportError("Missing left paranthesis", 0);
				return false;
			}

			// Pop the left parenthesis from the stack, but not onto the output queue.
			s--;
		}
		else // 
		{
			int p1= GetOperatorPrecedence(o1.t);
			if (IsRightAssociative(o1.t))
				while (s>0 && p1<GetOperatorPrecedence(vS[s-1].t))
					vT_out.Add(vS[--s]);
			else
				while (s>0 && p1<=GetOperatorPrecedence(vS[s-1].t))
					vT_out.Add(vS[--s]);
			vS.SetAt(s++, o1);
		}
	}

	// Pop entire stack to output
	while (s>0)
		vT_out.Add(vS[--s]);

	return true;
}

enum {
	ERR_NONE= 0,
	ERR_STACKUNDERFLOW,
	ERR_INCOMPATIBLETYPES,
	ERR_ILLEGALOP
};

inline void PushInt(expr_operand *pS, int &s, int i)
{
	pS[s].val_int= i;
	pS[s].t= TOKEN_INT;
	s++;
}

inline void PushStr(expr_operand *pS, int &s, const CString &str)
{
	pS[s].val_str= str;
	pS[s].t= TOKEN_STR;
	s++;
}

inline void PushDbl(expr_operand *pS, int &s, double x)
{
	pS[s].val_dbl= x;
	pS[s].t= TOKEN_DBL;
	s++;
}

inline void PushVar(expr_operand *pS, int &s, const CString &str)
{
	pS[s].val_str= str;

	const char *ach= str;
	const char *ach_end= ach+str.GetLength();
	if (Table_CheckFloatFormat(ach)==ach_end)
	{
		pS[s].t= TOKEN_STR_AND_DBL;
		pS[s].val_dbl= atof(ach);
	}
	else
		pS[s].t= TOKEN_STR;

	s++;
}

inline bool PopBool(expr_operand *pS, int &s)
{
	// note: recasts ints/dbl to bool as needed
	return pS[--s].GetBool();
}

inline double PopDbl(expr_operand *pS, int &s)
{
	// note: recasts ints to dbl if needed
	double x= pS[--s].GetDbl();
	return x;
}

CString PopStr(expr_operand *pS, int &s)
{
	return pS[--s].val_str;
}

inline void Swap(expr_operand *pS, int s)
{
	expr_operand tmp= pS[s-2];
	pS[s-2]= pS[s-1];
	pS[s-1]= tmp;
}

int DoBinary(int op, expr_operand *pS, int &s)
{
	if (s<2)
		return ERR_STACKUNDERFLOW;

	if (pS[s-1].HasNum() && pS[s-2].HasNum())
	{
		switch(op)
		{
			case TOKEN_LOGICAND : PushInt(pS, s, PopBool(pS, s) & PopBool(pS, s)); break;
			case TOKEN_LOGICOR  : PushInt(pS, s, PopBool(pS, s) | PopBool(pS, s)); break;
			case TOKEN_EQ  : PushInt(pS, s, PopDbl(pS, s) == PopDbl(pS, s)); break;
			case TOKEN_NEQ : PushInt(pS, s, PopDbl(pS, s) != PopDbl(pS, s)); break;
			case TOKEN_LT  : PushInt(pS, s, PopDbl(pS, s) >  PopDbl(pS, s)); break;
			case TOKEN_LET : PushInt(pS, s, PopDbl(pS, s) >= PopDbl(pS, s)); break;
			case TOKEN_GT  : PushInt(pS, s, PopDbl(pS, s) <  PopDbl(pS, s)); break;
			case TOKEN_GET : PushInt(pS, s, PopDbl(pS, s) <= PopDbl(pS, s)); break;
			case TOKEN_ADD : PushDbl(pS, s, PopDbl(pS, s) + PopDbl(pS, s)); break;
			case TOKEN_MUL : PushDbl(pS, s, PopDbl(pS, s) * PopDbl(pS, s)); break;
			case TOKEN_SUB : Swap(pS, s); PushDbl(pS, s, PopDbl(pS, s) - PopDbl(pS, s)); break;
			case TOKEN_DIV : Swap(pS, s); PushDbl(pS, s, PopDbl(pS, s) / PopDbl(pS, s)); break;
			default: 
				return ERR_ILLEGALOP;
		}
	}
	else if (pS[s-1].HasStr() && pS[s-2].HasStr())
	{
		// note: items are popped in reverse order
		CString s1= PopStr(pS, s);
		CString s0= PopStr(pS, s);
		switch(op)
		{
			case TOKEN_EQ    : PushInt(pS, s, s0 == s1); break;
			case TOKEN_NEQ   : PushInt(pS, s, s0 != s1); break;
			case TOKEN_LT    : PushInt(pS, s, s0 <  s1); break;
			case TOKEN_LET   : PushInt(pS, s, s0 <= s1); break;
			case TOKEN_GT    : PushInt(pS, s, s0 >  s1); break;
			case TOKEN_GET   : PushInt(pS, s, s0 >= s1); break;
			case TOKEN_ADD   : PushStr(pS, s, s0 +  s1); break;
			case TOKEN_LIKE  : PushInt(pS, s, strmatch(s0, s1)); break;
			case TOKEN_HAS   : PushInt(pS, s, int(s0.Find(s1)>=0)); break;
			case TOKEN_INSET : PushInt(pS, s, strinset(s0, s1, '|')>=0); break;
			default: 
				return ERR_ILLEGALOP;
		}
	}
	else 
	{
		 printf("type 1: '%s'\n", (const char *)pS[s-1].val_str);
		 printf("type 2: '%s'\n", (const char *)pS[s-2].val_str);
		return ERR_INCOMPATIBLETYPES;
	}

	return ERR_NONE;
}

CString GetErrString(int nErr)
{
	switch(nErr)
	{
		case ERR_ILLEGALOP : return "Illegal operation";
		case ERR_NONE : return "Success";
		case ERR_INCOMPATIBLETYPES : return "Incompatible types";
		case ERR_STACKUNDERFLOW : return "Operator lacks operand";
		default : return "Unknown error";
	}
}

int expr_Evaluate(const CTableIndex &ti, const CVector<expr_token> &vT, const CVector<CString> &vTstr, bool bRowID, int col_or_row, expr_operand &res, expr_operand *pStack)
{
	// returns type-token of result (if successful, value in corresponding val_??? variable), TOKEN_ERROR otherwise

	int s= 0; // stack pointer

	int nErr= 0;
	for (int i=0;nErr==0 && i<vT.GetSize();i++)
	{
		int op= vT[i].t;
		switch (op)
		{
			case TOKEN_INT : 
				PushInt(pStack, s, vT[i].i); break;
			case TOKEN_DBL : 
				PushDbl(pStack, s, vT[i].x); break;
			case TOKEN_STR : 
				PushStr(pStack, s, vTstr[vT[i].i]); break;
			case TOKEN_VAR :
				{
					// robust
					int r= bRowID ? vT[i].i : col_or_row;
					int c= bRowID ? col_or_row : vT[i].i;
					if (r>=0 && r<ti.GetRowCount() &&
						c>=0 && c<ti.GetColCount(r))
						PushVar(pStack, s, Table_GetAt(ti, r, c));
					else
						PushVar(pStack, s, "");
				}
				break;
			case TOKEN_LOGICNOT : 
				if (s==0)
				{
					ReportError("NOT without operand", 0);
					return TOKEN_ERROR;
				}
				if (!pStack[s-1].HasNum())
				{
					ReportError("NOT requires numeric operand", 0);
					return TOKEN_ERROR;
				}
				PushInt(pStack, s, int(!PopBool(pStack, s)));
				break;
			case TOKEN_INSET :
			case TOKEN_LIKE :
			case TOKEN_HAS : 
			case TOKEN_LOGICAND :
			case TOKEN_LOGICOR : 
			case TOKEN_EQ :
			case TOKEN_NEQ : 
			case TOKEN_LT :
			case TOKEN_LET :
			case TOKEN_GT :
			case TOKEN_GET :
			case TOKEN_ADD :
			case TOKEN_SUB : 
			case TOKEN_MUL :
			case TOKEN_DIV : 
				nErr= DoBinary(vT[i].t, pStack, s); 
				if (nErr!=ERR_NONE)
				{
					ReportError(GetErrString(nErr), 0);
					return TOKEN_ERROR;
				}
				break;
			case TOKEN_UNARYMINUS : 
				if (s==0)
				{
					ReportError("Unary minus without operand", 0); 
					return TOKEN_ERROR;
				}
				if (!pStack[s-1].HasNum())
				{
					ReportError("Unary minus requires numeric operand", 0);
					return TOKEN_ERROR;
				}
				pStack[s-1].val_int= -pStack[s-1].val_int;
				pStack[s-1].val_dbl= -pStack[s-1].val_dbl;
				break;
			case TOKEN_EXP : 
				ReportError("Operator '^' not implemented", 0); return TOKEN_ERROR;

			default: 
				ReportError("Token not implemented", ::Format("%d", int(op)));
				return TOKEN_ERROR;
		}
	}

	if (s==0)
	{
		ReportError("Operator without operand", 0);
		return TOKEN_ERROR;
	}
	if (s>1)
	{
		ReportError("Operand without operator", 0);
		return TOKEN_ERROR;
	}

	res= pStack[0];
	return pStack[0].t;
}

bool expr_Tokenize(const CTableIndex &ti, const char *ach, bool bRowID, CVector<expr_token> &vT, CVector<CString> &vTstr)
{
	// Tokenize expression string + convert to rpn
	while (*ach)
	{
		ach= Scan_SkipRichSpace(ach);
		const char *ach_err= ach;

		if (ach[0]=='=')
		{
			if (ach[1]=='=')
				vT.Add(expr_token(TOKEN_EQ, 0));
			else
			{
				ReportError("Use == instead of = to check equality", ach_err);
				return false;
			}
			ach+= 2;
		}
		if (ach[0]=='~' && ach[1]=='=')
		{
			vT.Add(expr_token(TOKEN_LIKE, 0));
			ach+=2;
		}
		else if (ach[0]=='!' && ach[1]=='=')
		{
			vT.Add(expr_token(TOKEN_NEQ, 0));
			ach+=2;
		}
		else if (ach[0]=='&' && ach[1]=='&')
		{
			vT.Add(expr_token(TOKEN_LOGICAND, 0));
			ach+=2;
		}
		else if (ach[0]=='|' && ach[1]=='|')
		{
			vT.Add(expr_token(TOKEN_LOGICOR, 0));
			ach+=2;
		}
		else if (ach[0]=='\'')
		{
			ach++;
			const char *ach0= ach;
			while (*ach && *ach!='\'')
				ach++;
			if (*ach==0)
			{
				ReportError("End quote expected", ach_err);
				return false;
			}
			vT.Add(expr_token(TOKEN_STR, vTstr.GetSize()));
			vTstr.Add(CString(ach0, ach-ach0));
			ach++;
		}
		else if (ach[0]=='(')
		{
			vT.Add(expr_token(TOKEN_PAROPEN, 0));
			ach++;
		}
		else if (ach[0]==')')
		{
			vT.Add(expr_token(TOKEN_PARCLOSE, 0));
			ach++;
		}
		else if (ach[0]=='-')
		{
			// Determine if minus or unary minus by looking at previous token
			if (vT.GetSize()>0 && (vT[vT.GetSize()-1].t<TOKEN_OPERAND || vT[vT.GetSize()-1].t==TOKEN_PARCLOSE))
				vT.Add(expr_token(TOKEN_UNARYMINUS, 0));
			else
				vT.Add(expr_token(TOKEN_SUB, 0));
			ach++;
		}
		else if (ach[0]=='+')
		{
			vT.Add(expr_token(TOKEN_ADD, 0));
			ach++;
		}
		else if (ach[0]=='*')
		{
			vT.Add(expr_token(TOKEN_MUL, 0));
			ach++;
		}
		else if (ach[0]=='/')
		{
			vT.Add(expr_token(TOKEN_DIV, 0));
			ach++;
		}
		else if (ach[0]=='^')
		{
			vT.Add(expr_token(TOKEN_EXP, 0));
			ach++;
		}
		else if (ach[0]=='<')
		{
			if (ach[1]=='=')
			{
				vT.Add(expr_token(TOKEN_LET, 0));
				ach+=2;
			}
			else
			{
				vT.Add(expr_token(TOKEN_LT, 0));
				ach++;
			}
		}
		else if (ach[0]=='>')
		{
			if (ach[1]=='=')
			{
				vT.Add(expr_token(TOKEN_GET, 0));
				ach+=2;
			}
			else
			{
				vT.Add(expr_token(TOKEN_GT, 0));
				ach++;
			}
		}
		else
		{
			// Must be a string, unsigned int, or float number
			CString s;
			ach= expr_ScanIdentifier(ach, s);
			if (!ach)
			{
				ReportError("Illegal identifer format", ach_err);
				return false;
			}
			const char *ach_s= s;
			const char *ach_s_end= ach_s+s.GetLength();

			if (Table_CheckUnsignedIntFormat(ach_s)==ach_s_end)
				vT.Add(expr_token(TOKEN_INT, atoi(s)));
			else if (Table_CheckFloatFormat(ach_s)==ach_s_end)
			{
				double x= atof(s);
				vT.Add(expr_token(TOKEN_DBL, 0, x));
			}
			else if (s=="LIKE")
				vT.Add(expr_token(TOKEN_LIKE, 0));
			else if (s=="HAS")
				vT.Add(expr_token(TOKEN_HAS, 0));
			else if (s=="IN")
				vT.Add(expr_token(TOKEN_INSET, 0));
			else if (s=="AND")
				vT.Add(expr_token(TOKEN_LOGICAND, 0));
			else if (s=="OR")
				vT.Add(expr_token(TOKEN_LOGICOR, 0));
			else if (s=="NOT")
				vT.Add(expr_token(TOKEN_LOGICNOT, 0));
			else if (s=="TRUE")
				vT.Add(expr_token(TOKEN_INT, 1));
			else if (s=="FALSE")
				vT.Add(expr_token(TOKEN_INT, 0));
			else if (s[0]=='@')
			{
				// User wants this to be coordinate (even if numeric)
				s= s.Mid(1, s.GetLength());
				const char *ach_s= s;
				const char *ach_s_end= ach_s+s.GetLength();
				if (Table_CheckUnsignedIntFormat(ach_s)==ach_s_end)
					vT.Add(expr_token(TOKEN_VAR, atoi(s)));
				else
				{
					ReportError("Illegal identifier", CString("@")+s);
					return 0;
				}
			}
			else if (s[0]=='$')
			{
				// User wants this to be a variable name (even if numeric)
				s= s.Mid(1, s.GetLength());
				int i= bRowID ? Table_FindRow(ti, s) : Table_FindCol(ti, s);
				if (i<0)
				{
					ReportError(bRowID ? "Row not found" : "Column not found", s);
					return false;
				}
				vT.Add(expr_token(TOKEN_VAR, i));
			}			
			else
			{
				// Variable specified by unquoted row/col name
				int i= bRowID ? Table_ResolveRow(ti, s) : Table_ResolveCol(ti, s);
				if (i<0)
				{
					ReportError(bRowID ? "Row not found" : "Column not found", s);
					return false;
				}
				vT.Add(expr_token(TOKEN_VAR, i));
			}
		}
	}

	// Convert from infix form to postfix/RPN form
	CVector<expr_token> vRPN;
	if (!ToRPN(vT, vRPN))
		return false;

	vT= vRPN;
	return true;
}

