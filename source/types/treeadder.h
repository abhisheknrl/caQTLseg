//
//	treeadder.h --	Template class for performing serial additions, as defined 
//					by the + operator, in a tree structure. Useful when the cost 
//					for adding large (complex) T is computationally expensive 
//					compared with adding small (simple) T.
//
//	Björn Nilsson, 2009
//

template<class T> // T must have + and = operators, and an identity element
class CRevolverAdder
{
protected:
	CVector<T> m_vX;	// Level sum
	int m_nN;			// Number of terms added
public:
	CRevolverAdder() { m_nN= 0; }

	void Init(int levels) // call after ctor before other functions
	{
		m_vX.SetSize(levels);
		m_nN= 0;
	}

	void SetSum(T x) 
	{
		m_vX.SetAt(0, x);
		m_nN= 1;
	}

	inline void Add(const T &x) 
	{ 
		int i= m_nN % m_vX.GetSize();
		if (i>=m_nN)
			m_vX.SetAt(i, x);
		else
			m_vX.SetAt(i, x + m_vX[i]);
		m_nN++;
	}

	T GetSum() const
	{
		ASSERT(m_nN>0); // Cumulator not initialized, call SetSum or Add first
		int n;
		if (m_nN<m_vX.GetSize())
			n= m_nN;
		else
			n= m_vX.GetSize();

		T s= m_vX[0];
		for (int i=1;i<n;i++)
			s= m_vX[i] + s;

		return s;
	}
};
/*
template<class T> // T must have + and = operators, and an identity element
class CTreeAdder
{
protected:
	CVector<T> m_vX;	// Level sum
	CVector<bool> m_vP;	// Level parity (true <--> vX valid at position)
	int m_nN;			// Number of terms added
public:
	CTreeAdder() { m_nN= 0; }

	void Init(int levels) // call after ctor before other functions
	{
		ASSERT(levels>1);
		m_vX.SetSize(levels);
		m_vP.SetSize(levels);
		m_vP.Fill(false);
		m_nN= 0;
	}

	void SetSum(T x) 
	{
		m_vP.Fill(false);
		m_vP.SetAt(0, true);
		m_vX.SetAt(0, x);
		m_nN= 1;
	}

	inline void Add(T x) 
	{ 
		int level= 0;
		while (m_vP[level])
		{
			printf("level %d\n", level);
			x= m_vX[level] + x;
			g_mul++;
			m_vP.SetAt(level, false);
			level++;
		}
		m_vX.SetAt(level, x);
		m_vP.SetAt(level, true);
		m_nN++;
	}

	T GetSum() const
	{
		ASSERT(m_nN>0); // Cumulator not initialized, call SetSum or Add first
		int levels=1;
		int n= m_nN;
		while ((n= n >> 1)>0)
			levels++;

		int i=0;
		while (i<levels && !m_vP[i]) 
			i++;
		T s= m_vX[i];
		i++;
		for (;i<levels;i++)
			if (m_vP[i])
				s= s + m_vX[i];
		return s;
	}
};
*/