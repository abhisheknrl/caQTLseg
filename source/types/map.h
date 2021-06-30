/***************************************************************************
 *
 *   map.h -- Japanized MFC Style Subset Map Classes
 *
 *   Erik Persson, Sketchware, 2000
 */

#ifndef TYPES_MAP_H
#define TYPES_MAP_H

#define MFC_CLASSES

#ifdef MFC_CLASSES
#include <afx.h>
#include <afxcoll.h>
#else
#include <types/string.h>

// Map iterator
struct __POSITION { };
typedef __POSITION* POSITION;

#define BEFORE_START_POSITION ((POSITION)-1L)

#endif // MFC_CLASSES

/****************************** CMemBlock **********************************/

struct CMemBlock     // warning variable length structure
{
	CMemBlock* pNext;
#if (_AFX_PACKING >= 8)
	DWORD dwReserved[1];    // align on 8 byte boundary
#endif
	// BYTE data[maxNum*elementSize];

	void* data() { return this+1; }

	static CMemBlock* Create(CMemBlock*& head, UINT nMax, UINT cbElement);
			// like 'calloc' but no zero fill
			// may throw memory exceptions

	void FreeDataChain();       // free this one and links
};

/**************************** CString Helpers *******************************/

template<class T>
inline void ConstructElement(T *pElement)
{
  // First set to zero, then call the constructor
  memset((void*)pElement, 0, sizeof(T));
	::new((void*)pElement) T;
}

template<class T>
inline void DestructElement(T *pElement)
{
	pElement->~T();
}

template<class T>
inline void CopyElement(T *pSrc, T *pDest)
{
	*pSrc= *pDest;
}

/***************************************************************************/

template<class VALUE, class ARG_VALUE>
class CStringMap
{
public:
	typedef CString KEY;
	typedef const TCHAR *ARG_KEY;

protected:
	// Association
	struct CAssoc
	{
		CAssoc* pNext;
		unsigned int nHashValue;  // needed for efficient iteration
		KEY key;
		VALUE value;
	};

public:
  // Construction
	CStringMap(int nBlockSize = 10);

// Attributes
	// number of elements
	int GetCount() const;
	bool IsEmpty() const;

	// Lookup
	bool Lookup(ARG_KEY key, VALUE& rValue) const;
	bool LookupKey(ARG_KEY key, ARG_KEY& rKey) const;

// Operations
	// Lookup and add if not there
	VALUE& operator[](ARG_KEY key);

	// add a new (key, value) pair
	void SetAt(ARG_KEY key, ARG_VALUE newValue);

	// removing existing (key, ?) pair
	bool RemoveKey(ARG_KEY key);
	void RemoveAll();

	// iterating all (key, value) pairs
	POSITION GetStartPosition() const;
	void GetNextAssoc(POSITION& rNextPosition, KEY& rKey, VALUE& rValue) const;

	// advanced features for derived classes
	unsigned int GetHashTableSize() const;
	void InitHashTable(unsigned int hashSize, bool bAllocNow = true);

// Overridables: special non-virtual (see map implementation for details)
	// Routine used to user-provided hash keys
	unsigned int HashKey(ARG_KEY key) const;

// Implementation
protected:
	CAssoc** m_pHashTable;
	unsigned int m_nHashTableSize;
	int m_nCount;
	CAssoc* m_pFreeList;
	struct CMemBlock* m_pBlocks;
	int m_nBlockSize;

	CAssoc* NewAssoc();
	void FreeAssoc(CAssoc*);
	CAssoc* GetAssocAt(ARG_KEY, unsigned int&) const;

public:
	~CStringMap();
};

template<class VALUE, class ARG_VALUE>
inline int CStringMap<VALUE, ARG_VALUE>::GetCount() const
{ return m_nCount; }

template<class VALUE, class ARG_VALUE>
inline bool CStringMap<VALUE, ARG_VALUE>::IsEmpty() const
{ return m_nCount == 0; }

template<class VALUE, class ARG_VALUE>
inline void CStringMap<VALUE, ARG_VALUE>::SetAt(ARG_KEY key, ARG_VALUE newValue)
{ (*this)[key] = newValue; }

template<class VALUE, class ARG_VALUE>
inline POSITION CStringMap<VALUE, ARG_VALUE>::GetStartPosition() const
{ return (m_nCount == 0) ? NULL : BEFORE_START_POSITION; }

template<class VALUE, class ARG_VALUE>
inline unsigned int CStringMap<VALUE, ARG_VALUE>::GetHashTableSize() const
{ return m_nHashTableSize; }

template<class VALUE, class ARG_VALUE>
CStringMap<VALUE, ARG_VALUE>::CStringMap(int nBlockSize)
{
	ASSERT(nBlockSize > 0);

	m_pHashTable = NULL;
	m_nHashTableSize = 17;  // default size
	m_nCount = 0;
	m_pFreeList = NULL;
	m_pBlocks = NULL;
	m_nBlockSize = nBlockSize;
}

template<class VALUE, class ARG_VALUE>
inline unsigned int CStringMap<VALUE, ARG_VALUE>::HashKey(ARG_KEY key) const
{
	unsigned int nHash = 0;
	while (*key)
		nHash = (nHash<<5) + nHash + *key++;
	return nHash;
}

template<class VALUE, class ARG_VALUE>
void CStringMap<VALUE, ARG_VALUE>::InitHashTable(
	unsigned int nHashSize, bool bAllocNow)
//
// Used to force allocation of a hash table or to override the default
//   hash table size of (which is fairly small)
{
	ASSERT(m_nCount == 0);
	ASSERT(nHashSize > 0);

	if (m_pHashTable != NULL)
	{
		// free hash table
		delete[] m_pHashTable;
		m_pHashTable = NULL;
	}

	if (bAllocNow)
	{
		m_pHashTable = new CAssoc* [nHashSize];
		memset(m_pHashTable, 0, sizeof(CAssoc*) * nHashSize);
	}
	m_nHashTableSize = nHashSize;
}

template<class VALUE, class ARG_VALUE>
void CStringMap<VALUE, ARG_VALUE>::RemoveAll()
{
	if (m_pHashTable != NULL)
	{
		// destroy elements
		for (unsigned int nHash = 0; nHash < m_nHashTableSize; nHash++)
		{
			CAssoc* pAssoc;
			for (pAssoc = m_pHashTable[nHash]; pAssoc != NULL;
			  pAssoc = pAssoc->pNext)
			{
				DestructElement(&pAssoc->key);  // free up string data
				DestructElement(&pAssoc->value);
			}
		}

		// free hash table
		delete [] m_pHashTable;
		m_pHashTable = NULL;
	}

	m_nCount = 0;
	m_pFreeList = NULL;
	m_pBlocks->FreeDataChain();
	m_pBlocks = NULL;
}

template<class VALUE, class ARG_VALUE>
CStringMap<VALUE, ARG_VALUE>::~CStringMap()
{
	RemoveAll();
	ASSERT(m_nCount == 0);
}

/////////////////////////////////////////////////////////////////////////////
// Assoc helpers
// same as CList implementation except we store CAssoc's not CNode's
//    and CAssoc's are singly linked all the time

template<class VALUE, class ARG_VALUE>
CStringMap<VALUE, ARG_VALUE>::CAssoc*
CStringMap<VALUE, ARG_VALUE>::NewAssoc()
{
	if (m_pFreeList == NULL)
	{
		// add another block
		CMemBlock* newBlock = CMemBlock::Create(m_pBlocks, m_nBlockSize,
							sizeof(CStringMap<VALUE, ARG_VALUE>::CAssoc));
		// chain them into free list
		CStringMap<VALUE, ARG_VALUE>::CAssoc* pAssoc =
				(CStringMap<VALUE, ARG_VALUE>::CAssoc*) newBlock->data();
		// free in reverse order to make it easier to debug
		pAssoc += m_nBlockSize - 1;
		for (int i = m_nBlockSize-1; i >= 0; i--, pAssoc--)
		{
			pAssoc->pNext = m_pFreeList;
			m_pFreeList = pAssoc;
		}
	}
	ASSERT(m_pFreeList != NULL);  // we must have something

	CStringMap<VALUE, ARG_VALUE>::CAssoc* pAssoc = m_pFreeList;
	m_pFreeList = m_pFreeList->pNext;
	m_nCount++;
	ASSERT(m_nCount > 0);  // make sure we don't overflow
	
	ConstructElement(&pAssoc->key);
	ConstructElement(&pAssoc->value);
	return pAssoc;
}

template<class VALUE, class ARG_VALUE>
void CStringMap<VALUE, ARG_VALUE>::FreeAssoc(CStringMap<VALUE, ARG_VALUE>::CAssoc* pAssoc)
{
	DestructElement(&pAssoc->key);  // free up string data
	DestructElement(&pAssoc->value);

	pAssoc->pNext = m_pFreeList;
	m_pFreeList = pAssoc;
	m_nCount--;
	ASSERT(m_nCount >= 0);  // make sure we don't underflow

	// if no more elements, cleanup completely
	if (m_nCount == 0)
		RemoveAll();
}

template<class VALUE, class ARG_VALUE>
CStringMap<VALUE, ARG_VALUE>::CAssoc*
CStringMap<VALUE, ARG_VALUE>::GetAssocAt(ARG_KEY key, unsigned int& nHash) const
// find association (or return NULL)
{
	nHash = HashKey(key) % m_nHashTableSize;

	if (m_pHashTable == NULL)
		return NULL;

	// see if it exists
	CAssoc* pAssoc;
	for (pAssoc = m_pHashTable[nHash]; pAssoc != NULL; pAssoc = pAssoc->pNext)
	{
		if (pAssoc->key == key)
			return pAssoc;
	}
	return NULL;
}

/////////////////////////////////////////////////////////////////////////////

template<class VALUE, class ARG_VALUE>
bool CStringMap<VALUE, ARG_VALUE>::Lookup(ARG_KEY key, VALUE& rValue) const
{
	unsigned int nHash;
	CAssoc* pAssoc = GetAssocAt(key, nHash);
	if (pAssoc == NULL)
		return false;  // not in map

	rValue = pAssoc->value;
	return true;
}

template<class VALUE, class ARG_VALUE>
bool CStringMap<VALUE, ARG_VALUE>::LookupKey(ARG_KEY key, ARG_KEY& rKey) const
{
	unsigned int nHash;
	CAssoc* pAssoc = GetAssocAt(key, nHash);
	if (pAssoc == NULL)
		return false;  // not in map

	rKey= pAssoc->key;
	return true;
}

template<class VALUE, class ARG_VALUE>
VALUE& CStringMap<VALUE, ARG_VALUE>::operator[](ARG_KEY key)
{
	unsigned int nHash;
	CAssoc* pAssoc;
	if ((pAssoc = GetAssocAt(key, nHash)) == NULL)
	{
		if (m_pHashTable == NULL)
			InitHashTable(m_nHashTableSize);

		// it doesn't exist, add a new Association
		pAssoc = NewAssoc();
		pAssoc->nHashValue = nHash;
		pAssoc->key = key;
		// 'pAssoc->value' is a constructed object, nothing more

		// put into hash table
		pAssoc->pNext = m_pHashTable[nHash];
		m_pHashTable[nHash] = pAssoc;
	}
	return pAssoc->value;  // return new reference
}

template<class VALUE, class ARG_VALUE>
bool CStringMap<VALUE, ARG_VALUE>::RemoveKey(ARG_KEY key)
// remove key - return true if removed
{
	if (m_pHashTable == NULL)
		return false;  // nothing in the table

	CAssoc** ppAssocPrev;
	ppAssocPrev = &m_pHashTable[HashKey(key) % m_nHashTableSize];

	CAssoc* pAssoc;
	for (pAssoc = *ppAssocPrev; pAssoc != NULL; pAssoc = pAssoc->pNext)
	{
		if (pAssoc->key == key)
		{
			// remove it
			*ppAssocPrev = pAssoc->pNext;  // remove from list
			FreeAssoc(pAssoc);
			return true;
		}
		ppAssocPrev = &pAssoc->pNext;
	}
	return false;  // not found
}

/////////////////////////////////////////////////////////////////////////////
// Iterating

template<class VALUE, class ARG_VALUE>
void CStringMap<VALUE, ARG_VALUE>::GetNextAssoc(POSITION& rNextPosition,
	KEY& rKey, VALUE& rValue) const
{
	ASSERT(m_pHashTable != NULL);  // never call on empty map

	CAssoc* pAssocRet = (CAssoc*)rNextPosition;
	ASSERT(pAssocRet != NULL);

	if (pAssocRet == (CAssoc*) BEFORE_START_POSITION)
	{
		// find the first association
		for (unsigned int nBucket = 0; nBucket < m_nHashTableSize; nBucket++)
			if ((pAssocRet = m_pHashTable[nBucket]) != NULL)
				break;
		ASSERT(pAssocRet != NULL);  // must find something
	}

	// find next association
	//ASSERT(AfxIsValidAddress(pAssocRet, sizeof(CAssoc)));
	CAssoc* pAssocNext;
	if ((pAssocNext = pAssocRet->pNext) == NULL)
	{
		// go to next bucket
		for (unsigned int nBucket = pAssocRet->nHashValue + 1;
		  nBucket < m_nHashTableSize; nBucket++)
			if ((pAssocNext = m_pHashTable[nBucket]) != NULL)
				break;
	}

	rNextPosition = (POSITION) pAssocNext;

	// fill in return data
	rKey = pAssocRet->key;
	rValue = pAssocRet->value;
}

/***************************************************************************/

#ifndef MFC_CLASSES
typedef CStringMap<CString, const TCHAR *> CMapStringToString;
#endif

#endif // TYPES_MAP_H
