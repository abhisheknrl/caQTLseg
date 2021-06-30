/***************************************************************************
 *
 *   map.cpp -- Japanized MFC Style Subset Map Classes
 *
 *   Erik Persson, 2000
 */

#include <types/map.h>

#ifndef MFC_CLASSES
typedef char BYTE;
#endif

/***************************** CMemBlock ***********************************/

CMemBlock* CMemBlock::Create(CMemBlock*& pHead, UINT nMax, UINT cbElement)
{
	ASSERT(nMax > 0 && cbElement > 0);
	CMemBlock* p = (CMemBlock*) new BYTE[sizeof(CMemBlock) + nMax * cbElement];
			// may throw exception
	p->pNext = pHead;
	pHead = p;  // change head (adds in reverse order for simplicity)
	return p;
}

void CMemBlock::FreeDataChain()     // free this one and links
{
	CMemBlock* p = this;
	while (p != NULL)
	{
		BYTE* bytes = (BYTE*) p;
		CMemBlock* pNext = p->pNext;
		delete[] bytes;
		p = pNext;
	}
}
