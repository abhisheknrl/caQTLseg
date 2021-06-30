/***********************************************************************************
 *
 *   handybitmap.cpp -- A Bitmap Class With Draw capability
 *
 *   Erik Persson, private, 1999
 */

#include <afxwin.h>
#include "handybitmap.h"

CHandyBitmap::CHandyBitmap(int w, int h, CDC *pCompatibleDC /* = NULL */)
{
	// There must always exist a DC upon drawing, so we might just as well create it now

	m_w= w;
	m_h= h;

	if (CreateCompatibleBitmap(pCompatibleDC,w,h))
	{
		if (m_DC.CreateCompatibleDC(pCompatibleDC))
		{	m_DC.SelectObject(this);
			return;
		}
	}
	AfxThrowMemoryException();
}

CHandyBitmap::~CHandyBitmap()
{
	m_DC.DeleteDC();
}

void CHandyBitmap::Detach()
{
	m_DC.SelectObject((HBITMAP) NULL);
	Detach();
};

CDC *CHandyBitmap::GetDC()
{
	return &m_DC;
}

int CHandyBitmap::ReleaseDC(CDC *pDC)
{
	return TRUE; // success
}

void CHandyBitmap::Draw(CDC *pDC, int x0, int y0, DWORD dwRop /* = SRCCOPY */)
{
	(void) pDC->BitBlt(x0,y0,m_w,m_h, &m_DC, 0,0, dwRop );
}

void CHandyBitmap::Get(CDC *pDC, int x0, int y0)
{
	(void) m_DC.BitBlt(0,0,m_w,m_h, pDC, x0,y0, SRCCOPY );
}


