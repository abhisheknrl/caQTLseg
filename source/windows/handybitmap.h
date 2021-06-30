/***********************************************************************************
 *
 *   handybitmap.h -- A Handy Bitmap Class
 *
 *   Erik Persson, private, 1999
 */

#ifndef WINDOWS_HANDYBITMAP_H
#define WINDOWS_HANDYBITMAP_H

class CHandyBitmap : public CBitmap
{
	int m_w, m_h;
	CDC m_DC;

public:
	CHandyBitmap(int w, int h, CDC *pCompatibleDC= NULL);
	~CHandyBitmap();

	CDC *GetDC();
	int ReleaseDC(CDC *pDC);

	void Detach();

	void Draw(CDC *pDC, int x0, int y0, DWORD dwRop= SRCCOPY );
	void Get(CDC *pDC, int x0, int y0);

	int GetWidth() { return m_w; }
	int GetHeight() { return m_h; }
};

#endif // WINDOWS_HANDYBITMAP_H