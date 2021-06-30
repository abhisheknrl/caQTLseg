/***********************************************************************************
 *
 *   Clipboard.h -- Clip Board Handling
 *
 *   Erik Persson, private, 1999
 */

#ifndef CLIPBOARD_H
#define CLIPBOARD_H

class CClipboard
{
private:
	HWND m_hMainWnd;

	UINT m_uiFileNameFormat;
	bool HasFormat(UINT uiClipFormat);

public:
	CClipboard(HWND hMainWnd);
	~CClipboard();

	bool HasFileName();
	CString GetFileName();
	void SetFileName(CString s);

	// CF_TEXT functions
	bool SetText(CString s);
	bool GetText(CString &s);
};

#endif CLIPBOARD_H
