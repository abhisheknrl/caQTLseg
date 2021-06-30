/***********************************************************************************
 *
 *   Clipboard.cpp -- Clip Board Handling
 *
 *   Erik Persson, private, 1999
 */

#include <afxwin.h>
#include <shlobj.h>
#include "clipboard.h"

/***********************************************************************************/

CClipboard::CClipboard(HWND hMainWnd)
{
	m_uiFileNameFormat= RegisterClipboardFormat(CFSTR_FILENAME);
	m_hMainWnd= hMainWnd;
}

CClipboard::~CClipboard()
{
}

bool CClipboard::HasFormat(UINT uiClipFormat)
{
	// See if there the specifiedis a filename on a clipboardlist available.

	OpenClipboard(m_hMainWnd);
	UINT uiFmt = 0;
	while (uiFmt= EnumClipboardFormats(uiFmt)) 
	{
		if (uiFmt == uiClipFormat) 
		{
			CloseClipboard();
			return true;
		}
	}
	CloseClipboard();    
	return false;
}

bool CClipboard::HasFileName()
{
    // See if there is a filename on a clipboardlist available.
	return HasFormat(m_uiFileNameFormat);
}

CString CClipboard::GetFileName()
{
	CString s;

	if (OpenClipboard(m_hMainWnd))
	{
		HGLOBAL hGlobText= GetClipboardData(m_uiFileNameFormat); 
		LPSTR aGlobText= (LPSTR) :: GlobalLock(hGlobText);
		
		if (aGlobText)
		{
			int len= 0;
			while (aGlobText[len]) len++; 
			s= CString(aGlobText,len);
		}

		::GlobalUnlock(hGlobText);

		CloseClipboard();	
	}
	return s;
}

void CClipboard::SetFileName(CString s)
{
	if (OpenClipboard(m_hMainWnd))
	{	
		if (EmptyClipboard())
		{
			HGLOBAL hGlobText= ::GlobalAlloc(GMEM_SHARE,s.GetLength()+1);
			LPSTR aGlobText= (LPSTR) :: GlobalLock(hGlobText);
			ASSERT(aGlobText);
			
			strcpy(aGlobText, s);
			::GlobalUnlock(hGlobText);

			if (::SetClipboardData( m_uiFileNameFormat, hGlobText ) == NULL )
			{
				// AfxMessageBox( "Unable to set Clipboard data" );
			}
		}
		CloseClipboard();
	}
}

bool CClipboard::SetText(CString s)
{
    if (!OpenClipboard(m_hMainWnd))
        return false;
    EmptyClipboard();

	// Convert LF -> CR/LF
	CString sSrc;
	char *ps= s.GetBuffer(s.GetLength()+1);
	ps[s.GetLength()]= (char) 0;
	while (*ps)
	{
		if (*ps==10 && *(ps+1) != 13)
			sSrc += (char) 13;
		sSrc += *ps;
		ps++;
	}

	// Allocate global memory handle.
    HGLOBAL hClip= GlobalAlloc(GMEM_MOVEABLE, (sSrc.GetLength()+1) * sizeof(char)); 
    if (!hClip)
    { 
        CloseClipboard(); 
        return false; 
    }

    // Lock handle and copy text to buffer.
    char *pClip= (char *)GlobalLock(hClip);
	char *pSrc= sSrc.GetBuffer(sSrc.GetLength());
    memcpy(pClip, pSrc, sSrc.GetLength() * sizeof(char));
    pClip[sSrc.GetLength()]= (char) 0; // null terminator
    GlobalUnlock(hClip);

    // Place the handle on the clipboard.
    bool bOk= SetClipboardData(CF_TEXT, hClip) != NULL;
	CloseClipboard();
	return bOk;
}

bool CClipboard::GetText(CString &s)
{
	s= "";

    if (!IsClipboardFormatAvailable(CF_TEXT))
		return false;
    if (!OpenClipboard(m_hMainWnd))
        return false;

	bool bOk= false;
    HANDLE hClip= GetClipboardData(CF_TEXT);
    if (hClip)
    {
        char *pClip= (char *)GlobalLock(hClip);
        if (pClip)
        {
			int L= strlen(pClip);
			char *pBuf= s.GetBuffer(L);
			int i=0;
			while (char c= *pClip++)
			{
				if (c!=13)
					pBuf[i++]= c;
			}
			s.ReleaseBuffer(i);

			bOk= true;
            GlobalUnlock(hClip);
        }
    }

    CloseClipboard();
	return bOk;
}
