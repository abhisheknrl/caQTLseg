/*******************************************************************************
 *
 *   uxThemeWrapper.cpp -- Wraps uxtheme.dll
 *
 *   Björn Nilsson, 2003
 *
 *	 Based on code by Pierre ARNAUD, OPaC bright ideas, Switzerland :
 *		"This code may be reused for any use (commercial or free software)."
 *
 */

// #define UXTHEMEWRAPPER_IMPLEMENTATION

#include "uxThemeWrapper.h"

/*****************************************************************************/
// The schema file expects these to be defined by the user.

#define TMT_ENUMDEF 8
#define TMT_ENUMVAL TEXT('A')
#define TMT_ENUM	TEXT('B')

#include "tmschema.h"
#define SCHEMA_STRINGS
#include "tmschema.h" 

/*****************************************************************************
 * 
 *  Pierre ARNAUD wrote: 
 *
 *	The code below has been adapted from Microsoft's sample "ThemeExplorer"
 *	which was no longer available online when I implemented this (June 2002).
 *
 *	I therefore took it from a CodeProject sample :
 *
 *	  "Add XP Visual Style Support to OWNERDRAW Controls"
 */


/*****************************************************************************/
// The one and only VisualStyles handler.
// This was declared static in the dll version. /bn

VisualStylesXP xpStyle;
// static VisualStylesXP xpStyle;

/*****************************************************************************/

VisualStylesXP::VisualStylesXP(void)
{
	m_hThemeDll = LoadLibrary("UxTheme.dll");
}

VisualStylesXP::~VisualStylesXP(void)
{
	if (m_hThemeDll!=NULL)
		FreeLibrary(m_hThemeDll);
	m_hThemeDll = NULL;
}

void*
VisualStylesXP::GetProc(LPCSTR szProc, void* pfnFail)
{
	void* pRet = pfnFail;
	if (m_hThemeDll != NULL)
		pRet = GetProcAddress(m_hThemeDll, szProc);
	return pRet;
}

HTHEME
VisualStylesXP::OpenThemeData(HWND hwnd, LPCWSTR pszClassList)
{
	PFNOPENTHEMEDATA pfnOpenThemeData = (PFNOPENTHEMEDATA)GetProc("OpenThemeData", (void*)OpenThemeDataFail);
	return (*pfnOpenThemeData)(hwnd, pszClassList);
}

HRESULT
VisualStylesXP::CloseThemeData(HTHEME hTheme)
{
	PFNCLOSETHEMEDATA pfnCloseThemeData = (PFNCLOSETHEMEDATA)GetProc("CloseThemeData", (void*)CloseThemeDataFail);
	return (*pfnCloseThemeData)(hTheme);
}

HRESULT
VisualStylesXP::DrawThemeBackground(HTHEME hTheme, HDC hdc,
									int iPartId, int iStateId,
									const RECT *pRect, const RECT *pClipRect)
{
	PFNDRAWTHEMEBACKGROUND pfnDrawThemeBackground = 
		(PFNDRAWTHEMEBACKGROUND)GetProc("DrawThemeBackground", (void*)DrawThemeBackgroundFail);
	return (*pfnDrawThemeBackground)(hTheme, hdc, iPartId, iStateId, pRect, pClipRect);
}


HRESULT
VisualStylesXP::DrawThemeText(HTHEME hTheme, HDC hdc,
							  int iPartId, int iStateId,
							  LPCWSTR pszText, int iCharCount, DWORD dwTextFlags,
							  DWORD dwTextFlags2, const RECT *pRect)
{
	PFNDRAWTHEMETEXT pfn = (PFNDRAWTHEMETEXT)GetProc("DrawThemeText", (void*)DrawThemeTextFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, pszText, iCharCount, dwTextFlags, dwTextFlags2, pRect);
}

HRESULT
VisualStylesXP::GetThemeBackgroundContentRect(HTHEME hTheme,  HDC hdc,
											  int iPartId, int iStateId,
											  const RECT *pBoundingRect, RECT *pContentRect)
{
	PFNGETTHEMEBACKGROUNDCONTENTRECT pfn = (PFNGETTHEMEBACKGROUNDCONTENTRECT)GetProc("GetThemeBackgroundContentRect", (void*)GetThemeBackgroundContentRectFail);
	return (*pfn)(hTheme,  hdc, iPartId, iStateId,  pBoundingRect, pContentRect);
}

HRESULT
VisualStylesXP::GetThemeBackgroundExtent(HTHEME hTheme,  HDC hdc,
										 int iPartId, int iStateId,
										 const RECT *pContentRect, RECT *pExtentRect)
{
	PFNGETTHEMEBACKGROUNDEXTENT pfn = (PFNGETTHEMEBACKGROUNDEXTENT)GetProc("GetThemeBackgroundExtent", (void*)GetThemeBackgroundExtentFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, pContentRect, pExtentRect);
}

HRESULT
VisualStylesXP::GetThemePartSize(HTHEME hTheme, HDC hdc,
								 int iPartId, int iStateId,
								 RECT * pRect, enum THEMESIZE eSize, SIZE *psz)
{
	PFNGETTHEMEPARTSIZE pfnGetThemePartSize = 
		(PFNGETTHEMEPARTSIZE)GetProc("GetThemePartSize", (void*)GetThemePartSizeFail);
	return (*pfnGetThemePartSize)(hTheme, hdc, iPartId, iStateId, pRect, eSize, psz);
}

HRESULT
VisualStylesXP::GetThemeTextExtent(HTHEME hTheme, HDC hdc,
								   int iPartId, int iStateId,
								   LPCWSTR pszText, int iCharCount, DWORD dwTextFlags,
								   const RECT *pBoundingRect, RECT *pExtentRect)
{
	PFNGETTHEMETEXTEXTENT pfn = (PFNGETTHEMETEXTEXTENT)GetProc("GetThemeTextExtent", (void*)GetThemeTextExtentFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, pszText, iCharCount, dwTextFlags,  pBoundingRect, pExtentRect);
}

HRESULT
VisualStylesXP::GetThemeTextMetrics(HTHEME hTheme,  HDC hdc,
									int iPartId, int iStateId,
									TEXTMETRIC* ptm)
{
	PFNGETTHEMETEXTMETRICS pfn = (PFNGETTHEMETEXTMETRICS)GetProc("GetThemeTextMetrics", (void*)GetThemeTextMetricsFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId,  ptm);
}

HRESULT
VisualStylesXP::GetThemeBackgroundRegion(HTHEME hTheme,  HDC hdc,
										 int iPartId, int iStateId,
										 const RECT *pRect,  HRGN *pRegion)
{
	PFNGETTHEMEBACKGROUNDREGION pfn = (PFNGETTHEMEBACKGROUNDREGION)GetProc("GetThemeBackgroundRegion", (void*)GetThemeBackgroundRegionFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, pRect, pRegion);
}

HRESULT
VisualStylesXP::HitTestThemeBackground(HTHEME hTheme,  HDC hdc,
									   int iPartId, int iStateId,
									   DWORD dwOptions, const RECT *pRect, HRGN hrgn, POINT ptTest,  WORD *pwHitTestCode)
{
	PFNHITTESTTHEMEBACKGROUND pfn = (PFNHITTESTTHEMEBACKGROUND)GetProc("HitTestThemeBackground", (void*)HitTestThemeBackgroundFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, dwOptions, pRect, hrgn, ptTest, pwHitTestCode);
}

HRESULT
VisualStylesXP::DrawThemeEdge(HTHEME hTheme, HDC hdc,
							  int iPartId, int iStateId,
							  const RECT *pDestRect, UINT uEdge, UINT uFlags, RECT *pContentRect)
{
	PFNDRAWTHEMEEDGE pfn = (PFNDRAWTHEMEEDGE)GetProc("DrawThemeEdge", (void*)DrawThemeEdgeFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, pDestRect, uEdge, uFlags, pContentRect);
}

HRESULT
VisualStylesXP::DrawThemeIcon(HTHEME hTheme, HDC hdc,
							  int iPartId, int iStateId,
							  const RECT *pRect, HIMAGELIST himl, int iImageIndex)
{
	PFNDRAWTHEMEICON pfn = (PFNDRAWTHEMEICON)GetProc("DrawThemeIcon", (void*)DrawThemeIconFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, pRect, himl, iImageIndex);
}

BOOL
VisualStylesXP::IsThemePartDefined(HTHEME hTheme,
								   int iPartId, int iStateId)
{
	PFNISTHEMEPARTDEFINED pfn = (PFNISTHEMEPARTDEFINED)GetProc("IsThemePartDefined", (void*)IsThemePartDefinedFail);
	return (*pfn)(hTheme, iPartId, iStateId);
}

BOOL
VisualStylesXP::IsThemeBackgroundPartiallyTransparent(HTHEME hTheme,
													  int iPartId, int iStateId)
{
	PFNISTHEMEBACKGROUNDPARTIALLYTRANSPARENT pfn = (PFNISTHEMEBACKGROUNDPARTIALLYTRANSPARENT)GetProc("IsThemeBackgroundPartiallyTransparent", (void*)IsThemeBackgroundPartiallyTransparentFail);
	return (*pfn)(hTheme, iPartId, iStateId);
}

HRESULT
VisualStylesXP::GetThemeColor(HTHEME hTheme,
							  int iPartId, int iStateId,
							  int iPropId, COLORREF *pColor)
{
	PFNGETTHEMECOLOR pfn = (PFNGETTHEMECOLOR)GetProc("GetThemeColor", (void*)GetThemeColorFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, pColor);
}

HRESULT
VisualStylesXP::GetThemeMetric(HTHEME hTheme, HDC hdc,
							   int iPartId, int iStateId,
							   int iPropId,  int *piVal)
{
	PFNGETTHEMEMETRIC pfn = (PFNGETTHEMEMETRIC)GetProc("GetThemeMetric", (void*)GetThemeMetricFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, iPropId, piVal);
}

HRESULT
VisualStylesXP::GetThemeString(HTHEME hTheme,
							   int iPartId, int iStateId,
							   int iPropId,  LPWSTR pszBuff, int cchMaxBuffChars)
{
	PFNGETTHEMESTRING pfn = (PFNGETTHEMESTRING)GetProc("GetThemeString", (void*)GetThemeStringFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, pszBuff, cchMaxBuffChars);
}

HRESULT
VisualStylesXP::GetThemeBool(HTHEME hTheme,
							 int iPartId, int iStateId,
							 int iPropId,  BOOL *pfVal)
{
	PFNGETTHEMEBOOL pfn = (PFNGETTHEMEBOOL)GetProc("GetThemeBool", (void*)GetThemeBoolFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, pfVal);
}

HRESULT
VisualStylesXP::GetThemeInt(HTHEME hTheme,
							int iPartId, int iStateId,
							int iPropId,  int *piVal)
{
	PFNGETTHEMEINT pfn = (PFNGETTHEMEINT)GetProc("GetThemeInt", (void*)GetThemeIntFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, piVal);
}

HRESULT
VisualStylesXP::GetThemeEnumValue(HTHEME hTheme,
								  int iPartId, int iStateId,
								  int iPropId,  int *piVal)
{
	PFNGETTHEMEENUMVALUE pfn = (PFNGETTHEMEENUMVALUE)GetProc("GetThemeEnumValue", (void*)GetThemeEnumValueFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, piVal);
}

HRESULT
VisualStylesXP::GetThemePosition(HTHEME hTheme,
								 int iPartId, int iStateId,
								 int iPropId,  POINT *pPoint)
{
	PFNGETTHEMEPOSITION pfn = (PFNGETTHEMEPOSITION)GetProc("GetThemePosition", (void*)GetThemePositionFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, pPoint);
}

HRESULT
VisualStylesXP::GetThemeFont(HTHEME hTheme,  HDC hdc,
							 int iPartId, int iStateId,
							 int iPropId,  LOGFONT *pFont)
{
	PFNGETTHEMEFONT pfn = (PFNGETTHEMEFONT)GetProc("GetThemeFont", (void*)GetThemeFontFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, iPropId, pFont);
}

HRESULT
VisualStylesXP::GetThemeRect(HTHEME hTheme,
							 int iPartId, int iStateId,
							 int iPropId,  RECT *pRect)
{
	PFNGETTHEMERECT pfn = (PFNGETTHEMERECT)GetProc("GetThemeRect", (void*)GetThemeRectFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, pRect);
}

HRESULT
VisualStylesXP::GetThemeMargins(HTHEME hTheme, HDC hdc,
								int iPartId, int iStateId,
								int iPropId,  RECT *prc,  MARGINS *pMargins)
{
	PFNGETTHEMEMARGINS pfn = (PFNGETTHEMEMARGINS)GetProc("GetThemeMargins", (void*)GetThemeMarginsFail);
	return (*pfn)(hTheme, hdc, iPartId, iStateId, iPropId, prc, pMargins);
}

HRESULT
VisualStylesXP::GetThemeIntList(HTHEME hTheme,
								int iPartId, int iStateId,
								int iPropId,  INTLIST *pIntList)
{
	PFNGETTHEMEINTLIST pfn = (PFNGETTHEMEINTLIST)GetProc("GetThemeIntList", (void*)GetThemeIntListFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, pIntList);
}

HRESULT
VisualStylesXP::GetThemePropertyOrigin(HTHEME hTheme,
									   int iPartId, int iStateId,
									   int iPropId,  enum PROPERTYORIGIN *pOrigin)
{
	PFNGETTHEMEPROPERTYORIGIN pfn = (PFNGETTHEMEPROPERTYORIGIN)GetProc("GetThemePropertyOrigin", (void*)GetThemePropertyOriginFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId, pOrigin);
}

HRESULT
VisualStylesXP::SetWindowTheme(HWND hwnd, LPCWSTR pszSubAppName, LPCWSTR pszSubIdList)
{
	PFNSETWINDOWTHEME pfn = (PFNSETWINDOWTHEME)GetProc("SetWindowTheme", (void*)SetWindowThemeFail);
	return (*pfn)(hwnd, pszSubAppName, pszSubIdList);
}

HRESULT
VisualStylesXP::GetThemeFilename(HTHEME hTheme,
								 int iPartId, int iStateId,
								 int iPropId,  LPWSTR pszThemeFileName, int cchMaxBuffChars)
{
	PFNGETTHEMEFILENAME pfn = (PFNGETTHEMEFILENAME)GetProc("GetThemeFilename", (void*)GetThemeFilenameFail);
	return (*pfn)(hTheme, iPartId, iStateId, iPropId,  pszThemeFileName, cchMaxBuffChars);
}

COLORREF
VisualStylesXP::GetThemeSysColor(HTHEME hTheme, int iColorId)
{
	PFNGETTHEMESYSCOLOR pfn = (PFNGETTHEMESYSCOLOR)GetProc("GetThemeSysColor", (void*)GetThemeSysColorFail);
	return (*pfn)(hTheme, iColorId);
}

HBRUSH
VisualStylesXP::GetThemeSysColorBrush(HTHEME hTheme, int iColorId)
{
	PFNGETTHEMESYSCOLORBRUSH pfn = (PFNGETTHEMESYSCOLORBRUSH)GetProc("GetThemeSysColorBrush", (void*)GetThemeSysColorBrushFail);
	return (*pfn)(hTheme, iColorId);
}

BOOL
VisualStylesXP::GetThemeSysBool(HTHEME hTheme, int iBoolId)
{
	PFNGETTHEMESYSBOOL pfn = (PFNGETTHEMESYSBOOL)GetProc("GetThemeSysBool", (void*)GetThemeSysBoolFail);
	return (*pfn)(hTheme, iBoolId);
}

int
VisualStylesXP::GetThemeSysSize(HTHEME hTheme, int iSizeId)
{
	PFNGETTHEMESYSSIZE pfn = (PFNGETTHEMESYSSIZE)GetProc("GetThemeSysSize", (void*)GetThemeSysSizeFail);
	return (*pfn)(hTheme, iSizeId);
}

HRESULT
VisualStylesXP::GetThemeSysFont(HTHEME hTheme, int iFontId,  LOGFONT *plf)
{
	PFNGETTHEMESYSFONT pfn = (PFNGETTHEMESYSFONT)GetProc("GetThemeSysFont", (void*)GetThemeSysFontFail);
	return (*pfn)(hTheme, iFontId, plf);
}

HRESULT
VisualStylesXP::GetThemeSysString(HTHEME hTheme, int iStringId, LPWSTR pszStringBuff, int cchMaxStringChars)
{
	PFNGETTHEMESYSSTRING pfn = (PFNGETTHEMESYSSTRING)GetProc("GetThemeSysString", (void*)GetThemeSysStringFail);
	return (*pfn)(hTheme, iStringId, pszStringBuff, cchMaxStringChars);
}

HRESULT
VisualStylesXP::GetThemeSysInt(HTHEME hTheme, int iIntId, int *piValue)
{
	PFNGETTHEMESYSINT pfn = (PFNGETTHEMESYSINT)GetProc("GetThemeSysInt", (void*)GetThemeSysIntFail);
	return (*pfn)(hTheme, iIntId, piValue);
}

BOOL
VisualStylesXP::IsThemeActive()
{
	PFNISTHEMEACTIVE pfn = (PFNISTHEMEACTIVE)GetProc("IsThemeActive", (void*)IsThemeActiveFail);
	return (*pfn)();
}

BOOL
VisualStylesXP::IsAppThemed()
{
	PFNISAPPTHEMED pfnIsAppThemed = (PFNISAPPTHEMED)GetProc("IsAppThemed", (void*)IsAppThemedFail);
	return (*pfnIsAppThemed)();
}

HTHEME
VisualStylesXP::GetWindowTheme(HWND hwnd)
{
	PFNGETWINDOWTHEME pfn = (PFNGETWINDOWTHEME)GetProc("GetWindowTheme", (void*)GetWindowThemeFail);
	return (*pfn)(hwnd);
}

HRESULT
VisualStylesXP::EnableThemeDialogTexture(HWND hwnd, DWORD dwFlags)
{
	PFNENABLETHEMEDIALOGTEXTURE pfn = (PFNENABLETHEMEDIALOGTEXTURE)GetProc("EnableThemeDialogTexture", (void*)EnableThemeDialogTextureFail);
	return (*pfn)(hwnd, dwFlags);
}

BOOL
VisualStylesXP::IsThemeDialogTextureEnabled(HWND hwnd)
{
	PFNISTHEMEDIALOGTEXTUREENABLED pfn = (PFNISTHEMEDIALOGTEXTUREENABLED)GetProc("IsThemeDialogTextureEnabled", (void*)IsThemeDialogTextureEnabledFail);
	return (*pfn)(hwnd);
}

DWORD
VisualStylesXP::GetThemeAppProperties()
{
	PFNGETTHEMEAPPPROPERTIES pfn = (PFNGETTHEMEAPPPROPERTIES)GetProc("GetThemeAppProperties", (void*)GetThemeAppPropertiesFail);
	return (*pfn)();
}

void
VisualStylesXP::SetThemeAppProperties(DWORD dwFlags)
{
	PFNSETTHEMEAPPPROPERTIES pfn = (PFNSETTHEMEAPPPROPERTIES)GetProc("SetThemeAppProperties", (void*)SetThemeAppPropertiesFail);
	(*pfn)(dwFlags);
}

HRESULT
VisualStylesXP::GetCurrentThemeName(LPWSTR pszThemeFileName, int cchMaxNameChars,
									LPWSTR pszColorBuff, int cchMaxColorChars,
									LPWSTR pszSizeBuff, int cchMaxSizeChars)
{
	PFNGETCURRENTTHEMENAME pfn = (PFNGETCURRENTTHEMENAME)GetProc("GetCurrentThemeName", (void*)GetCurrentThemeNameFail);
	return (*pfn)(pszThemeFileName, cchMaxNameChars, pszColorBuff, cchMaxColorChars, pszSizeBuff, cchMaxSizeChars);
}

HRESULT
VisualStylesXP::GetThemeDocumentationProperty(LPCWSTR pszThemeName,
											  LPCWSTR pszPropertyName,
											  LPWSTR pszValueBuff, int cchMaxValChars)
{
	PFNGETTHEMEDOCUMENTATIONPROPERTY pfn = (PFNGETTHEMEDOCUMENTATIONPROPERTY)GetProc("GetThemeDocumentationProperty", (void*)GetThemeDocumentationPropertyFail);
	return (*pfn)(pszThemeName, pszPropertyName, pszValueBuff, cchMaxValChars);
}


HRESULT
VisualStylesXP::DrawThemeParentBackground(HWND hwnd, HDC hdc,  RECT* prc)
{
	PFNDRAWTHEMEPARENTBACKGROUND pfn = (PFNDRAWTHEMEPARENTBACKGROUND)GetProc("DrawThemeParentBackground", (void*)DrawThemeParentBackgroundFail);
	return (*pfn)(hwnd, hdc, prc);
}

HRESULT
VisualStylesXP::EnableTheming(BOOL fEnable)
{
	PFNENABLETHEMING pfn = (PFNENABLETHEMING)GetProc("EnableTheming", (void*)EnableThemingFail);
	return (*pfn)(fEnable);
}

/*****************************************************************************/

bool Themes_FindVisualStyle(const wchar_t* class_name, 
								   const wchar_t* find_part, 
								   const wchar_t* find_state,
								   int & part_id, 
								   int & state_id)
{
	const TMSCHEMAINFO *pSchemaInfo = GetSchemaInfo();
	const int nSchCnt = pSchemaInfo->iPropCount;
	const wchar_t szParts[] = L"PARTS";
	const wchar_t szStates[] = L"STATES";
	const TMPROPINFO* pPropTable = pSchemaInfo->pPropTable;
	
	wchar_t findClassName[35];

	wcscpy (findClassName, class_name);
	wcscat (findClassName, szParts);
	
	int i = 0;
	
	//	Move past the items at the beginning of the file.
	while ((i < nSchCnt)
		&& (!wcsstr (pSchemaInfo->pPropTable[i].pszName, szParts)))
	{
		i++;
	}

	if (i == nSchCnt)
	{	
		return FALSE;	// No parts were found
	}
	
	for (int j = i; j < nSchCnt; j++)
	{
		if (wcscmp (findClassName, pPropTable[j].pszName) == 0)
		{
			// Get the theme handle if it exists for this item.
			// Open a handle to the theme data.
			
			HTHEME theme = NULL;
			
			if (xpStyle.IsAppThemed ())	{
				theme = xpStyle.OpenThemeData(NULL, class_name);
			}
			
			if (theme)  // If the application does not find a handle, return.
			{
				// Retrieve the children.
				
				j++;

				int nPartNum = 0;  // This equates to the enum value for the part.
				
				while (!(wcsstr(pPropTable[j].pszName, szParts)) &&
					   !(wcsstr(pPropTable[j].pszName, szStates)))
				{
					nPartNum++;  // Part numbers start with number 1.
						
					// Retrieve the states.
					int k = j++;
					
					const wchar_t* current_name = pPropTable[k].pszName;
					
					wchar_t temp[45];
					wcscpy (temp, find_part);
					wcscat (temp, szStates);  // temp now equals partSTATES.
					
					if (find_state)
					{
						while ((k < nSchCnt)
							&& (!(wcsstr(pPropTable[k++].pszName, temp))))
						{
							//	skip...
						}
					}
					else
					{
						k = nSchCnt;
					}
					
					if (k < nSchCnt)  //The item is found.
					{
						// Fill in states.
						int nStateNum = 0;
						while ((k < nSchCnt)
							&& (!(wcsstr(pPropTable[k].pszName, szParts)))
							&& (!(wcsstr(pPropTable[k].pszName, szStates))))
						{
							nStateNum++;
							
							if (wcscmp (find_part, pPropTable[j-1].pszName) == 0)
							{
								if (wcscmp (find_state, pPropTable[k].pszName) == 0)
								{
									part_id  = nPartNum;
									state_id = nStateNum;
									xpStyle.CloseThemeData (theme);
									return true;
								}
							}
							k++;
						}
					}
					else
					{
						if (wcscmp (find_part, current_name) == 0)
						{
							part_id  = nPartNum;
							state_id = 0;
							xpStyle.CloseThemeData (theme);
							return true;
						}
					}
				}
				
				j--;  // The item is incremented in the while statements.
				
				xpStyle.CloseThemeData (theme);
				theme = NULL;
			}
		}
	}		
	return false;
}

/*****************************************************************************/

bool __stdcall Themes_IsAppThemed ()
{
	return xpStyle.IsAppThemed () ? true : false;
}

bool __stdcall Themes_DrawBackground (const wchar_t* name,
									   const wchar_t* part_name, const wchar_t* state_name,
									   void* hdc,
									   int ox, int oy, int dx, int dy,
									   int clip_ox, int clip_oy, int clip_dx, int clip_dy)
{
	HDC hDC = (HDC) hdc;
	bool ok = false;
	
	if (xpStyle.IsAppThemed ())
	{
		HTHEME theme = xpStyle.OpenThemeData (NULL, name);
		
		if (theme != NULL)
		{
			RECT rect;
			RECT clip;
			
			rect.left   = ox;					//	box of full tab page background
			rect.right  = ox+dx;
			rect.top    = oy;
			rect.bottom = oy+dy;
			
			clip = rect;
			
			POINT point;
			
			//	Temporarily shift the coordinate system back into that of the parent (if needed),
			//	to make sure uxTheme draws the background with the same 'phase' as the original
			//	tab page.
			
			GetViewportOrgEx (hDC, &point);
			
			if ((point.x == 0) && (point.y == 0))
			{
				if ((clip_ox != 0) || (clip_oy != 0))
				{
					/*
					rect.left   -= clip_ox;
					rect.right  -= clip_ox;
					rect.top    -= clip_oy;
					rect.bottom -= clip_oy;
					*/
				}
			}
			
			int part  = 0;
			int state = 0;
			
			try
			{
				if (Themes_FindVisualStyle (name, part_name, state_name, part, state))
				{
					ok = (S_OK == xpStyle.DrawThemeBackground (theme, hDC, part, state, &rect, &clip));
				}
			}
			catch (...)
			{
				//	Swallow any exceptions - just in case, so we don't crash the caller if something
				//	gets really wrong.
			}
			
			xpStyle.CloseThemeData (theme);
		}
	}
	
	return ok;
}

bool __stdcall Themes_DrawThemeParentBackground (void* hwnd, void* hdc)
{
	return (S_OK == xpStyle.DrawThemeParentBackground ((HWND)hwnd, (HDC)hdc, NULL));
}

bool __stdcall Themes_DrawThemeParentBackgroundRect (void* hwnd, void* hdc, int ox, int oy, int dx, int dy)
{
	RECT rect;
	
	rect.left   = ox;
	rect.right  = ox+dx;
	rect.top    = oy;
	rect.bottom = oy+dy;
	
	return (S_OK == xpStyle.DrawThemeParentBackground ((HWND)hwnd, (HDC)hdc, &rect));
}

bool __stdcall Themes_GetTextColor (const wchar_t* name,
									 const wchar_t* part_name, const wchar_t* state_name,
									 int* r, int* g, int* b)
{
	bool ok = false;
	
	if (xpStyle.IsAppThemed ())
	{
		HTHEME theme = xpStyle.OpenThemeData (NULL, name);
		
		if (theme != NULL)
		{
			try
			{
				int part;
				int state;
				
				if (Themes_FindVisualStyle (name, part_name, state_name, part, state))
				{
					COLORREF color;
					int prop = TMT_TEXTCOLOR;
					
					if (S_OK == xpStyle.GetThemeColor (theme, part, state, prop, &color))
					{
						*r = GetRValue (color);
						*g = GetGValue (color);
						*b = GetBValue (color);
						ok = true;
					}
				}
			}
			catch (...)
			{
				//	Swallow any exceptions - just in case, so we don't crash the caller if something
				//	gets really wrong.
			}
			
			xpStyle.CloseThemeData (theme);
		}
	}
	
	return ok;
}

/*****************************************************************************/
