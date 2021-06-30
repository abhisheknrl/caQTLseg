/*******************************************************************************
 *
 *   registry.cpp -- Hides registry operations
 *
 *   Björn Nilsson, 2002
 */

#include <afxwin.h>
#include <atlbase.h>
#include <types/stringfun.h>
#include <types/taggedtree.h>
#include <format/treeparser.h>
#include <format/treecomposer.h>
#include <system/filesystem.h>
#include <system/environment.h>
#include <system/stdout.h>
#include <system/stderr.h>
#include "registry.h"

static CRef<CTaggedTree> g_rRegistry;
static CString g_sUserFile;
static CString g_UserPath;

CRef<CTaggedTree> GetRegTree(CString sSubTree)
{
	if (!g_rRegistry)
		return CRef<CTaggedTree>(NULL);
	if (sSubTree.IsEmpty())
		return g_rRegistry;
	else
		return g_rRegistry->GetSubTree(sSubTree);
}

// Helper routine
void SplitKey(CString sFullKey, CString &sSubTree, CString &sProp) 
{
	int i= sFullKey.GetLength();
	while (--i>=0 && sFullKey[i] != '\\') {}

	if (i>=0)
	{
		sProp= sFullKey.Right(sFullKey.GetLength()-1-i);
		sSubTree= sFullKey.Left(i);
	}
	else
	{
		sProp= sFullKey;
		sSubTree= "";
	}
}

bool OpenRegistry()
{	
	if (g_rRegistry)
		return true;

	// Allocate tree.
	try
	{
		g_rRegistry.CreateNew();
		CTaggedTree &tree= g_rRegistry.GetExclusive(); 
		tree.m_sKey= "registry";
		return true;
	}
	catch(...)
	{
		return false;
	}
}

/***********************************************************************************/

// Returns false if app configuration could not be 
// initialized. Succesful loading of the user
// file is not necessary.
bool LoadRegistry(CString sApp, CString sUser)
{
	if (!OpenRegistry())
		return false;

	bool bOk= false;
	CRef<CTaggedTree> rApp= ::LoadTreeFile(sApp);
	if (!rApp || rApp->m_sKey != "app")
		return false;

	CRef<CTaggedTree> rUser= ::LoadTreeFile(sUser);
	g_sUserFile= sUser;

	CTaggedTree &tree= g_rRegistry.CreateExclusive();
	if (rApp)
		bOk= tree.SetSubTreeRef("app", rApp, true);
	else
		bOk= false;

	if (rUser)
		tree.SetSubTreeRef("user", rUser, true);

	return bOk;
}

bool SaveRegistry()
{
	bool bOk= false;
	if (g_rRegistry)
	{
		CString sErr;
		CRef<CTaggedTree> rUser= g_rRegistry->GetSubTree("user");
		if (rUser && !g_sUserFile.IsEmpty())
		{
			bOk= ::SaveTreeFile(rUser, g_sUserFile, &sErr);
			if (!bOk)
				eprintf("Could not write %s (%s)\n",
				(const char *) g_sUserFile, (const char *) sErr);
		}
	}

	return bOk;
}

// Get registry string
CString GetRegString(CString sKey)
{
	/*
	if (USE_WINDOWSREGISTRY)
	{
		sKey = GetRegistryPath() + sKey;
		sKey= ToBackSlash(sKey);

		// Use the Windows registry
		sKey= ToBackSlash(sKey);
		CString sSubTree, sField;
		SplitKey(sKey, sSubTree, sField);
		CRegKey k;
		if (k.Open(HKEY_CURRENT_USER, sSubTree)==ERROR_SUCCESS)
		{
			_TCHAR buf[10000];
			unsigned long n= sizeof(buf);
			if (k.QueryValue(&buf[0], sField, &n)==ERROR_SUCCESS)
			{
				buf[n+1]= 0; // This is merely a safety precaution. buf[n] should already be zero.
				return CString(&buf[0]);
			}
		}
		return "";
	}
	else
	*/
	{
		// Use our own tree.
		sKey= ToForwardSlash(sKey);
		return g_rRegistry->GetProperty(sKey);
	}
}

// Check if registry string has a value
bool HasRegString(CString sKey)
{
	// Use our own tree.
	sKey= ToForwardSlash(sKey);
	return g_rRegistry->HasProperty(sKey);
}

// Set registry string
bool SetRegString(CString sKey, CString sValue)
{
	/*
	if (USE_WINDOWSREGISTRY)
	{
		sKey = GetRegistryPath() + sKey;
		sKey= ToBackSlash(sKey);

		// Use the Windows registry
		sKey= ToBackSlash(sKey);
		CString sSubTree, sField;
		SplitKey(sKey, sSubTree, sField);

		CRegKey k;
		if (k.Open(HKEY_CURRENT_USER, "")==ERROR_SUCCESS)
		{
			if (k.SetKeyValue(sSubTree, sValue, sField)==ERROR_SUCCESS)
				return true;
		}

		return false;
	}
	else
	*/
	{
		// Use our own tree.
		sKey= ToForwardSlash(sKey);
		CTaggedTree &tree= g_rRegistry.CreateExclusive();
		return tree.SetProperty(sKey, sValue);
	}
}

