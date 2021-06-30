/*******************************************************************************
 *
 *   registry.h -- Load/save user settings
 *
 *   Björn Nilsson, 2002-2008
 */

#ifndef REGISTRY_H
#define REGISTRY_H

#include <types/taggedtree.h>
#include <types/ref.h>

CString GetRegString(CString sKey);
bool SetRegString(CString sKey, CString sValue);
bool HasRegString(CString sKey);
CRef<CTaggedTree> GetRegTree(CString sKey);

// Init etc
bool LoadRegistry(CString sApp, CString sUser);
bool SaveRegistry();

#endif
