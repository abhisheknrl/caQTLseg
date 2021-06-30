/*******************************************************************************
 *
 *	 TaggedTree -- Implements a tree with directory-like addressing
 *
 *   Björn Nilsson, 2002-2008
 */

#ifndef TAGGEDTREE_H
#define TAGGEDTREE_H

// #include <afxcoll.h>
#include "types/vector.h"
#include "types/ref.h"
#include "types/claimable.h"
#include "types/taggedlist.h"

class CTaggedTree : public CClaimable
{
public:
	CString	m_sKey; // Must be non-empty!
	CTaggedList m_vTagList;
protected:
	CVector< CRef<CTaggedTree> > m_vrSubTree;

	void SplitKey(CString sFullKey, CString *sSubTree, CString *sProp) const;
	bool SetPropertyRecursive(const char *aFullKey, CString sValue);

public:
	CTaggedTree();
	CTaggedTree(CTaggedTree &tt);
	CRef<CTaggedTree> CreateCopy() const;

	CRef<CTaggedTree> GetSubTree(const char *aName) const;
	CRef<CTaggedTree> GetSubTree(int nIndex) const
	{	return m_vrSubTree[nIndex]; }
	int GetSubTreeCount() const
	{	return m_vrSubTree.GetSize(); }

	void AddSubTree(CRef<CTaggedTree> rNew)
	{ m_vrSubTree.Add(rNew); }

	bool HasProperty(const char *aFullKey) const;
	CString GetProperty(const char *aFullKey) const;
	bool SetProperty(const char *aFullKey, CString sValue);
	bool DeleteSubTree(const char *aName);
	bool DeleteSubTree(int nIndex);
	bool SetSubTreeRef(const char *aName, CRef<CTaggedTree> rTree, bool bAutoCreate);
	bool RenameSubTree(const char *aName, const char *aNewName);
	bool Merge(const CTaggedList &rSrc, CString sPath);
protected:
	~CTaggedTree(); // Protected. Use release!
};

#endif // TAGGEDTREE_H
