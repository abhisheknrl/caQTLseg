/*******************************************************************************
 *
 *	 TaggedTree -- Tagged Tree Data Type
 *
 *   Björn Nilsson, 2002
 */

#include <string.h>
#include "taggedtree.h"

CTaggedTree::CTaggedTree()
{
	m_sKey = "";
}

// Copy constructor.
CTaggedTree::CTaggedTree(CTaggedTree &tt)
{
	m_sKey= tt.m_sKey;
	m_vTagList= tt.m_vTagList;
	
	// Note that we are not using ref-counted copying of
	// the vector since this will not increase the ref-count
	// of the references objects.
	m_vrSubTree.SetSize(tt.m_vrSubTree.GetSize());
	for (int i=m_vrSubTree.GetSize();--i>=0;)
		m_vrSubTree.SetAt(i, tt.m_vrSubTree[i]); 
}

// Copy constructor from reference
CRef<CTaggedTree> CTaggedTree::CreateCopy() const
{
	try
	{
		CRef<CTaggedTree> rCopy;
		CTaggedTree &tree= rCopy.CreateNew();
		tree.m_sKey= m_sKey;
		tree.m_vTagList= m_vTagList;

		// Note that we are not using ref-counted copying of
		// the vector since this will not increase the ref-count
		// of the references objects.
		tree.m_vrSubTree.SetSize(m_vrSubTree.GetSize());
		for (int i= m_vrSubTree.GetSize(); --i>=0 ;)
			tree.m_vrSubTree.SetAt(i, m_vrSubTree[i]); 

		return &tree;
	}
	catch(...)
	{
	}

	return NULL;
}

CTaggedTree::~CTaggedTree()
{
}

// Code inherited from SetSubTreeRef.
//
// NOTE: aName is a full path, aNewName is the new name 
//		 for the last branch, e.g RenameSubTree("a/b/c/d", "e") 
//		 will yield "a/b/c/e".
// NOTE: Renaming the *this tree is done by setting 
//		 this->m_sKey (cannot be done with RenameSubTree).
// NOTE: RenameSubTree provides a simple means of using wild-cards:
//		 For ex: RenameSubTree("a/*/b/d", "e") will rename all
//		 subtrees of "a" that in turn have a subtree "b/d".
//		 This is very useful when renaming components.
bool CTaggedTree::RenameSubTree(const char *aName, const char *aNewName)
{
	if (aNewName[0]==0)
	{
		// Must have a proper name.
		ASSERT(false);
		return false;
	}

	if (!aName || aName[0]==0)
		// Null key or illegal key format.
		return false;

	int nLen = 0;
	while ((aName[nLen] != 0) && ( aName[nLen] != '/'))
		nLen++;

	CString sKey= CString(aName, nLen);
	bool bOk= true;
	for (int i=0;bOk && (i<m_vrSubTree.GetSize()); i++)
	{
		// NOTE: Very primitive wild-card matching.
		if ((sKey[0]=='*') || (m_vrSubTree[i]->m_sKey == sKey))
		{
			CRef<CTaggedTree> rSubTree= m_vrSubTree[i];
			m_vrSubTree.SetAt(i, NULL); // Temporarily release the old reference.
			CTaggedTree &tree= rSubTree.CreateExclusive();
			m_vrSubTree.SetAt(i, rSubTree);
			if (aName[nLen] == '/')
				bOk &= tree.RenameSubTree(&aName[nLen+1], aNewName);
			else
				tree.m_sKey= aNewName;
		}
	}

	return true;
}

// Code inherited from GetSubTree and SetProperty
// Sets the subtree reference at a specified depth.
bool CTaggedTree::SetSubTreeRef(const char *aName, CRef<CTaggedTree> rSet, bool bAutoCreate) 
{
	if (aName && aName[0])
	{
		int nLen = 0;
		while ((aName[nLen] != 0) && ( aName[nLen] != '/'))
			nLen++;

		CString sKey= CString(aName, nLen);
		for (int i=0;i<m_vrSubTree.GetSize();i++)
		{
			if (m_vrSubTree[i]->m_sKey == sKey)
			{
				if (aName[nLen] == '/')
				{
					// Recurse
					CRef<CTaggedTree> rSubTree= m_vrSubTree[i];
					m_vrSubTree.SetAt(i, NULL); // Temporarily release the old reference.
					CTaggedTree &tree= rSubTree.CreateExclusive();
					m_vrSubTree.SetAt(i, rSubTree);
					return tree.SetSubTreeRef(&aName[nLen+1], rSet, bAutoCreate);
				}
				else
				{
					// Set reference
					m_vrSubTree.SetAt(i, rSet);
					if (!rSet)
						m_vrSubTree.Delete(i);

					return true;
				}
			}
		}

		// Subtree not found.
		if (aName[nLen] == '/') 
		{
			if (bAutoCreate)
			{
				CRef<CTaggedTree> rSubTree;
				CTaggedTree &tree= rSubTree.CreateNew();
				tree.m_sKey= sKey;
				AddSubTree(rSubTree);

				return tree.SetSubTreeRef(&aName[nLen+1], rSet, bAutoCreate);
			}
		}
		else
		{
			if (rSet)
				AddSubTree(rSet);
			return true;
		}
	}

	// Null key or illegal key format
	return false;
}

bool CTaggedTree::DeleteSubTree(const char *aName) 
{
	return SetSubTreeRef(aName, NULL, false);
}

bool CTaggedTree::DeleteSubTree(int nIndex)
{
	if (nIndex>=0 && nIndex<m_vrSubTree.GetSize())
	{
		m_vrSubTree.SetAt(nIndex, NULL);
		m_vrSubTree.Delete(nIndex);
		return true;
	}

	ASSERT(false);
	return false;
}

CRef<CTaggedTree> CTaggedTree::GetSubTree(const char *aName) const
{
	if (aName && aName[0])
	{
		int nLen = 0;
		while ((aName[nLen] != 0) && ( aName[nLen] != '/'))
			nLen++;

		CString sKey= CString(aName, nLen);
		for (int i=0;i<m_vrSubTree.GetSize();i++)
			if (m_vrSubTree[i]->m_sKey == sKey)
			{
				if (aName[nLen] == '/')
					return m_vrSubTree[i]->GetSubTree(&aName[nLen+1]);
				else
					return m_vrSubTree[i];
			}

		return NULL;
	}

	CRef<CTaggedTree> rRoot = (CTaggedTree *) this;
	return rRoot;
}

// Splits a full key specification into subtree and property keys.
void CTaggedTree::SplitKey(CString sFullKey, CString *sSubTree, CString *sProp) const
{
	int i= sFullKey.GetLength();
	while (--i>=0 && sFullKey[i] != '/') {}

	if (i>=0)
	{
		*sProp = sFullKey.Right(sFullKey.GetLength()-1-i);
		*sSubTree = sFullKey.Left(i);
	}
	else
	{
		*sProp = sFullKey;
		*sSubTree = "";
	}
}

bool CTaggedTree::HasProperty(const char *aFullKey) const
{
	// Returns true if prop exists
	CString sSubTree;
	CString sProp;
	SplitKey(aFullKey, &sSubTree, &sProp);
	
	CRef<CTaggedTree> rTree= GetSubTree(sSubTree);
	if (rTree)
		return rTree->m_vTagList.FindKey(sProp)>=0;
	return false;
}

CString CTaggedTree::GetProperty(const char *aFullKey) const
{
	CString sSubTree;
	CString sProp;
	SplitKey(aFullKey, &sSubTree, &sProp);
	
	CRef<CTaggedTree> rTree= GetSubTree(sSubTree);
	if (rTree)
		return rTree->m_vTagList.GetProperty(sProp);
	else
		return "";
}

bool CTaggedTree::SetProperty(const char *aFullKey, CString sValue)
{
	if (HasProperty(aFullKey) && GetProperty(aFullKey)==sValue)
		return true;

	bool bOk= SetPropertyRecursive(aFullKey, sValue);
	ASSERT(GetProperty(aFullKey)==sValue);
	return bOk;
}

// Sets a property. Creates all the required subtrees.
bool CTaggedTree::SetPropertyRecursive(const char *aFullKey, CString sValue)
{
	if (aFullKey && aFullKey[0])
	{
		int nLen = 0;
		while ((aFullKey[nLen] != 0) && ( aFullKey[nLen] != '/'))
			nLen++;

		if ( aFullKey[nLen] == '/' )
		{			
			// We have a subtree key.
			CString sKey= CString(aFullKey, nLen);

			for (int i=0;i<m_vrSubTree.GetSize();i++)
			{
				if (m_vrSubTree[i]->m_sKey == sKey)
				{
					CRef<CTaggedTree> rOld= m_vrSubTree[i];
					// TODO: Någonsin exclusiv pga kopieringen till rOld?
					if (rOld.IsExclusive())
					{
						// If exclusive, do not copy.
						CTaggedTree &subtree= rOld.GetExclusive();
						return subtree.SetPropertyRecursive(&aFullKey[nLen+1], sValue);
					}
					else
					{
						// Subtree not exclusive. Create copy.
						CTaggedTree &subtree= rOld.CreateCopy();
						if (subtree.SetPropertyRecursive(&aFullKey[nLen+1], sValue))
						{
							CRef<CTaggedTree> rCopy= &subtree;
							m_vrSubTree.SetAt(i, rCopy);
							return true;
						}
						else
							return false;
					}
				}
			}

			// Subtree not found. Create it.
			CRef<CTaggedTree> rSubTree;
			try
			{			
				CTaggedTree *pTree= new CTaggedTree();
				
				ASSERT(!sKey.IsEmpty());
				pTree->m_sKey= sKey;
				rSubTree= pTree;
				pTree->Release();
			}
			catch(...)
			{
				return false;
			}

			CTaggedTree &subtree = rSubTree.GetExclusive();
			if (subtree.SetPropertyRecursive(&aFullKey[nLen+1], sValue))
			{
				// TRACE("Adding subtree %s to tree %s\n", rSubTree->m_sKey, m_sKey);
				AddSubTree(rSubTree);
				return true;
			}
			else
				return false;
		}
	}

	// No subtree key.
	return m_vTagList.SetProperty(aFullKey, sValue);
}

bool CTaggedTree::Merge(const CTaggedList &rSrc, CString sPath)
{
	// Copies all keys from rSrc that do not already 
	// exist in the subtree sPath of this tree. 
	// Note that existance is equivalent to non-emptiness here.
	sPath += "/";
	int I=rSrc.GetCount();
	for (int i=0;i<I;i++)
	{
		CString sProp= sPath + rSrc.GetKey(i);
		CString s= GetProperty(sProp);
		if (s.IsEmpty())
		{
			if (!SetProperty(sProp, rSrc.GetValue(i)))
				return false;
		}
	}

	return true;
}
