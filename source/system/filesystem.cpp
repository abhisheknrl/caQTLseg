/***********************************************************************************
 *
 *   filesystem.cpp -- Interface to file system
 *
 *   Erik Persson and Björn Nilsson, 2000-2003
 *   Björn Nilsson, 2003-
 * 
 */

#include "filesystem.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <io.h>
#endif

#include "environment.h"
#include "types/stringfun.h"
#include "types/vector.h"

//POSIX
#ifndef _WIN32
#include <errno.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <memory.h>
// #include <utime.h>
#endif

// WINDOWS
#ifdef _WIN32
#ifdef MFC_CLASSES
#include <afxwin.h>
#else
#include <windows.h> // For various things
#endif
#endif

#define STRBUFSIZE (0x7fff)

CFileSystem g_FileSystem;

/***********************************************************************************/

// I belive this is a free value
// #define ERROR_CUSTOM (0x0000f000) /ep

CString CFileSystem::GetLastError() const
{
#ifdef _WIN32
	int err= ::GetLastError();
	if (err == NO_ERROR)
		return m_sLastCustomError;
	else
	{
		CString s;

		/*if (err==ERROR_FILE_NOT_FOUND)
		{
			// Windows has annoyingly lengthy description for that
			s= "File not found";
		}
		else*/
		{
		void *lpMsgBuf;
		FormatMessage(
			FORMAT_MESSAGE_ALLOCATE_BUFFER |
			FORMAT_MESSAGE_FROM_SYSTEM |
			FORMAT_MESSAGE_IGNORE_INSERTS,
			NULL,
			err,
			MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
			(LPTSTR) &lpMsgBuf,
			0,
			NULL
		);
			s= (const char *) lpMsgBuf;

		LocalFree( lpMsgBuf );

			// Trim EOLs and periods
			while (s.GetLength() && (s[s.GetLength()-1]==13
				 || s[s.GetLength()-1]==10
				 || s[s.GetLength()-1]=='.' ))
				s= CString(s,s.GetLength()-1);
		}
		// The call to ::GetLastError appears to 'consume' the error
		// so we need to cache it
		m_sLastCustomError= s;
		return s;
	}
#else
	// Unix/POSIX errors
	switch (errno)
	{
          case E2BIG: return "Argument list too long";
          case EACCES: return "Permission denied";
          case EADDRINUSE: return "Address in use";
          case EADDRNOTAVAIL: return "Address not available";
          case EAFNOSUPPORT: return "Address family not supported";
          case EAGAIN: return "Resource unavailable";
          case EALREADY: return "Connection already in progress";
          case EBADF: return "Bad file descriptor";
          case EBADMSG: return "Bad message";
          case EBUSY: return "Device or resource busy";
//          case ECANCELED: return "Operation canceled";
          case ECHILD: return "No child processes";
          case ECONNABORTED: return "Connection aborted";
          case ECONNREFUSED: return "Connection refused";
          case ECONNRESET: return "Connection reset";
          case EDEADLK: return "Resource deadlock would occur";
          case EDESTADDRREQ: return "Destination address required";
          case EDOM: return "Mathematics argument out of domain of function";
          case EDQUOT: return "Reserved";
          case EEXIST: return "File exists";
          case EFAULT: return "Bad address";
          case EFBIG: return "File too large";
          case EHOSTUNREACH: return "Host is unreachable";
          case EIDRM: return "Identifier removed";
          case EILSEQ: return "Illegal byte sequence";
          case EINPROGRESS: return "Operation in progress";
          case EINTR: return "Interrupted function";
          case EINVAL: return "Invalid argument";
          case EIO: return "I/O error";
          case EISCONN: return "Socket is connected";
          case EISDIR: return "Is a directory";
          case ELOOP: return "Too many levels of symbolic links";
          case EMFILE: return "Too many open files";
          case EMLINK: return "Too many links";
          case EMSGSIZE: return "Message too large";
          case EMULTIHOP: return "Reserved";
          case ENAMETOOLONG: return "Filename too long";
          case ENETDOWN: return "Network is down";
          case ENETRESET: return "Connection aborted by network";
          case ENETUNREACH: return "Network unreachable";
          case ENFILE: return "Too many files open in system";
          case ENOBUFS: return "No buffer space available";
          case ENODATA: return "No message is available on the STREAM head read queue";
          case ENODEV: return "No such device";
          case ENOENT: return "No such file or directory";
          case ENOEXEC: return "Executable file format error";
          case ENOLCK: return "No locks available";
          case ENOLINK: return "Reserved [ENOLINK]";
          case ENOMEM: return "Not enough space";
          case ENOMSG: return "No message of the desired type";
          case ENOPROTOOPT: return "Protocol not available";
          case ENOSPC: return "No space left on device";
          case ENOSR: return "No STREAM resources";
          case ENOSTR: return "Not a STREAM";
          case ENOSYS: return "Function not supported";
          case ENOTCONN: return "The socket is not connected";
          case ENOTDIR: return "Not a directory";
          case ENOTEMPTY: return "Directory not empty";
          case ENOTSOCK: return "Not a socket";
          case ENOTSUP: return "Operation not supported on socket";
          case ENOTTY: return "Inappropriate I/O control operation";
          case ENXIO: return "No such device or address";
          case EOVERFLOW: return "Value too large to be stored in data type";
          case EPERM: return "Operation not permitted";
          case EPIPE: return "Broken pipe";
          case EPROTO: return "Protocol error";
          case EPROTONOSUPPORT: return "Protocol not supported";
          case EPROTOTYPE: return "Protocol wrong type for socket";
          case ERANGE: return "Result too large";
          case EROFS: return "Read-only file system";
          case ESPIPE: return "Invalid seek";
          case ESRCH: return "No such process";
          case ESTALE: return "Reserved";
          case ETIME: return "Stream ioctl() timeout";
          case ETIMEDOUT: return "Connection timed out";
          case ETXTBSY: return "Text file busy";
//          case EWOULDBLOCK: return "Operation would block";
          case EXDEV: return "Cross-device link";
	}
	return "";
#endif
}

void CFileSystem::SetLastError(const char *ach)
{
#ifdef _WIN32
	m_sLastCustomError= ach;
	::SetLastError(NO_ERROR);
#else
	// TODO: UNIX: Implement
#endif
}

/***********************************************************************************/

CString CFileSystem::GetFullName(const char *ach) const
{
	char aBuf[STRBUFSIZE];
#ifdef _WIN32
	char *pFilePart;

	int n= GetFullPathName(ach, STRBUFSIZE, aBuf, &pFilePart); //Win32 API

	// Fit inside buffer?
	if (n>=0 && n<STRBUFSIZE)
		return CString(aBuf,n);
	return CString("");
#else
	if (realpath(ach, &aBuf[0])!=NULL)
		return CString(aBuf);
	else
	{
		// printf("realpath failed\n");
		return CString(ach); // return original as a best guess...
		// return "";
	}
#endif
}

bool CFileSystem::IsEqualFilename(const char *aFilename0, const char *aFilename1)
{
#ifdef _WIN32
	return strcmp_nocase(g_FileSystem.GetFullName(aFilename0),g_FileSystem.GetFullName(aFilename1))==0;
#else
	return strcmp(g_FileSystem.GetFullName(aFilename0),g_FileSystem.GetFullName(aFilename1))==0;
#endif
}

/***********************************************************************************/

CString PathPrefixFromFullName(const char *aFullName)
{
	size_t n= strlen(aFullName);

	// Get past trailing FILENAME_SLASH in filename itself
	if (n>1 && aFullName[n-1]==FILENAME_SLASH && aFullName[n-2]!=':')
		n--;

	// Skip chars until next FILENAME_SLASH
	while (n>0 && aFullName[n-1]!=FILENAME_SLASH)
		n--;

	return CString(aFullName,n);
}

CString StrippedFileNameFromFullName(const char *aFullName)
{
	// Return the filename itself
	// Or enpty string ("") if s is root ("e:\")

	// Delete trailing backFILENAME_SLASH (except if root)
	size_t n= strlen(aFullName);
	if (n>1 && aFullName[n-1]==FILENAME_SLASH && aFullName[n-2]!=':')
		n--;

	// Get part till next backFILENAME_SLASH
	size_t i= n;
	while (i>1 && !(aFullName[i-1]==FILENAME_SLASH))
		i--;

	return CString(aFullName+i,n-i);
}

const char *ScanSharePrefix(const char *ach)
{
#ifdef _WIN32
	// Double backslash --> Windows share
	if (ach[0]=='\\' && ach[1]=='\\')
		return ach+2;
#endif
	// TODO: UNIX: Do share names ever occur? For ex, when trying to reach Windows over SMB?
	return NULL;
}

// Changed into returning the prefix without glue operator "\" /EP
// When the input is a simple file name, returns non-advanced ach with pPrefix=""
const char *ScanMinimalPathPrefix(const char *ach, CString *pPrefix)
{
	const char *ach_share= ScanSharePrefix(ach);
	if (ach_share!=NULL)
	{
		*pPrefix= "\\"; // Non-empty return value needed for SubtractFileNames
		return ach_share;
	}

	if (ach[0]==FILENAME_SLASH)
		ach++;
	for (int i=0; ach[i]; i++)
		if (ach[i]==FILENAME_SLASH)
		{
			*pPrefix= CString(ach,i); // without slash
			return ach+i+1; // include slash
		}

	// Simple file name
	*pPrefix= "";
	return ach;
}

/**********************************************************************************/

CString CFileSystem::SubtractFullNames(const char *aFullFileName, const char *aFullDirName) const
{
	// Subtraction of normalized file names (as returned by GetFullName())
	const char *aF= aFullFileName;
	const char *aD= aFullDirName;
	const char *aF1, *aD1;
	CString fpre, dpre;
	int similar=0;

	// Advance aF and aD as long as we find equal path prefixes
	aF1= ScanMinimalPathPrefix(aF,&fpre);
	aD1= ScanMinimalPathPrefix(aD,&dpre);
	while (fpre != "" && fpre==dpre)
	{
		// Similar prefixes. Proceed
		similar++;
		aF= aF1;
		aD= aD1;
		aF1= ScanMinimalPathPrefix(aF,&fpre);
		aD1= ScanMinimalPathPrefix(aD,&dpre);
	}

	if (similar==0)
		// Dissimilar drives. Return full file name
		return aFullFileName;

	// ScanMinimalPathPrefix does not accept the last part of a path
	// So we have to check it separately
	if (!strcmp(aF,aD))
	{
		// <filename> = <dirname> (Names cancel fully)
		// Return "." which is the zero of filename algebra
		return ".";
	}
	else if (fpre == CString(aD))
	{
		// <filename> = <dirname> \ <f1> (Dirname is ancestor of Filename)
		// Simple return what's left after <dirname> has been removed
		return aF1;
	}
	else
	{
		// There is a residual in <dirname> that does not cancel out
		CString sResult;

		if (dpre == CString(aF))
		{
			// <dirname> = <filename> \ <d1> (Filename is ancestor of Dirname)
			// Grab away one prefix-part (more are added in while-loop below)
			sResult= ".."; // at least one level
			aD1= ScanMinimalPathPrefix(aD1,&dpre); // scan away one prefix-part to match the ".."
		}
		else
		{
			// <filename> = <common> \ <f> (none is ancestor of the other)
			// <dirname>  = <common> \ <d>
			// Each has an unique remainder after canceling out common part
			sResult= ".." + CString(FILENAME_SLASH) + CString(aF);
		}

		// Prepend extra parenting levels
		while (dpre != "")
		{
			sResult= ".." + CString(FILENAME_SLASH) + sResult;
			aD1= ScanMinimalPathPrefix(aD1,&dpre);
		}

		return sResult;
	}
}

CString CFileSystem::AddFileNames(const char *aFullName, const char *aRelName) const
{
	if (IsAbsPath(aRelName)) // Absolute path
		return aRelName;
	if (!aFullName || aFullName[0]==0) // Nothing to prepend
		return aRelName;

	// Decompose relative path
	CString sResult= aFullName;
	const char *aR= aRelName;
	CString rpre;
	aR= ScanMinimalPathPrefix(aR, &rpre);
	while (rpre != "")
	{
		// rpre will not contain a trailing FILENAME_SLASH
		if (rpre=="..")
			sResult= g_FileSystem.GetParentName(sResult);
		else if (rpre != ".")
		{
			size_t slen= sResult.GetLength();
			if (slen<1 || sResult[slen-1]!=FILENAME_SLASH)
				// Glue needed (Some paths - "C:\" etc already have glue)
				sResult += CString(FILENAME_SLASH);
			sResult += rpre;
		}
		aR= ScanMinimalPathPrefix(aR,&rpre);
	}

	// Last part of relative name
	rpre= aR;
	if (rpre=="..")
	{
		sResult= g_FileSystem.GetParentName(sResult);
	}
	else if (rpre != ".")
	{
		size_t slen= sResult.GetLength();
		if (slen<1 || sResult[slen-1]!=FILENAME_SLASH)
			// Glue needed (Some paths - "C:\" etc already have glue)
			sResult += CString(FILENAME_SLASH);
		sResult += rpre;
	}
	return sResult;
}

/************************* PathPrefix and StrippedFileName *************************/
/*
	PathPrefix + StrippedFileName = FileName
	for all file names without a trailing FILENAME_SLASH
*/

CString CFileSystem::GetPathPrefix(const char *ach) const
{
	return PathPrefixFromFullName(ach);
}

CString CFileSystem::GetParentName(const char *ach) const
{
	// Return the file name of the parent directory.
	// Only differs from GetPathPrefix by the absence of terminal FILENAME_SLASH

	CString sPath= GetPathPrefix(ach);
	size_t n= sPath.GetLength();

#ifdef _WIN32
	if (n>2 && sPath[n-1]==FILENAME_SLASH && sPath[n-2]==':')
		return sPath;
#endif
	if (n>0 && sPath[n-1]==FILENAME_SLASH)
		return sPath.Left(n-1);
	else
		return sPath;
}

CString CFileSystem::GetStrippedName(const char *ach) const
{
	// Return the filename itself (=part after terminal FILENAME_SLASH)
	CString sPrefix;
	const char *ach0= NULL;
	while (true)
	{
		ach0= ScanMinimalPathPrefix(ach, &sPrefix);
		if (!ach0 || sPrefix.IsEmpty())
			return CString(ach);
		ach= ach0;
	}
}

CString CFileSystem::GetRelativeName(const char *aFileName, const char *aDirName) const
{
	CString a= GetFullName(aFileName);
	CString b= GetFullName(aDirName);
	return SubtractFullNames(a,b);
}

// Concatenate the directory aDirName with the relative filename aRelName
CString CFileSystem::GetCombinedName(const char *aDirName, const char *aRelName) const
{
	CString sD= GetFullName(ToSystemSlash(aDirName));
	CString sR= ToSystemSlash(aRelName);
	return AddFileNames(sD, sR);
}

/***************************************************************************/

bool CFileSystem::IsRoot(const char *ach) const
{
	// TODO: Extend this to handle share names
	CString fn= GetFullName(ach);
	size_t fnlen= fn.GetLength();
	if (fnlen>=2 && fn[fnlen-1]==FILENAME_SLASH && fn[fnlen-2]==':')
		return true;
	return false;
}

bool CFileSystem::IsAbsPath(const char *ach) const
{
	if (ach==NULL || ach[0]==0)
	{
		ASSERT(false);
		return false;
	}

	// Share name --> Absolute
	if (ScanSharePrefix(ach)!=NULL)
		return true;

	// First char FILENAME_SLASH --> Absolute
	if (ach[0]==FILENAME_SLASH)
		return true;
#ifdef _WIN32
	// Windows alias followed by :\ --> Absolute
	size_t len= strlen(ach);
	for (size_t i=0;i+1<len;i++)
		if (ach[i]==':' && ach[i+1]=='\\')
			return true;
#endif
	return false;
}

/*
// What is this (old) code supposed to do? BN 2009-02-10
bool CFileSystem::IsSimpleFileName(const char *ach) const
{
	CString fn= GetFullName(ach);
	int fnlen= fn.GetLength();
	for (int i= 0; i<fnlen-1; i++)
		if (fn[i]=='\\')
			return true;

	return false;
}
*/

/***************************************************************************/

/*
bool IsNormalChar(char c)
{
	return c >= 'a' && c <= 'z' ||
		    c >= 'A' && c <= 'Z' ||
			 c >= '0' && c <= '9';
}
*/

CString CFileSystem::GetNameWithoutExtension(const char *ach) const
{
	if (!ach)
	{
		ASSERT(false);
		return "";
	}

	size_t len= strlen(ach);
	for (size_t i=len;i>0 && ach[--i]!=FILENAME_SLASH;)
	{
		if (ach[i]=='.')
			return CString(ach, i); // strip extension
	}
	return ach; // no extension, return unmodified

	/*
	int len= strlen(ach);
	if (len>2 && IsNormalChar(ach[len-1]))
	{
		int i1= len-1;
		while (i1>2 && IsNormalChar(ach[i1-1]))
			i1--;

		if (ach[i1-1]=='.')
		{
			// Extension - strip
			return CString(ach,i1-1);
		}
	}
	// No extension - return unmodified
	return ach;
	*/
}

CString CFileSystem::GetFileExtension(const char *ach) const
{
	if (!ach)
	{
		ASSERT(false);
		return "";
	}

	size_t len= strlen(ach);
	for (size_t i=len;i>0 && ach[--i]!=FILENAME_SLASH;)
	{
		if (ach[i]=='.')
			return CString(ach+i+1); // return extension
	}
	return ""; // no extension, return empty string

	/*
	int len= strlen(ach);
	if (len>2 && IsNormalChar(ach[len-1]))
	{
		int i1= len-1;
		while (i1>2 && IsNormalChar(ach[i1-1]))
			i1--;

		if (ach[i1-1]=='.')
		{
			// Extension
			return CString(ach+i1,len-i1);
		}
	}
	// No extension - return empty string
	return "";
	*/
}

/***********************************************************************************/

#ifdef _WIN32
bool CFileSystem::CreateDirRecursive(const char *ach)
{
	// Empty strings and null pointers reflect unintended use
	if (!ach || !ach[0])
	{
		ASSERT(false);
		return false;
	}

	// Recursive directory creation
	// If already exists, the function will succeed
	// GetLastError tells the cause of failure
	int dwAttr= GetAttributes(ach);
	if (dwAttr==-1)
	{
		int dwErr= ::GetLastError();

		if (dwErr!=ERROR_FILE_NOT_FOUND &&
			 dwErr!=ERROR_PATH_NOT_FOUND)
			 // General error
			 return false;

		if (!IsRoot(ach))
		{
			CString sParent= g_FileSystem.GetParentName(ach);
			if (!g_FileSystem.CreateDirRecursive(sParent))
				return false; // fail
		}
		return ::CreateDirectory(ach,NULL) == TRUE;
	}
	else if (dwAttr & FILE_ATTRIBUTE_DIRECTORY)
	{
		// Already exists
		return true;
	}
	else
	{
		// Not a directory
		::SetLastError(ERROR_BAD_PATHNAME);
		return false;
	}
}
#endif

/***********************************************************************************/

CLoadBuf::CLoadBuf()
{
	m_aLoadBuf= 0;
	m_BufSize= 0;
	m_bUsed= false;
}

CLoadBuf::~CLoadBuf()
{
	if (m_aLoadBuf)
		free(m_aLoadBuf);
}

bool CLoadBuf::SetSize(size_t size)
{
	if (!m_aLoadBuf)
	{
		if ((m_aLoadBuf= (char *) malloc(size)))
			m_BufSize= size;
		return m_aLoadBuf != 0;
	}
	if (char *a2= (char *) realloc(m_aLoadBuf,size))
	{
		m_aLoadBuf= a2;
		m_BufSize= size;
		return true;
	}
	return false;
}

#ifdef _WIN32
bool CLoadBuf::QuickLoad(const char *aFileName)
{
	long long nSize;
	if (!g_FileSystem.GetFileSize(aFileName, &nSize))
		return false;
	m_nFileSize= nSize;

	HANDLE hFile= CreateFile(aFileName,
		GENERIC_READ,
		FILE_SHARE_READ,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL);

	if (!hFile)
		return false;

#ifndef _WIN64
	if (nSize>=0x100000000-2)
	{
		printf("Error: File size > 32 bits, cannot load file: %s\n", aFileName);
		return false; // cannot load >32 bit-sized files on Win32, only on Win64
	}
#endif
	
	DWORD nBuffer= 100000000;
	bool bOk= false;
	{
		bOk= SetSize(size_t(nSize) + 2); // extra bytes to make sure there is a terminating null (end-of-file in 8- and 16-bit text files)
		if (bOk)
		{
			char *a_buf= m_aLoadBuf;
			while (bOk && nSize>0)
			{
				DWORD nRequested= DWORD(__min(nBuffer, nSize));
				DWORD nRead;
				bOk= false;
				if (ReadFile(hFile, a_buf, nRequested, &nRead, NULL))
				{
					if(nRead==nRequested)
					{
						nSize -= nRead;
						a_buf += nRead;
						bOk= true;
					}
				}
			}

			if (bOk)
			{
				a_buf[0]= 0;
				a_buf[1]= 0;
			}
			else 
				SetSize(0);
		}
	}

	CloseHandle(hFile);
	return bOk;
}
#endif

//#define LOW_LEVEL_IO
#ifdef LOW_LEVEL_IO
#include <io.h>
#include <fcntl.h>
#endif
// Load file into buffer and null terminate
bool CLoadBuf::Load(const char *aFileName)
{
#ifdef _WIN32
	// Use WinAPI directly for quicker loads
	return QuickLoad(aFileName);
#endif

	// Use old code on other systems.
#ifdef LOW_LEVEL_IO
	int hf;
#else
	FILE *pf;
#endif

	size_t nread;
	size_t quant= 1000000;
	size_t accsize= 0;

#ifdef LOW_LEVEL_IO
  if ((hf= _open( aFileName, _O_RDONLY | _O_BINARY)))
#else
	if ((pf= fopen(aFileName,"rb")))
#endif
	{
		if (m_BufSize < accsize+quant)
			if (!SetSize(accsize+quant))
			{
				g_FileSystem.SetLastError("Out of memory");
				return false;
			}

		if (quant < m_BufSize)
			quant= m_BufSize; // Read large chunk

		char *aBuf= m_aLoadBuf;
#ifdef LOW_LEVEL_IO
		int _read( int handle, void *buffer, unsigned int count );
		while (( nread= _read(hf, aBuf, quant) ) > 0)
#else
		while (( nread= fread(aBuf, 1, quant, pf) ) > 0)
#endif
		{
			accsize += nread;
			aBuf += nread;

			if (m_BufSize < accsize+quant)
			{
				if (!SetSize(accsize+quant))
				{
					g_FileSystem.SetLastError("Out of memory");
					return false;
				}
				aBuf= m_aLoadBuf+accsize;
			}
		}

		quant += (quant>>1);
		aBuf[0]= 0; // null terminator

#ifdef LOW_LEVEL_IO
		_close(hf);
#else
		fclose(pf);
#endif
		m_nFileSize= accsize;
		return true;
	}
	// Copy libc error string into FileSystem
	g_FileSystem.SetLastError(strerror(errno));
	m_nFileSize= 0;
	return false;
}


/***********************************************************************************/
/*
void CFileSystem::Job()
{
	int nExitCode= ERROR_SUCCESS;

	// Post message about job termination
	m_bIsRunning= false;
	m_pView->PostMessage(WM_JOB_TERMINATED,(WPARAM) nExitCode,0);
}

UINT CFileSystem_JobThread( LPVOID pParam )
{
	//Controlling function
	CFileSystem *pS= (CFileSystem *) pParam;
	pS->Job();
	return 0; // means OK
}

bool CFileSystem::StartJob()
{
	if (m_pJobThread= AfxBeginThread(CFileSystem_JobThread, (LPVOID) this,
		THREAD_PRIORITY_NORMAL,
		0, // stacksize,
		CREATE_SUSPENDED// dwCreateFlags
		))
	{
		m_pJobThread->m_bAutoDelete= false;
		m_pJobThread->ResumeThread();

		m_bIsRunning= true;
		return true; // success
	}
	return false; // failure
}

void CFileSystem::AddCmd(char *aCmd)
{
	CSingleLock access_lock(&m_cmd_access_mutex, false);
	(void) access_lock.Lock();
	{
		m_cmd_list.push_back(string(aCmd));

		// Set cmd_available
		(void) m_cmd_available_evt.SetEvent();
	}
	access_lock.Unlock();
}
*/

CFileSystem::CFileSystem()
/*
  m_ListAccessMutex(FALSE, // Not initially own
		NULL, // No name
		NULL), // No attrs
  m_FileDoneEvent(FALSE,   // Not initially own
		TRUE,				   // Manually reset
		NULL,             // No Name
		NULL)				   // No security attrs

  */
{
	m_LoadCnt= 0;
}

CFileSystem::~CFileSystem()
{
	// Delete buffer objects
	for (int i= 0; i<m_vpLoadBuf.GetSize(); i++)
	{
		ASSERT(!m_vpLoadBuf[i]->m_bUsed);
		delete m_vpLoadBuf[i];
	}

	//	if (m_pJobThread) delete m_pJobThread;
}

/*******************************************************************************/

CLoadBuf *CFileSystem::AllocBuf()
{
	// Find a free buffer object
	for (int i= 0; i<m_vpLoadBuf.GetSize(); i++)
	{
		if (!m_vpLoadBuf[i]->m_bUsed)
		{
			m_vpLoadBuf[i]->m_bUsed= true;
			return m_vpLoadBuf[i];
		}
	}

	// No free buffer objects due to recursive load (or first time)
	CLoadBuf *pNew= NULL;
	try
	{
		pNew= new CLoadBuf;
	}
	catch(...) { } // new CLoadBuf failed

	if (pNew)
	{
		pNew->m_bUsed= true;
		try
		{
			m_vpLoadBuf.Add(pNew);
			return pNew;
		}
		catch(...) {} // Memory allocation failure during Add

		delete pNew;
	}
	return NULL;
}

void CFileSystem::FreeBuf(CLoadBuf *pLoadBuf)
{
	// Free load buffer object after use

	ASSERT(pLoadBuf != NULL && pLoadBuf->m_bUsed); // Encourage the programmer to be thorough
	if (!pLoadBuf) // ...but be robust in Release mode
		return;

	// New style: Do real unalloc
	int i=0;
	for (;i<m_vpLoadBuf.GetSize();i++)
		if (m_vpLoadBuf[i]==pLoadBuf)
			break;
	if (i<m_vpLoadBuf.GetSize())
	{
		delete m_vpLoadBuf[i];
		m_vpLoadBuf.Delete(i);
	}

	/*
	// Old style: Keep the buffer to avoid the penalty of constantly
	// allocating a new buffer. Just marked as unused. This modum vitae
	// appropriate in Schematic Studio, but not here.
	if (pLoadBuf)
		pLoadBuf->m_bUsed= false;
	*/
}

CLoadBuf *CFileSystem::Load(const char *aFileName)
{
	// Use Free(...) to release the buffer after use
	if (CLoadBuf *pLoadBuf= AllocBuf())
	{
		if (pLoadBuf->Load(aFileName))
		{
			// Keep track of the number of files that have been loaded
			// for performance analysis
			m_LoadCnt++;
			return pLoadBuf;
		}
		FreeBuf(pLoadBuf);
	}
	else
		SetLastError("Out of memory");
	return NULL;
}

bool CFileSystem::IsEqual(const char *aFileName, const char *aData)
{
	// Compare file contents to data string, return true if equal
	bool bEqual= false;
	if (CLoadBuf *pLB= Load(aFileName))
	{
		if (const char *aCh= GetFileData(pLB))
			bEqual= !strcmp(aCh,aData);

		FreeBuf(pLB);
	}
	// else possibly nonexistent file - treat as different

	return bEqual;
}

char *CFileSystem::GetFileData(CLoadBuf *pLoadBuf)
{
	return pLoadBuf->m_aLoadBuf;
}

/***********************************************************************************/

CString CFileSystem::GetCurrentDir() const
{
	char aBuf[STRBUFSIZE];
#ifdef _WIN32
	int n= GetCurrentDirectory(STRBUFSIZE,aBuf); //Win32 API
	if (n>=0 && n<STRBUFSIZE)	// Fit inside buffer?
		return CString(aBuf,n);
#else
	if (getcwd(&aBuf[0], STRBUFSIZE)!=0)
		return CString(aBuf);
#endif
	return CString("");
}

#ifdef _WIN32
void CFileSystem::SetCurrentDir(const char *aDir)
{
	SetCurrentDirectory(aDir);
}
#endif

void CFileSystem::SetExecutableFileName(const char *ach)
{
	m_sExeName= GetCombinedName(GetCurrentDir(), ach);
	/*
#ifdef _WIN32
	m_sExeName= GetFullName(ach);
#else
	// TODO: Allow use of GetFullName on both Windows and Unix
	char aBuf[STRBUFSIZE];
	if (getcwd(aBuf, STRBUFSIZE)==0)
		m_sExeName= "";
	else
		m_sExeName= GetCombinedName(aBuf, argv[0]);
#endif
	*/
}

CString CFileSystem::GetExecutableFileName() const
{
#ifdef _WIN32
	if (m_sExeName.IsEmpty())
	{
		// Windows
		if (HMODULE h= GetModuleHandle(NULL))
		{
			char aBuf[STRBUFSIZE];
			int n= ::GetModuleFileName(h,aBuf,STRBUFSIZE);

			// Success and fit inside buffer?
			if (n>0 && n < STRBUFSIZE)
				return CString(aBuf,n);
		}
	}
#endif
	return m_sExeName;
}

/*******************************************************************************/

CString CFileSystem::GetUserDir() const
{
#ifdef _WIN32
	// Windows
	return Environment::GetEnv("USERPROFILE");
#else
	// Unix
	return "~";
#endif
}

/*******************************************************************************/

bool CFileSystem::CreateBinary(const char *aFileName, const char *aBinary, size_t nBinary)
{
	// Binary file I/O is always used for maximum efficency
	ASSERT(nBinary 	>= 0);
	bool bOk= false;

	// No filename or empty file name? (fopen crashes otherwise)
	if (!aFileName || !aFileName[0])
		return false;

	// Create an empty file on disk
	if (FILE *pf= fopen(aFileName,"wb"))
	{
		bOk= true;
		if (nBinary)
		{	
			size_t nWritten= fwrite(aBinary, 1, nBinary, pf);
			if (nWritten != nBinary)
				bOk= false;
		}
		fclose(pf);
	}
	return bOk;
}

bool CFileSystem::CreateBinary(const char *aFileName, const char *aStr)
{
	return CreateBinary(aFileName,aStr,strlen(aStr));
}

/***********************************************************************************/

bool CFileSystem::FileExists(const char *ach) const
{
	// Works on both Windows and Unix
	return GetAttributes(ach)>=0;
}

#ifdef _WIN32
bool CFileSystem::IsProtected(const char *ach) const
{
	// Return if the file is protected for deletion or writing
	int dwAttr= GetAttributes(ach);

	if (dwAttr==-1)
		return false; // File does not exist

	if (dwAttr & FILE_ATTRIBUTE_READONLY)
		return true;
	return false;
}
#endif

/*
#ifdef _WIN32
CTime CFileSystem::GetModifiedTime(const char *aFileName) const
{
	// Get the time when the file was last modified
	WIN32_FIND_DATA findFileData;
	HANDLE hFind = FindFirstFile((LPTSTR)aFileName, &findFileData);
	if (hFind == INVALID_HANDLE_VALUE)
		return CTime(0);
	FindClose(hFind);
	return CTime(findFileData.ftLastWriteTime);
}
#else
#error TODO: Implement CFileSystem::GetModifiedTime under Unix
#endif
*/

bool CFileSystem::ListFiles(const char *aDir, const char *aFilter, CVector<CString> &v_out)
{
	// aDir= Foldername
	// aFilter= Stripped filename filter, possibly with * and ? wildcards
	v_out.SetSize(0);
	CString sDir;
	if (aDir && aDir[0]!=0)
		sDir= aDir;
	else
		sDir= GetCurrentDir();

#ifdef _WIN32
	// Windows
	_finddata_t f;
	intptr_t hFile= _findfirst(GetCombinedName(sDir, "*.*"), &f);
	if (hFile == -1)
	{
		//printf("sDir=='%s'\n", sDir);
		//printf("hFile==-1\n");
		return false; // Avoid _findclose if failure (dirs are never empty, at least . and .. should be there)
	}
	do
	{
		//printf("%s\n", f.name);
		if ((f.attrib & _A_SUBDIR)==0 && strmatch_nocase(f.name, aFilter))
			v_out.Add(g_FileSystem.GetCombinedName(aDir, f.name));
     } while( _findnext( hFile, &f ) == 0 );
	_findclose( hFile );
#else
	// Unix
    DIR *dirp= opendir(sDir);
	if (!dirp)
		return false;
    struct dirent *dp;
	while ((dp= readdir(dirp)))
	{
		CString sFileName= dp->d_name;
		int nAttr= GetAttributes(g_FileSystem.GetCombinedName(aDir, sFileName));
        if ((nAttr & S_IFREG)!=0 && strmatch(sFileName, aFilter))
			v_out.Add(g_FileSystem.GetCombinedName(aDir, sFileName));
    }
    (void) closedir(dirp);
#endif
    return true;
}

bool CFileSystem::ListSubdirectories(const char *aDir, const char *aFilter, CVector<CString> &v_out)
{
	// aDir= Foldername
	// aFilter= Stripped filename filter, possibly with * and ? wildcards
	v_out.SetSize(0);
	CString sDir;
	if (aDir && aDir[0]!=0)
		sDir= aDir;
	else
		sDir= GetCurrentDir();

#ifdef _WIN32
	// Windows
	_finddata_t f;
	intptr_t hFile= _findfirst(GetCombinedName(sDir, "*.*"), &f);
	if (hFile==-1)
		return false; // Avoid _findclose if failure (dirs are never empty, at least . and .. should be there)
	do
	{
		if ((f.attrib & _A_SUBDIR)!=0 && strmatch_nocase(f.name, aFilter) && strcmp(".", f.name)!=0 && strcmp("..", f.name))
			v_out.Add(g_FileSystem.GetCombinedName(aDir, f.name));
     } while( _findnext( hFile, &f ) == 0 );
	_findclose( hFile );
#else
	// Unix
    DIR *dirp= opendir(sDir);
	if (!dirp)
		return false;
    struct dirent *dp;
	while ((dp= readdir(dirp)))
	{
		CString sFileName= dp->d_name;
		int nAttr= GetAttributes(g_FileSystem.GetCombinedName(aDir, sFileName));
        if ((nAttr & S_IFDIR)!=0 && strmatch(sFileName, aFilter) && strcmp(".", sFileName)!=0 && strcmp("..", sFileName))
			v_out.Add(g_FileSystem.GetCombinedName(aDir, sFileName));
    }
    (void) closedir(dirp);
#endif
    return true;
}

/*
// NOTE: This currently returns full names for all the files
// Relative file names would seem tighter
// NOTE: This uses Windows pattern matching where a.scl1 matches *scl
// - weird! so we should rather use a different pattern matching method
#include <system/stdout.h>
CVector<CString> CFileSystem::ListFiles(CString sDir, CString sPattern,
													 bool bDirs,
													 bool bRecursive)
{
	CVector<CString> vsList;

	CString sFilePattern= GetCombinedName(sDir,sPattern);

	WIN32_FIND_DATA fd;
	HANDLE hSearch = FindFirstFile(sFilePattern, &fd);
	if (hSearch != INVALID_HANDLE_VALUE)
	{
		bool bMore= true;
		while (bMore)
		{
			const char *ach= fd.cFileName;

			if (!strcmp(ach,".") || !strcmp(ach,".."))
			{
				// ignore Windows dummy dirs "." and ".."
			}
			else
			{
				bool bIsDir= (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
				if (bDirs == bIsDir)
					vsList.Add( GetCombinedName(sDir, ach) ); // caller wants dirs/files
			}
			bMore= (FindNextFile(hSearch, &fd) != 0);
		}
	}

	// Close the search handle.
	if (hSearch != INVALID_HANDLE_VALUE)
		(void) FindClose(hSearch);

	// Breadth-first recursion into subdirs
	if (bRecursive)
	{
		// Get the subdirs (non-recursive without pattern-matching)
		CVector<CString> vsDir= ListFiles(sDir,"*",true,false);

		for (int i= 0; i<vsDir.GetSize(); i++)
			vsList.Append( ListFiles(vsDir[i], sPattern, bDirs, bRecursive ));
	}

	return vsList;
}
*/

#ifdef MFC_CLASSES
bool CFileSystem::SaveBitmap(HDC hDC, CBitmap &bm, CString sFileName)
{
	BITMAP b;
	int nRes= bm.GetObject(sizeof(BITMAP), &b);

	bool bOk= false;
	try
	{
		BITMAPINFOHEADER h;
		h.biSize= sizeof(BITMAPINFOHEADER);
		h.biWidth= b.bmWidth;
		h.biHeight= b.bmHeight;
		h.biPlanes= 1;
		h.biBitCount= 24;
		h.biCompression= BI_RGB;
		h.biSizeImage= 0;
		h.biXPelsPerMeter= 3000;
		h.biYPelsPerMeter= 3000;
		h.biClrUsed= 0;
		h.biClrImportant= 0;

		BITMAPINFO bi;
		bi.bmiHeader= h;

		int nSize= h.biHeight*h.biWidth*3;
		int nFileSize= nSize + sizeof(BITMAPINFOHEADER) + 14;
		char *pBin= new char[nFileSize];
		pBin[0]= 'B';
		pBin[1]= 'M';

		BITMAPINFOHEADER *pHdr= (BITMAPINFOHEADER *)(pBin + 14);
		*pHdr= bi.bmiHeader;

		char *pBits= pBin + sizeof(BITMAPINFOHEADER) + 14;

		int nLines= GetDIBits(hDC,
			(HBITMAP)bm,
			0,
			b.bmHeight,
			pBits,
			&bi,
			DIB_RGB_COLORS);

		long *pInt= (long *)(pBin+2);
		*pInt++ = nFileSize;
		*pInt++ = 0;
		*pInt++ = pBits-pBin;

		bOk= g_FileSystem.CreateBinary(sFileName, pBin, nFileSize);

		delete [] pBin;
	}
	catch (...)
	{
	}

	return bOk;
}
#endif

#ifdef _WIN32
bool CFileSystem::GetFileSize(const char *aFilename, long long *p_sz)
{
	if (aFilename==NULL || *aFilename==0 || p_sz==NULL)
		return false;

	HANDLE hFile= CreateFile(aFilename,
		GENERIC_READ,
		FILE_SHARE_READ,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	if (!hFile)
		return false;

	LARGE_INTEGER sz;
	bool bOk= GetFileSizeEx(hFile, &sz)!=0;
	*p_sz= (long long (sz.HighPart) << 32) | long long (sz.LowPart);

	(void) CloseHandle(hFile);
	return bOk;
}
#endif

int CFileSystem::GetAttributes(const char *aFilename) const
{
	// Returns attr or -1 for both Windows and Unix
	// NOTE: Meaning of returned value is different though!
#ifdef _WIN32
	return ::GetFileAttributes(aFilename);
#else
	struct stat s;
	if (stat(aFilename, &s)==0)
		return s.st_mode;
	return -1;
#endif
}

CString CFileSystem::ToSystemSlash(const char *ach) const
{
	CString s_out= ach;
#ifdef _WIN32
	for (size_t i=0;i<s_out.GetLength();i++)
		if (s_out[i]=='/' || s_out[i]=='\\')
			s_out.SetAt(i, FILENAME_SLASH);
#endif
	return s_out;
}

CString CFileSystem::GetTmpFilename(const char *aDir) const
{
	// Invents a temporary file name in a given dir
	CString sName;
	do
	{
		sName= GetCombinedName(aDir, ::Format("__%d%d__.tmp", rand(), rand()));
	}
	while (FileExists(sName));
	return sName;
}

CString CFileSystem::GetSafeStrippedName(CString s)
{
	s= GetStrippedName(s);
	for (size_t i=0;i<s.GetLength();i++)
	{
#ifdef _WIN32
		if (s[i]==FILENAME_SLASH || s[i]=='/' || s[i]==' ' || s[i]=='$')
#else
		if (s[i]==FILENAME_SLASH || s[i]==' ')
#endif
		{
			s.SetAt(i, '_');
		}
	}
	return s;
}

#ifdef _WIN32
bool CFileSystem::Delete(const char *aFilename)
{
	return ::DeleteFile(aFilename)!=0;
}
#endif

#ifdef _WIN32	
bool CFileSystem::GetFileTime(const char *aFilename, FILETIME *lpCreate, FILETIME *lpAccess, FILETIME *lpWrite)
{
	// WINDOWS
	if (!aFilename)
		return false;
	HANDLE hFile= CreateFile(aFilename,
		GENERIC_READ, 
		FILE_SHARE_READ, 
		NULL, 
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	if (!hFile)
		return false;
	bool bOk= ::GetFileTime(hFile, lpCreate, lpAccess, lpWrite)!=0;
	CloseHandle(hFile);
	return bOk;
}
#else
bool CFileSystem::GetFileTime(const char *aFilename, time_t &atime, time_t &mtime)
{
	// POSIX
	if (!aFilename)
		return false;
	struct stat sb;
	if (stat(aFilename, &sb)!=0)
		return false;
	atime= sb.st_atime;
	mtime= sb.st_mtime;
	return true;
}
#endif

#ifdef _WIN32
bool CFileSystem::SetFileTime(const char *aFilename, FILETIME *lpCreate, FILETIME *lpAccess, FILETIME *lpWrite)
{
	// WINDOWS
	if (!aFilename)
		return false;
	HANDLE hFile= CreateFile(aFilename, 
		GENERIC_WRITE,
		FILE_SHARE_READ, 
		NULL, 
		OPEN_EXISTING,
		0,
		NULL);
	if (!hFile)
		return false;	
	bool bOk= ::SetFileTime(hFile, lpCreate, lpAccess, lpWrite)!=0;

	CloseHandle(hFile);
	return bOk;
}
#else
bool CFileSystem::SetFileTime(const char *aFilename, const time_t atime, const time_t mtime)
{
	// POSIX
	if (!aFilename)
		return false;
	
	utimbuf tb;
	tb.actime= atime;
	tb.modtime= mtime;
	return (utime(aFilename, &tb)==0);
}
#endif

#ifdef _WIN32
bool timecmp(const FILETIME *lpT0, const FILETIME *lpT1)
{
	return (lpT0->dwHighDateTime==lpT1->dwHighDateTime && lpT0->dwLowDateTime==lpT1->dwLowDateTime);
}
#endif

#ifdef _WIN32
// WINDOWS
bool CFileSystem::HasFileTime(const CString &sFilename, const FILETIME *lpCreate, const FILETIME *lpAccess, const FILETIME *lpWrite)
{
	// Returns true if sFilename has specific time stamps (NULL pointers mean the field won't be compared)
	if (!g_FileSystem.FileExists(sFilename))
		return false;
	FILETIME tCreate, tAccess, tWrite;
	if (!g_FileSystem.GetFileTime(sFilename, &tCreate, &tAccess, &tWrite))
		return false;
	if (lpCreate!=NULL && !timecmp(&tCreate, lpCreate))
		return false;
	if (lpAccess!=NULL && !timecmp(&tAccess, lpAccess))
		return false;
	if (lpWrite!=NULL && !timecmp(&tWrite, lpWrite))
		return false;
	return true;
}
#else
// POSIX
bool CFileSystem::HasFileTime(const char *aFilename, const time_t *p_atime, const time_t *p_mtime)
{
	if (!g_FileSystem.FileExists(aFilename))
		return false;

	time_t atime, mtime;
	if (!g_FileSystem.GetFileTime(aFilename, atime, mtime))
		return false;

	if (p_atime!=0 && *p_atime!=atime)
		return false;
	if (p_mtime!=0 && *p_mtime!=mtime)
		return false;

	return true;
}
#endif

#ifdef _WIN32
int CFileSystem::CompareFiles(const char *aFile0, const char *aFile1, CString *psError)
{
	// returns:
	//		0 if files differ in content
	//		1 if files equal
	//		-1 if failed (sets psError in that case)

	long long size0, size1;
	if (!g_FileSystem.GetFileSize(aFile0, &size0))
	{
		if (psError)
			*psError= ::Format("Cannot open '%s'", aFile0);
		return -1;
	}
	if (!g_FileSystem.GetFileSize(aFile1, &size1))
	{
		if (psError)
			*psError= ::Format("Cannot open '%s'", aFile1);
		return -1;
	}
	if (size0!=size1)
		return 0;

	HANDLE h0= CreateFile(aFile0,
		GENERIC_READ,
		FILE_SHARE_READ,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	if (!h0)
	{
		if (psError)
			*psError= ::Format("Cannot open '%s'", aFile0);
		return -1;
	}

	HANDLE h1= CreateFile(aFile1,
		GENERIC_READ,
		FILE_SHARE_READ,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	if (!h1)
	{
		if (psError)
			*psError= ::Format("Cannot open '%s'", aFile1);
		return -1;
	}

	const DWORD nBuf= 10000000;
	CVector<char> vBuf0(nBuf), vBuf1(nBuf);
	char *pBuf0= vBuf0.GetBuffer();
	char *pBuf1= vBuf1.GetBuffer();

	bool bOk= true;
	bool bEqual= true;
	while (bOk && bEqual && size0>0)
	{
		DWORD nRequested= DWORD(__min(nBuf, size0));
		DWORD nRead0,  nRead1;
		ReadFile(h0, pBuf0, nRequested, &nRead0, NULL);
		ReadFile(h1, pBuf1, nRequested, &nRead1, NULL);
		if (nRead0!=nRead1)
		{
			if (psError)
				*psError= "ReadFile error";
			bOk= false;
		}

		if (bOk)
		{
			for (DWORD i=0;bEqual && i<nRead0;i++)
				bEqual= pBuf0[i]==pBuf1[i];
			size0 -= nRead0;
		}
	}

	CloseHandle(h0);
	CloseHandle(h1);

	if (!bOk)
		return -1;
	if (bEqual)
		return 1;
	else
		return 0;
}

bool CFileSystem::CopyFile(const CString &sSrc, const CString &sDst, CString *psError)
{
	// copies source to destination, sets the destination's create and write times to the source's

	bool bOk= true;

	// open src
	HANDLE h_src= CreateFile(sSrc,
		GENERIC_READ,
		FILE_SHARE_READ,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	if (!h_src)
	{
		if (psError)
			*psError= ::Format("Cannot open '%s'", sSrc);
		return false;
	}
	long long nSize_src;
	if (!g_FileSystem.GetFileSize(sSrc, &nSize_src))
	{
		if (psError)
			*psError= ::Format("Cannot determine size of '%s'", sSrc);
		return false;
	}

	// open dst
	HANDLE h_dst= CreateFile(sDst, 
		GENERIC_WRITE,
		FILE_SHARE_WRITE,
		NULL,
		CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	if (!h_dst)
	{
		if (psError)
			*psError= ::Format("Cannot create '%s'", sDst);
		CloseHandle(h_src);
		return false;
	}

	const DWORD nBufSize= 10000000;
	CVector<char> vBuf(nBufSize);
	char *pBuf= vBuf.GetBuffer();

	while (bOk && nSize_src>0)
	{
		DWORD nRequested= DWORD(__min(nSize_src, nBufSize));
		DWORD nRead;
		bOk= ReadFile(h_src, pBuf, nRequested, &nRead, NULL)!=0;
		if (bOk)
		{
			if (nRequested==nRead)
			{
				DWORD nWritten;
				bOk= WriteFile(h_dst, pBuf, nRead, &nWritten, NULL)!=0;
				if (bOk)
				{
					if (nWritten==nRead)
					{
						nSize_src -= nRead;
					}
					else
					{
						if (psError)
							*psError= "WriteFile error";
						bOk= false;
					}
				}		
			}
			else
			{
				if (psError)
					*psError= "ReadFile error";
				bOk= false;
			}
		}
	}

	CloseHandle(h_src);
	CloseHandle(h_dst);

	// set time
	if (bOk)
	{
		FILETIME tWrite_src;
		bOk= g_FileSystem.GetFileTime(sSrc, NULL, NULL, &tWrite_src);
		if (bOk)
			bOk= g_FileSystem.SetFileTime(sDst, NULL, NULL, &tWrite_src);
	}

	return bOk;
}

#endif
