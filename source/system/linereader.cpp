//
// linereader.cpp --
//

#include "linereader.h"
#include "system/filesystem.h"

CLineReader::CLineReader()
{
	m_pf= 0;
}

CLineReader::~CLineReader()
{
	if (m_pf)
		(void) Close();
}

bool CLineReader::Open(const char *aFilename)
{
	if (m_pf)
		(void) Close();

	m_pf= fopen(aFilename, "rb");
	if (!m_pf)
		return false;
	m_vBuf.SetSize(0);
	m_aBuf= NULL;
	m_nRead= 0;
	m_nPos= 0;

#ifdef _WIN32
	m_nRead_sum= 0;
	if (!g_FileSystem.GetFileSize(aFilename, &m_nRead_filesize))
		return false;
#endif
	
	return true;
}

bool CLineReader::Close()
{
	bool bOk= true;
	if (m_pf)
	{
		bOk= fclose(m_pf)==0;
		m_pf= NULL;
	}
	return bOk;
}

bool CLineReader::GetLine(char *aLineBuf, size_t nBufSize, CVector<const char *> &vItems, int &nItems)
{
	return GetLine(aLineBuf, nBufSize, vItems, nItems, 9, 9, false);
}

bool CLineReader::GetLine(char *aLineBuf, size_t nBufSize, CVector<const char *> &vItems, int &nItems, char sep1, char sep2, bool bMulti)
{
	nItems= 0;

	if (m_vBuf.GetSize()==0)
	{
		// Init buffer on first call
		try
		{
			m_vBuf.SetSize(LINEREADER_QUANTA+1);
		}
		catch(...)
		{
			// Out of memory
			return false;
		}

		m_aBuf= m_vBuf.GetBuffer(m_vBuf.GetSize());
		m_nRead= fread(m_aBuf, 1, LINEREADER_QUANTA, m_pf);
		m_nPos= 0;
	}

	// Skip empty line
	*aLineBuf= 0;
	while ((m_aBuf[m_nPos]==10 || m_aBuf[m_nPos]==13) && m_nRead>0)
	{
		while ((m_aBuf[m_nPos]==10 || m_aBuf[m_nPos]==13) && m_nPos<m_nRead)
		{
			m_nPos++;
			m_nRead_sum++;
		}

		if (m_nPos>=m_nRead)
		{
			m_nRead= fread(m_aBuf, 1, LINEREADER_QUANTA, m_pf);
			m_nPos= 0;
		}
	}
	if (!m_nRead)
		return false; // EOF reached

	int state= -1; // -1= start, 0= tab, 1= char
	while (m_aBuf[m_nPos]!=10 && m_aBuf[m_nPos]!=13 && m_nRead>0)
	{
		const char *ach0= m_aBuf + m_nPos;
		while (m_aBuf[m_nPos]!=10 && m_aBuf[m_nPos]!=13 && m_nPos<m_nRead)
		{
			m_nPos++;
			m_nRead_sum++;
		}
		const char *ach1= m_aBuf + m_nPos;
		
		if (bMulti)
		{
			// treat multiple connected separators as one separator
			while (ach0<ach1 && nBufSize>0)
			{
				char c= *ach0++;
				if (c==sep1 || c==sep2)
				{
					state= 0;
					c= 0;
				}
				else
				{
					if (state<=0)
					{
						vItems.SetAt(nItems, aLineBuf);
						nItems++;
					}
					state= 1;
				}

				*aLineBuf++ = c;
				nBufSize--;
			}
		}
		else
		{
			if (nItems==0)
			{
				vItems.SetAt(0, aLineBuf); // assumes that vItems has been initialized -- TODO: make this safe
				nItems++;
			}
			
			while (ach0<ach1 && nBufSize>0)
			{
				char c= *ach0++;
				if ((c==sep1 || c==sep2) && nItems<vItems.GetSize()) // separator (sep1 and sep2 can be same character)
				{
					*aLineBuf= 0;
					//printf("nitems= %d\n", nItems);
					//printf("vitems= '%s'\n", vItems[nItems-1]);
				
					vItems.SetAt(nItems, aLineBuf+1);  // assumes that vItems has been initialized -- TODO: make this safe

					nItems++;
					c=0; // output a null terminator
				}
				*aLineBuf++ = c;
				nBufSize--;
			}

			/*
			// treat each separator 
			while (ach0<ach1 && nBufSize>0)
			{
				char c= *ach0++;
				if (state<=0 && nItems<vItems.GetSize())
				{
					if (nItems>0)
						*(aLineBuf-1)= 0;
					vItems.SetAt(nItems, aLineBuf);
					nItems++;
				}
				if (c==sep1 || c==sep2)
					state= 0;
				else
					state= 1;

				*aLineBuf++ = c;
				nBufSize--;
			}*/
		}
		*aLineBuf= 0; // final null terminator
		vItems.SetAt(nItems, aLineBuf+1); // extra dummy pointer to first char after last null

		if (m_nPos>=m_nRead)
		{
			m_nRead= fread(m_aBuf, 1, LINEREADER_QUANTA, m_pf);
			m_nPos= 0;
		}
	}

	// trim leading and trailing spaces
	for (int i=0;i<nItems;i++)
	{
		const char *p0= vItems[i];
		if (*p0==32)
		{
			while (*p0==32)
				p0++;
			vItems.SetAt(i, p0);
		}
		char *p1= (char *)vItems[i+1];
		if (p1>p0 && *(p1-1)==0)
			p1--;
		while (p1>p0 && *(p1-1)==32)
			p1--;
		*p1= 0;
	}

	if (!m_nRead)
		return false; // EOF reached
	return true;
}

float CLineReader::GetProgress()
{
	if (m_pf)
		return float(m_nRead_sum)/m_nRead_filesize;
	else
		return 0.0f;
}
