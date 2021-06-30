//
// job.h -- Interface for jobs with progress monitoring
//

#ifndef __IJOB_H__
#define __IJOB_H__

#include "types/string.h"
#include "types/vector.h"

class IJob
{
protected:
	unsigned int m_fedebeda;
	IJob() { m_fedebeda= 0xFEDEBEDA; }
	virtual ~IJob() {};

public:
	CVector<CString> m_vFilesCreated;

	// Progress bar
	virtual void SetProgress(float x)= 0;

	// Set text below progress bar
	virtual void SetMessage(const char *ach)= 0;

	// Add text to text box
	virtual void PutText(const char *ach)= 0;

	// Debugging
	bool IsValid() { return m_fedebeda==0xFEDEBEDA; }
};

#endif
