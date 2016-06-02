/**
 * \file      read.h
 * \brief     Read type declaration
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \date      Modified on Tue 2013-12-17 06:25:06 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __READ_H__
#define __READ_H__

#include <algorithm>

typedef unsigned ReadId;
typedef unsigned CoverageType;

class Read 
{
 public:
	CoverageType size;
	ReadId *rid;
        
	Read() 
	{ 
		size = 0; 
	}
	
	Read(CoverageType s) 
    { 
		rid = new ReadId[s]; 
		size = 0; 
	}

	~Read() 
	{ 
		delete[] rid; 
	}

	void add(ReadId &r) 
	{
		rid[size++] = r; 
	} 

	void sort() 
	{
		std::sort( rid, rid+size );
	}
};

#endif
