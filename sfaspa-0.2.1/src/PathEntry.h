/** 
 * @file       PathEntry.h
 * @brief      A single path summary
 * @details    A summary of single path which used clustering and extension during object initialization.
 * @date       Modified on Tue 2013-12-17 06:40:18 PM
 * @author     Youngik Yang
 * @version    0.02
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __PATH_ENTRY_H__
#define __PATH_ENTRY_H__

#include <string>
#include <fstream>
#include <cassert>
//#include <tr1/unordered_map>
#include "path.h"

/**
 * \brief Single path summary which is an entry of path clustering.
 */
struct PathEntry
{
	std::string seq;
	int lstop,rstop;
	int ltrim,rtrim;

	PathEntry();
	PathEntry(std::string &s, int, int, int, int);
	void dump( std::fstream &out );
	void load( std::fstream &in  );
};

typedef std::tr1::unordered_map<PathId, PathEntry> PathEntryMap;

#endif

