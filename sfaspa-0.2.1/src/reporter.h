/** 
 * @file       reporter.h
 * @brief      SPA assembly output generator.
 * @date       Modified on Thu 2014-01-16 05:33:19 PM
 * @author     Youngik Yang
 * @version    0.02
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __REPOERTER_H__
#define __REPOERTER_H__

#include "recruiter.h"



/**
 * \brief SPA assembly outout generator
 */
class Reporter
{
 private:
	size_t        npaths;
	char          **seqs;
    PathAligners  *Alns;
	NumToPathMap  *Id2Paths;
	PathToNumMap  *Path2Ids;
	FlagArray     *bad_paths;	
	MetaPath      *meta_paths;

 private:
	void write();
	void writePlacement();
	void writeAlignment();
	void writeSequences();
	void combine( std::string name );
	PathLengthMap getPathLengths();

 public:
	Reporter();
	~Reporter();
	void init(Loader *loader, MetaPath *mpath);
	void run();
};

#endif
