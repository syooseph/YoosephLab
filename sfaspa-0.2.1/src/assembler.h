/** 
 * @file       assembler.h
 * @brief      Main controller for SPA assembly stages.
 * @details    Wrapper class for all SPA assembly components.
 * @date       Fri 2014-01-03 01:07:00 PM
 * @date       Modified on Thu 2014-02-27 11:54:52 AM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

//#include <boost/filesystem.hpp>
#include "path.h"
#include "merger.h"
#include "ljoiner.h"
#include "sjoiner.h"
#include "placer.h"
#include "recruiter.h"
#include "reporter.h"
#include "rjoiner.h"
#include "rjoinerSA.h"

/** 
 * Assembly stages
 */
enum stage { EXTRACT, CLUSTER, RECLUSTER, LONGJOIN, SHORTJOIN, READJOIN, PLACE, RECRUIT, REPORT };

/** 
 *  Binary file of each assembly stage for object dump
 */
const std::string BinFiles[] = { 
	"path.bin", 
	"cluster.bin", 
	"metacluster.bin",
	"ljoin.bin", 
	"sjoin.bin", 
	"rjoin.bin",
	"place.bin",
	"recruit.bin"
};


/**
 * \class Assembler
 * \brief Top level controller of assembly modules
 */
class Assembler
{
 private:
	Loader      *loader;  ///< File loader
	Extracter   *finder;  ///< Path finder
	Merger      *merger;  ///< Merger for clustering paths
	Merger      *metacc;
	Connecter   *ljoins;  ///< Latcher for long overlapping paths
	Connecter   *sjoins;  ///< Latcher for short overlapping paths
	Connecter   *rjoins;  ///< Latcher for read bridging extension
	Placer      *placer;  ///< Read placer
	Recruiter   *radder;  ///< Read recruiter
	Reporter    *writer;  ///< Output generator
	
	PathId      *assem_reads;  ///< assembled reads pointer
	PathEntryMap assem_paths;  ///< assembled paths <PathId,Sequence>
	
 private:
	/**
	 * Retrieve file name of given assembly stage 
	 * \param stage assembly stage
	 */
	std::string getFileName(int stage);

	//============
	// major steps
	//============
	void load();     ///< Load files (reads, suffix arrays, graph input)
	void extract();  ///< Extract paths
	void cluster();  ///< Cluster paths
	void recluster();  ///< Cluster paths
	void connect();  ///< Extend paths
	void place();    ///< Place reads to paths
	void recruit();  ///< Recruit unassigned reads
	void report();   ///< Generate outputs

	//==================================
	// Specific types of path extensions
	//==================================
    void long_overlap_extend();  ///< Long overlap path extension
    void tiny_overlap_extend();  ///< Short overlap path extension
    void read_straddle_extend(); ///< Read briding path extendsion

	//===============
	// Reclaim memory
	//===============
	void clear();  ///< Clear objects
	void reset();  ///< Reset pointers

	//==========================
	// Dependent objects loaders
	//==========================
	void checkDependency(int stage); ///< Load dependent objects of given stage
	void reloadExtracter();          ///< Reload initial paths
	void initiateExtracter();        ///< Load initial paths
	void initiateMerger();           ///< Load clusters
	void initiateLongJoiner();       ///< Load long overlapping paths
	void initiateTinyJoiner();       ///< Load short overlapping paths
	void initiateReadJoiner();       ///< Load extended paths with bridging reads
	void initiatePlacer();           ///< Load read placements
	void initiateRecruiter();        ///< Load read recruitment
	void loadAssembledReads();       ///< Load assembled reads from initial paths

 public:
	Assembler(Param&); ///< constructor
	~Assembler();      ///< desctructor

	void run(); ///< Run assembler
};

#endif
