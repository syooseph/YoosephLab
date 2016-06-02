/** 
 * \file      connecter.h
 * \brief     Base class of Path extension
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \date      Modified on Fri 2014-02-28 03:43:19 PM
 * \bug       None.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __CONNECTER_H__
#define __CONNECTER_H__

#include "extender.h"
#include "merger.h"
#include "esectioner.h"
#include "log.h"
#include "ConnectEntries.h"

/** Kmer to PathId mapping */
typedef std::tr1::unordered_map<KmerId, std::set<PathId> > KmerToPathMap;

/** Sequence to PathId mapping */
typedef std::tr1::unordered_map<std::string, std::set<PathId> > StringToPathMap;

/** overlap size and StringToPathMap (head/tail sequence to PathId mapping) */
typedef std::tr1::unordered_map<size_t, StringToPathMap> EndStringMap;

/** Read flags */
typedef std::tr1::unordered_map<ReadId, bool> ReadFlagMap;

/**
 * Comparison of two entries.
 * Priority: size of k-mers
 */
inline bool cmp_entry( const JoinEntry &e1, const JoinEntry &e2 )
{
    if ( e1.nkmers >  e2.nkmers ) return true;    
	if ( e1.nkmers == e2.nkmers ) {
		if ( e1.seqlen > e2.seqlen ) return true;
		else if ( e1.seqlen < e2.seqlen ) return false; 
		else return e1.pid > e2.pid;
	}
    return false;
}

/** Types of path extensions */
enum LATCH_TYPE { LONG_OVERLAP_EXTEND, SHORT_OVERLAP_EXTEND, READ_BRIDGE_EXTEND };


/**
 * \class Connecter
 * \brief Path Extension
 */
class Connecter
{
 protected:
	PathIdSet     merged_paths;  ///< merged paths during path extension
	ExtenderMap   extenders;     ///< extended paths
	PathIdMap     id_map;        ///< path ID map
	Progress      progress;      ///< progess status
	LatchLog      log;           ///< Detailed log
	PathReadsMap  added_reads;   ///< added reads (during bridge reads extension)

	bool short_flag;             ///< short overlap extension?
	bool rjoin_flag;             ///< Read bridge join
	bool status;                 ///< extension completion status

	int extension_type;          ///< types of path extension
	int nreads;                  ///< no. of reads
	int max_libsize;             ///< maximum library size;
	int min_libsize;             ///< minimum library size;
	int lcount;                  ///< left extension success count
	int rcount;                  ///< right extension success count

	//=========
	// pointers
	//=========
	Loader*       loader;       ///< Loader object pointer
	PathEntryMap* ipaths;       ///< Assembled pths so far
	GSA*          gsa;          ///< GSA pointer
	char**        seqs;         ///< sequences
	ReadId*       pairs;        ///< paired-reads
	PathId*       preads;       ///< assembled reads so far
	


 protected:
	/** 
	 * Estimate maximum and minimum library size 
	 */
	void setLibraryRange();

	/** 
	 * Initialize path indices 
	 */
	virtual void initMaps() = 0;

	/** 
	 * Find all latchable paths and join.
	 */
	virtual void joinPaths();

	/** 
	 * Make a sequence length map
	 */
	PathLengthMap getPathLengths();

	/**  
	 * Find and latch similar paths 
	 * \param sbjct_pid current path ID
	 * \param direction latch direction
	 */
	void latchPath( PathId &sbjct_pid, int direction );

	/**
	 * Latch paths recursivley.
	 * Actual latching depends on overlap type (long or short).
	 */
	virtual bool latchOverlapPaths( PathId &spid, int direction ) = 0;

	/** 
	 * Connect two paths.
	 */
	void connect(PathId spid,
				 PathId mpid,
				 std::string &pivot,
				 std::string &match,
				 AlignSummary &summary,
				 int direction);

	/**
	 * Latch similar path to the current path.
	 * \param spid current path
	 * \param mpid similar path
	 * \param pivot current path string
	 * \param match matching path string
	 * \param summary alignment summary of two paths
	 * \param direction latch direction
	 */
	void extend( PathId &spid, 
				 PathId &mpid, 
				 std::string &pivot, 
				 std::string &match,
				 AlignSummary &summary, 
				 int direction );

	/**
	 * Generate extended sequence.
	 */
	bool updateSequence( PathId &spid, 
						 PathId &mpid, 
						 std::string &pivot, 
						 std::string &match,
						 AlignSummary &summary, 
						 int direction );

	/**
	 * Update path members
	 */
	void updateMembers(  PathId &spid, 
						 PathId &mpid, 
						 AlignSummary &summary, 
						 int direction );

	/**
	 * Update alignment positions of existing members
	 */
	void updateExistingMembers( const PathId &spid,
								const PathId &mpid,
								const int    &lgap,
								const int    &direction );
	
	/**
	 * Add new members in a merged path
	 */
	void addNewMembers( const PathId &spid,
						const PathId &mpid,
						const int    &lgap,
						const int    &direction );

	/** Update stop criteria */
	void updateStop( PathId spid,
					 PathId mpid,
					 int direction );

	/** Update trim criteria */
	void updateTrim( PathId spid,
					 PathId mpid,
					 int direction );

	/**
	 * Update index to seach further
	 */
	virtual void updateIndex( PathId spid,
							  std::string &old_pivot,
							  std::string &new_pivot,
							  int direction) = 0;
	
	/** Removed merged paths */
	void dropMergedPaths();


	/** Write out sequence after latching */
	void writeSequences();

	/** Reclaim memory */
	void clear();


	/** Align pair by comparing base by base */
	void compareBases( AlignSummary &summary, 
					   std::string &pivot, 
					   std::string &match);
	
	/** Align pair of seuqences */
	void alignPair( AlignSummary &summary,
					  std::string &pivot,
					  std::string &match,
					  int direction );

	/** Check end of each sequence to determine latchability */
	bool stopConflict( PathId spid, PathId mpid, int direction );

	/** Get maximun numerical path ID */
	PathId getMaxPathId();

	/**
	 * Find paired-end reads of between current and matching paths
	 */
	size_t findPairendReads( std::tr1::unordered_map<ReadId, bool> &pivot_reads,
							 std::string &pivot, 
							 std::string &match,
							 int direction );

	/** 
	 * Get rough sub-string from sequence based on library statistics
	 */
	std::string getStringInRange(std::string &seq, int direction, bool pivot);

	/**
	 * Find reads of given path by looking up suffix array.
	 */
	void findReads( std::string &query,
					std::tr1::unordered_map<ReadId, bool> &reads );
	
	/**
	 * Count pair-end reads of two set of reads
	 */
	size_t getPairedReads( std::tr1::unordered_map<ReadId, bool> & pivot_reads, 
						   std::tr1::unordered_map<ReadId, bool> & match_reads );


 public:
	/* Constructor */
	Connecter();

	/** Destructor */
	virtual ~Connecter();

	/** 
	 * Initalize objects.
	 * @param l Loader object (sequences, suffix arrays)
	 * @param e Assembled paths so far
	 * @param r Assembled reads so far
	 * @param t Connecter type
	 */
	void init(Loader *l, PathEntryMap *e, PathId *r, int t);


	/** Make path entries for path extension */
	void initMembership();


	//========
	// Getters
	//========
	/** \return extended paths */
	ExtenderMap*  getExtenders()  { return &extenders; }

	/** \return merted paths */
	PathIdSet     getMergePaths() { return merged_paths; } 

	/** \return path ID map */
	PathIdMap*    getPathIdMap()  { return &id_map; }

	/** Return const pointer */
	PathReadsMap* getAddedReads() { return &added_reads; }


	//==========
	// Functions
	//==========
	/** 
	 * Runner. 
	 * Single callable method to perform path extension.
	 */
	void run();

	/** Print extended object */
	void print(std::ostream &out);

	/** Dump binary objects to a file */
	void dump( std::string );

	/** Load binary objects from a file */
	void load( std::string );

	/** Remove temporary objects */
	virtual void purgeTemp() = 0;

	/** Return assembly status */
	bool getStatus() { return status; }

	/** Make path entries from extended objects */
	void makePathEntries( PathEntryMap &);
};

#endif
