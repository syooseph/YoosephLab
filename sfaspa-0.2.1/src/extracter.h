/** 
 * \file extracter.h
 * \brief Path extraction
 * \details This header defines path extraction class
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \date      Modified on Wed 2013-10-23 05:50:39 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __EXTRACTER_H__
#define __EXTRACTER_H__

#include <thread>
#include <cmath>
#include <type_traits>
#include "loader.h"
#include "ralign.h"
#include "log.h"
#include "seed.h"
#include "progress.h"
#include "PathEntry.h"
#include "ReadStartCount.h"
//#include "ReadEndCount.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

typedef std::list<PathId> PathIdList; 
typedef std::tr1::unordered_map<Vertex, PathIdList> VertexToPathMap;
typedef std::tr1::unordered_map<PathId, ReadAligner> ReadAlignerMap;
//typedef std::tr1::unordered_map<Vertex, bool> NodeFlags;
/** Word (string) to Paths map */
typedef std::tr1::unordered_map<std::string, std::set<PathId> > WordToPathMap; 
typedef std::tr1::unordered_map<PathId, int> PathIdPosMap;
typedef std::set<PathId> PathIdSet;
/** Trace strings during path extraction */
typedef std::tr1::unordered_map<std::string, bool> TraceStringMap;
typedef std::tr1::unordered_map<KmerId, bool> KmerFlags;

/* struct BoundLength */
/* { */
/* 	BoundType bound; */
/* 	size_t length; */
/* }; */

/* typedef std::unordered_map< */
/* BoundLength,  */
/* 	int, */
/* 	std::hash< */
	

/* typedef std::multimap<BoundType, UsedInfo> ReadStartMap; */
/* typedef std::multimap<BoundType, UsedInfo> ReadEndMap; */
//typedef std::multimap<BoundType, UsedInfo> BoundUsedMap;
//typedef std::map<UsedBound, size_t> UsedBoundMap;		


//typedef std::unordered_map< std::string, BoundArray > SearchHistory;

typedef std::pair<unsigned, BoundType> GsaBound;
//typedef std::set<GsaBound> GsaBoundSet;
typedef std::unordered_multimap< std::string, GsaBound> SearchHistory;

/**
 * \class Extracter
 * \brief Path extraction
 */
class Extracter
{
 private:
	//WordFreqMap     nback_freqs;   ///< back_trace+1 mer frequency map
	AdjacencyMap    ancestors;     ///< predecessor node map
	ReadAlignerMap  path_map;      ///< PathId and ReadAlignment mapping
	VertexToPathMap discovered;    ///< discovered nodes
	DeBruijnGraph   graph;         ///< k-mer graph
	CoverageMap     kmer_coverage; ///< k-mer depths
	VertexToKmerMap vertex_map;    ///< mapping between nodes and k-mers
	SearchLog       srch_log;      ///< assembly summary 
	ReadAlignLog    read_log;      ///< Read alignment log
	WordToPathMap   word_map;      ///< (k+2)-mer to PathId mapping 
	PathIdSet       check_pid;     ///< similar paths for (k+2)-mer to nback-mer match 
	PathIdPosMap    check_pos;     ///< check position of the check-pids
	TraceStringMap  traces_map;    ///< tracking of nback+1 mers for cycle detection
	ReadStartCount  *read_starts;
	//ReadEndCount    *read_ends;
	/* UsedBoundMap    check_starts;  ///< Read start count for one path */
	/* UsedBoundMap    check_ends;    ///< Read end count for one path */
	UsedBoundMap    used_starts;   ///< Read start count for all path
	UsedBoundMap    used_ends;     ///< Read end count for all path
	SearchHistory   rhist_map;    ///< Right extension history during each path extraction
	//SearchHistory   lhist_map;    ///< Left extension history during each path extraction

	NodeSet deleted_nodes;	
	
	size_t npaths;                 ///< no. of paths
	PathId *preads;                ///< read membership
	NodeFlags used_seeds;         ///< path nodes of failed path search
	NodeFlags srch_failed;         ///< path nodes of failed path search
	NodeFlags read_failed;         ///< path nodes of failed read alignment
	NodeFlags save_failed;         ///< path nodes of failed read alignment
	NodeFlags trimmed_nodes;
	KmerFlags trimmed_kmers;

Seed seed_nodes;               ///< seed object

	Progress progress;          

	/* Reference object for convenience */
	Loader   *loader;      ///< file loader object
	GSA      *gsa;         ///< generalized suffix arrays
	char     **seqs;       ///< sequence reads
	int      nreads;       ///< no. of sequences
	IntArray *bad_reads;   ///< invalid read IDs

	bool     status;

 private:
	/** 
	 * Constrcut graph 
	 */
	void buildGraph();

	/**
	 * Find minimum seed coverage when x% of k-mers are used for seeds.
	 * Skip when minimum seed coverage is given.
	 */
	void setStartMinCoverage();

	/**
	 * Set seed k-mer priority by successivly calling updatePathSearchPriority
	 */
	void initPathSearchPriority();

	void makeFailedSeedsPriority();

	/**
	 * Calculate priroity score for given vertex (k-mer).
	 * Score = coverage / exp(degrees)
	 */
	void updatePathSearchPriority(Vertex);

	/** 
	 * Extract all paths
	 */
	void extractPaths(double);

	/* void convert(); */
	void release();

	/** 
	 * Extract a single path starting from given seed kmer.
	 * Actual path search is done by calling extractGreedyBestPath.
	 */
	bool extractOnePath( Vertex seed, GraphPath &gpath );

	/**
	 * Determine whether a given seed is valid for path retrieval or not
	 */
	bool validStartNode( Vertex &seed );

	/**
	 * Return true if a path is short.
	 * False otherwise.
	 */
	bool shortPath( PathType &path );
	
	/**
	 * Extract a path from a given seed and save the path as a GraphPath object.
	 * A node is prepended/appended by calling extendNode
	 */
	bool extractGreedyBestPath( Vertex seed, GraphPath &graph_path );

	
	/**
	 * A single node exension.
	 * It checks whether current path found already. If so, extension stops.
	 * Otherwise, as long as it is not a terminal node, it searches the best neighboring node to extend judging by same reads in sub-path.
	 * If the neghiboring node is weakly supported, extension halts.
	 * Otherwise, it is prepended/appended to the path and repeat the steps.
	 */
	void extendNode( Vertex curr,
					 std::string &path_str,
					 //std::string &nback_str,
					 GraphPath &gpath,
					 //CyclicPath &cpath,
					 int direction,
					 int &reason );
	
	/**
	 * Terminal node check
	 * Stop codon, sink or source.
	 */
	bool isTerminalNode( Vertex curr,
						 GraphPath &gpath,
						 int direction,
						 int &reason );

	/**
	 * Find the best extension node.
	 * The max node is determined by the number of same reads across sub-path and the neighboring node.
	 * The max node search is done by calling getNeighboringNodes followed by findGreedyBestNeighbor.
	 */
	bool findMaxNeighbor( Vertex curr,
						  std::string &path_str,
						  //std::string &nback_str,
						  GraphPath &gpath,
						  int direction,
						  int &reason,
						  Vertex &max_node,
						  BoundArray &max_bound,
						  std::string &max_trace,
						  UsedBoundInfoList &max_useds );


	/**
	 * Check whether current path traps into same region.
	 * Cycle detection using suffix array.
	 */
	bool fromSameNode( std::string &max_trace,
					   PathType &path,
					   int direction );

	/**
	 * Get a sub-path of given length
	 */
	PathType getSubPath( PathType &path, int direction, int length );

	/**
	 * Get all neighboring nodes 
	 */
	NodeArray getNeighboringNodes( Vertex &curr,
								 //GraphPath &gpath,
								 int direction );

	/**
	 * Order nodes by read supports.
	 * Ordering is used for the best neighboring node search.
	 */
	void orderNodes( NodeArray &vset,
					 NodeArray &nodes,
					 std::vector<size_t> &covers );
	
	/** 
	 * Sanity check of neighboring nodes
	 */
	NodeArray findGoodNeighbors( Vertex, NodeArray&,int );

	/**
	 * Find the maximal node to extend using multiple CPUs
	 */
	void findGreedyBestNeighbor( GraphPath &,
								 NodeArray &, 
								 Vertex, 
								 std::string &path_str,
								 //std::string &nback_str,
								 //PathType &, 
								 Vertex &, 
								 int&, 
								 BoundArray &, 
								 std::string &,
								 UsedBoundInfoList &max_useds,
								 int);

	/**
	 * Find the maximal node to extend using single CPU.
	 */
	void findGreedyBestNeighborSP( NodeArray &, 
								 Vertex, 
								 std::string &nback_str,
								   //PathType &, 
								 Vertex &, 
								 int&, 
								 BoundArray &,
								   std::string &, 
								 int);

	//void combineMap( UsedBoundMap &master, UsedBoundMap &slave );
	void combineHist( SearchHistory &n, int direction );

	int getWeight( Vertex &curr_node, 
				   Vertex &v,
				   std::string &path_str,
				   std::string &trace_str,
				   char ch,
				   BoundArray &prev_bound,
				   BoundArray &curr_bound,
				   int direction );

	int getSubstringMatches( std::string &path_str,
							 std::string &trace_str,
							 char ch,
							 BoundArray &prev_bounds,
							 BoundArray &bound,
							 PathIdSet &Check_Pids,
							 PathIdPosMap &Check_Poss,
							 int direction );

	int getSizeOfHistSearch( char ch,
							 std::string &path_str, 
							 size_t max_size, 
							 int direction );

	int getStartReadsOnPath( std::string &path_str,
							 std::string &trace_str,
							 char ch,
							 std::vector<BoundType> &bound,
							 UsedBoundInfoList &Check_Starts,
							 SearchHistory &Check_Hists,
							 int direction );

	int getEndReadsOnPath( std::string &path_str,
						   std::string &trace_str,
						   char ch,
						   std::vector<BoundType> &bound,
						   UsedBoundInfoList &Check_Ends,
						   //SearchHistory &Check_Hists,
						   int direction );

	size_t getUsedCounts( UsedBoundMap *used_map,
						  BoundType &srch,
						  size_t srch_len,
						  size_t gsa_id,
						  bool start );

	BoundType refineRightBoundary(int p, 
								  BoundType &srch, 
								  std::string &srch_str,
								  size_t &srch_len,
								  int &Num, 
								  int &Used,
								  BoundArray &Bounds,
								  std::vector<int> &Ends);

	size_t getReadEndFromGsa( size_t left,
							  size_t right,
							  size_t length,
							  size_t sfa_no,
							  size_t rid_no );
	
	/** 
	 * Count recruited reads in other paths of given string size.
	 * String size must be between k+2 mer and nback-mer.
	 */
	int getIntervalCounts(std::string &trace_str, 
						  PathIdPosMap &check_pos,
						  PathIdSet &check_pid,
						  int direction);

	void findTraceStringFromPaths( std::string &trace_str, 
								   PathIdPosMap &check_pos,
								   PathIdSet &check_pid );

	/**
	 * Check validity of a path
	 */
	bool checkPath( GraphPath &gpath, ReadAligner &raln, ReadAlignLog &rlog );

	/**
	 * Valid path length?
	 */
	bool validLengthPath( ReadAligner &raln );

	/** 
	 * Trim nodes from the path after read placement
	 */
	void trimPath( GraphPath &p, ReadAligner &r );

	bool update( GraphPath &p, ReadAligner &r );

	/** 
	 * Save a discovered path and nodeos
	 */
	void savePath( GraphPath &p, ReadAligner &r, PathId pid );

	/** 
	 * Update kmer coverage
	 */
	bool updateCoverage( GraphPath &gp, ReadAligner &raln );

	/** 
	 * Update edge weight
	 */
	void updateGraph( GraphPath &gp, ReadAligner &raln );

	void dropWeakNodes( GraphPath &gpath );
	void dropIslands( GraphPath &gpath );

	/** 
	 * Update nback history
	 */
	/* void updateNback( ReadAligner &raln ); */

	void updateUsedBounds( ReadAligner &, bool ); 
	/* void updateReadStarts( GraphPath &p ); */

	/* void updateReadEnds(); */

	/**
	 * Update read membership
	 */
	void updateReads( ReadAligner &raln, PathId pid );

	/** 
	 * Update seed after a path retrieval
	 */
	void updateSeed( GraphPath &p );


	/** 
	 * Remove sub paths becuase of a cycle.
	 */
	void removeSubPath( GraphPath &gpath,
						int direction );
	

	/**
	 * Show path extraction summary
	 */
	void displaySummary();	

	/** 
	 * Write path sequences to a file
	 */
	void writePath();

	/**
	 * Erase correspoding entries from a map when graph is trimmed.
	 */
	void trimMap( );

	/**
	 * Record k-mers of the corresponding nodes in a path
	 */
	void setPathKmers( GraphPath &p );

	/**
	 * Release graph related objects
	 */
	void purgeGraph();

	/** 
	 * Assigned reason of failure to nodes from a path.
	 */
	void addFlags( PathType &path, NodeFlags &flags );

 public:
	/** 
	 * Constructor
	 */
	Extracter();

	/** 
	 * Constructor with loader object 
	 */
	Extracter(Loader *l);

	/**
	 * Destructor
	 */
	~Extracter();
	
	/**
	 * Initialize
	 */
	void init(Loader *l);

	/**
	 * Top level function for path extraction
	 */ 
	void run();

	/** 
	 * Print read length statistics.
	 */
	void stat();

	/** 
	 * Trim low supported nodes
	 */
	void trim();

	/** 
	 * Extrac paths
	 */
	void extract();

	/** 
	 * Clear resources
	 */
	void clear();

	//========
	// Getters
	//========
	SearchLog       getSearchLog() { return srch_log; } 
	ReadAlignerMap* getReadAlignerMap(){ return &path_map; }
	PathId*         getReadMembership() { return preads; }


	/**
	 * Dump binary graph path objects
	 */
	void dump( std::string );

	/**
	 * Load binary graph path objects
	 */
	void load( std::string );

	/** 
	 * Count number of assembled reads 
	 */
	int countAssembledReads();

	bool getStatus() { return status; }

	/* void setAssembledReads( PathIdArray & ); */
	/* void setPathEntries( PathEntryMap &); */
	
	void makePathEntries( PathEntryMap &);
};

#endif
