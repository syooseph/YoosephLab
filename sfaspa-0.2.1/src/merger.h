/** 
 * @file       merger.h
 * @brief      Path clustering
 * @date       Modified on Tue 2013-12-17 06:49:14 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __MERGER_H__
#define __MERGER_H__

#include "extracter.h"
#include "KmerPathIndex.h"
#include "PathKmerIndex.h"
#include "PathEntry.h"
#include "extender.h"
#include "msectioner.h"
#include "cluster.h"
//#include "prefilter.h"
#include "log.h"

typedef std::multimap<size_t, PathId> PathLengthMap;
//typedef PathBins Clusters;
//typedef std::tr1::unordered_map<PathId, PathEntry> PathEntryMap;

//typedef std::list<AlignSummary> AlignSummaryList;


typedef std::tr1::unordered_map<PathId, Cluster> Clusters;

/**
 * \brief Clustering paths
 */
class Merger
{
 private:
	PathEntryMap     *path_entries;
	PathIdMap        id_map;        ///< path ID map
	KmerPathIndex    kmer_paths;
	PathKmerIndex    path_kmers;
	PathKmerPosIndex path_kposs;
	Clusters         clusters;
	PathIdSet        merged_paths;
	Progress         progress;
	/* Prefilter       prefilter; */
	MergeLog         log;

	int kmer_size;

	/* // stat */
	/* size_t count_pat; // # sum of total path candidates by k-mer similarity. */
	/* size_t count_aln; */
	/* size_t count_suc; */
	/* size_t count_lgap; */
	/* size_t count_egap; */
	/* double align_time; */

	bool recall;
	bool status; /// run? loaded?
	
 private:
	void build();
	void initCluster();
	void cluster();
	PathLengthMap getPathLengths();
PathLengthMap getPivotLengths();

	bool findCluster(const KmerArray *query_kmers,PathId &query_pid, PathId &sbjct_pid, AlignSummary &summary);

	int getNewFilterKmer( size_t str_len );
	void findSimilarPaths( //PathFreqMap &path_freq,
						  int min_filter,
						  std::multimap<size_t, PathId> &freq_map,
						  PathId &query_pid,
						  std::string &query,
						   const KmerArray *query_kmers );

	bool mergible( int min_filter,
				   std::string &sbjct,
				   std::string &query, 
				   const KmerArray *skmers,
				   const KmerArray *qkmers,
				   const KmerPosMap *sposs,
				   const KmerPosMap *qposs,
				   AlignSummary &summary );

	bool doBaseAlign( AlignSummary &summary,
					  std::string &query,
					  std::string &sbjct,
					  Section &s );

	bool doPairAlign( AlignSummary &summary,
					  std::string &query,
					  std::string &sbjct,
					  Section &s );;
		
	void dropMergedPaths();
	void update();
	//void addSingleton();
	void writeClusters();
	//void makeSingletons();

	void rebuildIndex();
	
	std::string getPathSequence(PathId &pid);

	PathId getMaxPathId();

 public:
	Merger();
	~Merger();
	/* void init( Extracter *e ); */
	/* void init( ExtenderMap *e ); */
	/* void init( PathEntryMap &p ); */
	//void init( const PathEntryMap &p );
	void init( PathEntryMap *p );
	void release();
	Clusters *getClusters() { return &clusters; }
	//PathEntryMap *getPathEntries() { return &path_entries; }
	void run();
	void print();
	void clear();
	void setMergedPaths(PathIdSet &s) { merged_paths = s; }
	PathIdSet getMergedPaths() { return merged_paths; }
	void dump( std::string );
	void load( std::string );
	void makeClusters( Extracter *e );

	void setRecall( bool r ) { recall = r; }
	bool getStatus() { return status; }

	
	void makePathEntries( PathEntryMap  &);

	/** \return path ID map */
	PathIdMap*    getPathIdMap()  { return &id_map; }

};

#endif
