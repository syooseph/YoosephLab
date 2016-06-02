/** 
 * @file       recruiter.h
 * @brief      Read recruiter
 * @details    Previously unassigned reads are recruited.
 * @date       Thu 2014-01-02
 * @author     Youngik Yang
 * @version    0.02
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __RECRUITER_H__
#define __RECRUITER_H__

#include "metapath.h"
//#include "PathKmerIndex.h"
#include "PathKmerPosIndex.h"
#include "log.h"

/**
 * \brief Read recruiter
 */
class Recruiter : public MetaPath
{
 private:
	KmerPathIndex    kmer_paths; ///< kmer to path index
	PathKmerPosIndex path_kposs; ///< path to kmers

	RecruitLog recruit_log;
	AnchorLog  anchor_log;

 private:
	void clear();
	void recruit();
	/* void __recruitSP( ReadIdArray &); */
	/* void __recruitMP( ReadIdArray &); */
	void __recruit( ReadIdArray &);
	void makeIndex();

	PathId findPath( ReadId i, ReadPlacement &place );
	PathId findBestPathByAlignment( ReadId rid, std::string &read,
									std::multimap<size_t, PathId> &freq_map,
									ReadPlacement &place ,
									KmerArray &query_kmers,
									int min_filter );
	PathId findBestPathByAnchoring( ReadId rid, std::string &read,
									std::multimap<size_t, PathId> &freq_map,
									ReadPlacement &place,
									KmerArray &query_kmers,
									int min_filter);
	IntPair findSbjctRegion(std::string &read, 
							std::string &sbjct,
							KmerArray &query_kmers,
							KmerArray &sbjct_kmers);

	bool findRegion( std::multimap<int,int> &poss_pair,
					 std::multimap<int,int>::iterator &pt,
					 IntPair &rb,
					 IntPair &sb,
					 int &ct,
					 KmerArray &query_kmers );
	
	void makePosPairMap( std::multimap<int,int> &poss_pair,
						 KmerArray &query_kmers,
						 KmerArray &sbjct_kmers );
	void extendRegion( IntPair &rb,
					   IntPair &sb,
					   int qnkmers,
					   int snkmers );
	bool passKmerFilter( IntPair &rb,
						 IntPair &sb,
						 KmerArray &query_kmers,
						 KmerArray &sbjct_kmers );
	
	void findCandidatePaths( std::string &read,
							 KmerArray &query_kmers,
							 std::multimap<size_t, PathId> &freq_map,
							 int min_filter
							 );

	void printSummary( double t0 );

 public:
	Recruiter();
	~Recruiter();
	void run();
};

#endif
