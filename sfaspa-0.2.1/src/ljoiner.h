/** 
 * \file      ljoiner.h
 * \brief     Long overlapping path extension
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \bug       None.
 * \warning   None.
 * \date      Fri 2014-01-043
 * \copyright J. Craig Venter Institute.
 */

#ifndef __LONG_JOINER_H__
#define __LONG_JOINER_H__

#include "connecter.h"

/**
 * \class LongJoiner
 * \brief Path joiner for long overlapping paths.
 */
class LongJoiner : public Connecter
{
 private:
	KmerToPathMap    path_map; ///< kmer to path id mapping
	PathKmerIndex    path_kmers;
	PathKmerPosIndex path_kposs;

 private:
	/**
	 * Initialize kmer to path ID mapping
	 */
	void initMaps();

	/**
	 * Find overlapping paths and join them.
	 */
	bool latchOverlapPaths( PathId &spid, 
							int direction );

	/** 
	 * Find similar paths by k-mer counting
	 */
	void findPathsByKmerCount( std::vector<JoinEntry> &entries, 
							   PathId spid, 
							   int direction );
	/** 
	 * Join a path if possible.
	 */
	bool tryConnectOverlap( JoinEntry &entry, 
							PathId &spid, 
							std::string &pivot,
							const KmerArray *pkmers,
							const KmerPosMap *pposs, 
							std::tr1::unordered_map<ReadId, bool> &pivot_reads,
							int direction );

	bool checkByBases( AlignSummary &summary,
					   std::string &sub_pivot,
					   std::string &sub_match );

	bool checkByAlign( AlignSummary &summary,
					   std::string &sub_pivot,
					   std::string &sub_match,
					   int direction );

	/**
	 * Determine whether a path can be extendible.
	 */
	bool latchableCandidate( int &qs, 
							 int &ss, 
							 int &qe, 
							 int &se,
							 PathId spid, 
							 PathId mpid, 
							 std::string &pivot, 
							 std::string &match, 
							 const KmerArray *pkmers,
							 const KmerArray *mkmers,
							 const KmerPosMap *pposs,
							 const KmerPosMap *mposs,
							 int direction );

	/** 
	 * Make shared kmer counts of paths to current path
	 */
	void getKmerCountMap( std::map<PathId, int> &count_map,
						  PathId spid,
						  int direction );
	
	/**
	 * Update path-map of extended path & merged path
	 */
	void updateIndex( PathId spid,
					  std::string &old_pivot,
					  std::string &new_pivot,
					  int direction);

	/** 
	 * Add new entries of kmer from given path.
	 */
	void addPathIds( //KmerToPathMap &path_map,
					 PathId pid,
					 std::string sequence );

	/** 
	 * Drop kmer entries of path path.
	 */
	void dropPathIds( //KmerToPathMap &path_map,
					  PathId pid,
					  std::string sequence );


 public:
	LongJoiner();
	virtual ~LongJoiner();

	/** Release temporary storages */
	void purgeTemp();
};

#endif
