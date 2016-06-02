/** 
 * \file      KmerPathIndex.h
 * \brief     Kmer to path indexing
 * \details   This class creates an inverted index between kmer and path ids.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2012
 * \date      modified: Fri 2013-12-13 06:24:58 PM
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __KMER_PATH_INDEX_H__
#define __KMER_PATH_INDEX_H__

#include "core.h"
//#include <set>
//#include <unordered_map>
//#include <tr1/unordered_map>

/* #include <list> */
/* #include <boost/unordered_map.hpp> */
/* #include "kmer.h" */
/* #include "path.h" */

typedef std::pair<PathId, size_t> PathFreq;
typedef std::list<PathFreq> PathFreqList;
//typedef std::vector<PathFreq> PathFreqArray;
typedef std::tr1::unordered_map<KmerId, PathFreqList> KmerPathIndexMap;
//typedef std::tr1::unordered_map<KmerId, PathFreqArray> KmerPathIndexMap;
//typedef std::tr1::unordered_map<KmerId, std::set<PathId> > KmerPathIndexMap;
typedef std::tr1::unordered_map<KmerId, size_t> KmerFreqMap;
typedef std::tr1::unordered_map<PathId, size_t> PathFreqMap;

/**
 * \brief Path and its kmer index map
 */
class KmerPathIndex 
{
 private:
	KmerPathIndexMap index;

 public:

	/** Getter: inverted index */
	KmerPathIndexMap getKmerPathIndex();

	/** Reclaim memory */
	void clear();

	/** Kmer frequences */
	KmerFreqMap getKmerFreqMap(const KmerArray &kmers);

	void setKmerFreqMap( KmerFreqMap &count_map, KmerArray &kmers);

	/**
	 * Update indices of given kmers.
	 * key: kmer
	 * vlaue: path id and kmer frequencies
	 */
	void add( const KmerArray &kmers, PathId &pid );

	/** Does exist the key ? */
	bool has( KmerId &key ) ;

	/** Erase the key/value from index */
	void erase( KmerId &key );

	/** Return values by key */
	//PathFreqList getValue( KmerId &key ) ;
	
	/** Index size */
	int getSize() ;

	/**
	 * Count sum of all kmer occurrences of each path found in kmer count map
	 */
	void setPathFreqMap( PathFreqMap &path_freq, KmerFreqMap &count_map );

	/**
	 * Found of paths of given kmers
	 */
	void findPaths( PathFreqMap &path_freqs, KmerArray &kmers ) ;
	void findPaths( PathFreqMap &path_freqs, const KmerArray *kmers ) ;
};

#endif
