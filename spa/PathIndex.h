/** 
 * \file      KmerPathIndex.h
 * \brief     Kmer Path Index.
 * \details   This class creates an inverted index between kmer and path ids.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2012
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __KMER_PATH_INDEX_H__
#define __KMER_PATH_INDEX_H__


#include <list>
#include "kmer.h"
#include "path.h"

typedef std::pair<PathId, size_t> PathFreq;
typedef std::list<PathFreq> PathFreqList;
typedef std::tr1::unordered_map<KmerId, PathFreqList> KmerPathIndexMap;
typedef std::tr1::unordered_map<KmerId, size_t> KmerFreqMap;
typedef std::tr1::unordered_map<PathId, size_t> PathFreqMap;

class KmerPathIndex 
{
 private:
	KmerPathIndexMap index;

 public:

	/** Getter: inverted index */
	KmerPathIndexMap getKmerPathIndex()
	{
		return index;
	}

	/** Reclaim memory */
	void clear() 
	{
/* 		for ( KmerPathIndexMap::iterator it = index.begin(); it != index.end(); ++it )  */
/* 			delete it->second; */
		index.clear();
	}


	KmerFreqMap getKmerFreqMap(KmerArray &kmers)
	{
		//double t0 = mytime();
		KmerFreqMap count_map;
		for ( KmerArray::iterator it = kmers.begin(); it != kmers.end(); ++it ) 
			count_map[*it] = 0;
		for ( KmerArray::iterator it = kmers.begin(); it != kmers.end(); ++it ) 
			count_map[*it]++;
/* 		for ( KmerFreqMap::iterator it = count_map.begin(); it != count_map.end(); ++it ) */
/* 			std::cerr << it->first << ":" << it->second << " "; */
/* 		std::cerr << "\n"; */
		//std::cerr << "Frequency map:" << mytime()-t0 << "\n";
		return count_map;
	}
	
	/**
	 * Initialize index of given kmer
	 */
	void init( KmerArray &kmers )
	{
		KmerFreqMap count_map = getKmerFreqMap(kmers);
		for ( KmerFreqMap::iterator it = count_map.begin(); it != count_map.end(); ++it ) {
			index[it->first] = PathFreqList();
		}
	}
 
		
	void add( KmerArray &kmers, PathId &pid )
	{
		//double t0 = mytime();
		KmerFreqMap count_map = getKmerFreqMap(kmers);
		for ( KmerFreqMap::iterator it = count_map.begin(); it != count_map.end(); ++it ) {
			index[it->first].push_back( std::pair<PathId, size_t>( pid, it->second ) );
		}
		//std::cerr << "Add entry:" << mytime()-t0 << "\n";
	}


	/** Does exist the key ? */
	bool has( KmerId &key ) 
	{
		if ( index.find(key) == index.end() ) return 0;
		return 1;
	}

	/** Erase the key/value from index */
	void erase( KmerId &key )
	{
		if ( ! has(key) ) {
			std::cerr << "Key does not exist\n";  
			return;
		}
		//delete index[key];
		index.erase(key);
	}

	void remove( KmerId &key, PathId &path )
	{
		if ( ! has(key) ) {
			std::cerr << "Key does not exist\n";  
			return;
		}
		for ( PathFreqList::iterator it = index[key].begin(); it != index[key].end(); ) {
			if ( it->first == path ) {
				index[key].erase(it++);
				break;
			} else ++it;
		}
	}

	void remove( KmerArray &kmers, PathId &path )
	{
		for ( KmerArray::iterator kt = kmers.begin(); kt != kmers.end(); ++kt ) 
			remove( *kt, path );
	}

	
	/** Return values by key */
	PathFreqList getValue( KmerId &key ) 
	{
		if ( ! has(key) ) 
			return PathFreqList();
		return index[key];
	}

	/** Index size */
	int getSize() 
	{
		return index.size();
	}


	//PathFreqMap 
	void setPathFreqMap( PathFreqMap &path_freq, KmerFreqMap &count_map )
	{
		//double t0 = mytime();
		//PathFreqMap path_freq;

		for ( KmerFreqMap::iterator it = count_map.begin(); it != count_map.end(); ++it ) {
			KmerId kmer = it->first;
			size_t freq = it->second;
			
			KmerPathIndexMap::iterator jt = index.find(kmer);
			if ( jt == index.end() ) continue;
			
			for ( PathFreqList::iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt ) {
				PathId pid = kt->first;
				size_t num = kt->second;
				
				size_t min = freq > num ? num : freq;

				if ( path_freq.find(pid) == path_freq.end() ) 
					path_freq[pid] = 0;
				path_freq[pid] += min;
			}
		}

		//std::cerr << "Path Frequency map:" << mytime()-t0 << "\n";
		//return path_freq;
	}
	
	void findPaths( PathFreqMap &path_freqs, KmerArray &kmers ) 
	{		
		KmerFreqMap count_map = getKmerFreqMap(kmers);
		setPathFreqMap(path_freqs, count_map);
	}
};

#endif
