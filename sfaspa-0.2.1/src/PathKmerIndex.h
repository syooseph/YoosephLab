// Mon 2014-01-27 05:29:42 PM

#ifndef __PATH_KMER_INDEX_H__
#define __PATH_KMER_INDEX_H__

#include "core.h"

//#include <unordered_map>

typedef std::tr1::unordered_map<PathId, KmerArray> PathKmersMap;

class PathKmerIndex
{
 private:
	PathKmersMap index;

 public:

	/** Get entire index */
	PathKmersMap getMap() const;
	
	/** Add one entry */
	void add( const PathId pid, const KmerArray &kmers );

	/** Erase the key/value from index */
	void erase( const PathId pid );

	/** Update the key with values */
	void update( const PathId pid, const KmerArray &kmers );

	/** Index size */
	int getSize() const;

	/** Reclaim memory */
	void clear();

	/** Does exist the pid ? */
	bool has( const PathId pid ) const;

	/** Return values by key */
	KmerArray getValue( const PathId pid ) const;

	/** Return pointer by key */
	const KmerArray* getPointer( const PathId pid );

};


#endif
