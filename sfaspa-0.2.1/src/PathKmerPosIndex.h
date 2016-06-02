// Mon 2014-01-27 05:29:42 PM

#ifndef __PATH_KMER_POS_INDEX_H__
#define __PATH_KMER_POS_INDEX_H__

#include "core.h"

//#include <unordered_map>

typedef std::tr1::unordered_map<KmerId, std::vector<size_t> > KmerPosMap;
typedef std::tr1::unordered_map<PathId, KmerPosMap> PathKmerPosMap;

class PathKmerPosIndex
{
 private:
	PathKmerPosMap index;

 public:

	/** Get entire index */
	PathKmerPosMap getMap() const;
	
	/** Add one entry */
	void add( const PathId pid, const KmerPosMap &kmers );

	/** Erase the key/value from index */
	void erase( const PathId pid );

	/** Update the key with values */
	void update( const PathId pid, const KmerPosMap &kmers );

	/** Index size */
	int getSize() const;

	/** Reclaim memory */
	void clear();

	/** Does exist the pid ? */
	bool has( const PathId pid ) const;

	/* /\** Return values by key *\/ */
	KmerPosMap getValue( const PathId pid ) const;

	/** Return values by key */
	const KmerPosMap* getPointer( const PathId pid ) const;

};


#endif
