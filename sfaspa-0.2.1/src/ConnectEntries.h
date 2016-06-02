#include "core.h"

/**
 * \class JoinEntry
 * \brief Simplified version of path entry during path extension
 */
struct JoinEntry
{
	PathId pid;
	int nkmers;
	int seqlen;
	JoinEntry(PathId p, int n, int s );
};

/**
 * @brief   Read entry.
 * @details Read entry of rbriding read. It assumes 100% aligned.
 */
struct ReadEntry
{
	ReadId   read;
	unsigned rpos;
	unsigned ppos;
	ReadEntry();
	ReadEntry( ReadId i, unsigned r, unsigned p ) ;

	void dump( std::fstream &out );
	void load( std::fstream &in  );
};

typedef std::list<ReadEntry> ReadEntryList;
typedef std::tr1::unordered_map<PathId, ReadEntryList> PathReadsMap;
