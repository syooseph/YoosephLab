/** 
 * @file       metapath.h
 * @brief      Meta-path (augmented path) is a base class of Place and Recruiter
 * @date       Modified on Tue 2013-12-17 06:38:51 PM
 * @author     Youngik Yang
 * @version    0.02
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __META_PATH_H__
#define __META_PATH_H__

#include "merger.h"
#include "connecter.h"
#include "profile.h"
#include "anchor.h"
#include "progress.h"
//#include "summary.h"

typedef std::vector<bool> FlagArray;
typedef std::tr1::unordered_map<ReadId, std::vector<int> > ReadPosMatrixMap;
typedef std::vector<PathAligner> PathAligners;
typedef std::tr1::unordered_map<size_t, PathId> NumToPathMap;
typedef std::tr1::unordered_map<PathId, size_t> PathToNumMap;
typedef std::tr1::unordered_map<ReadId, size_t> ReadFreqMap;

const std::string PlacerFiles[] = { "placer.align.aac", "placer.profile.aac", "placer.fasta", "placer.place.txt" };
const std::string RecruiterFiles[] = { "recruiter.align.aac", "recruiter.profile.aac", "recruiter.fasta", "recruiter.place.txt" };

enum { DERIVED_PLACE, DERIVED_RECRUIT };

/**
 * \brief Read placement, MSA, and profile
 */
class MetaPath
{
 protected:
	int           nreads;
	size_t        npaths;
    PathAligners  Alns;
	NumToPathMap  Id2Paths;
	PathToNumMap  Path2Ids;
	FlagArray     bad_paths;
	Progress      progress;
	int           derived;
	bool          created;
	bool          verbose;
	bool          status;

	/* Reference loader objects */
	char    **seqs;


	PathId  *preads;

	/* Other assembly objects */
	ReadAlignerMap  *read_aligns;
	Clusters        *clusters;
	Clusters        *meta_ccs;
	//Summary         *summary;
	PathEntryMap    *path_entries;
	ExtenderMap     *read_exts;
	ExtenderMap     *tiny_exts;
	ExtenderMap     *long_exts;
	PathIdMap       *read_idmap;
	PathIdMap       *tiny_idmap;
	PathIdMap       *long_idmap;
	PathIdMap       *cluster_idmap;
	PathReadsMap    *added_reads;

 private:
	void writePlacement();
	void writeAlignment();
	void writeSequences();
	void combine( std::string name );
	
 public:
	MetaPath();
	virtual ~MetaPath();

	/**
	 * Copy constructor 
	 */
	MetaPath(const MetaPath &source);

	/**
	 * Operator overloading
	 */
	MetaPath& operator= (const MetaPath &source);

	void init( Loader *l, Extracter *p, Merger *m, Connecter *e, Connecter *s, Connecter *b, PathEntryMap *y );
	void init( Loader *l, Extracter *p, Merger *m, Merger *c, Connecter *e, Connecter *s, Connecter *b, PathEntryMap *y );
	virtual void run();
	void reset();
	void clear();

	void dump( std::string );
	void load( std::string );

	void write();
	void copy( const MetaPath& );

	void updatePathMembership();
	int countAssembledReads();

	PathLengthMap getPathLengths();
	PathAligners *getAlignments() { return &Alns; }
	NumToPathMap *getPathIds()    { return &Id2Paths; }
	PathToNumMap *getPathNums()   { return &Path2Ids; }
	FlagArray    *getBadPaths()   { return &bad_paths; }


	bool getStatus() { return status; }
};

#endif
