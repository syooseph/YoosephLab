#ifndef __SFA_READ_JOINER_H__
#define __SFA_READ_JOINER_H__

#include "rjoiner.h"


//typedef std::vector<char*> EndStrings;

struct PosInfo
{
	int      rpos;
	int      ppos;
	int      mlen;
	/* PosInfo() */
	/* { */
	/* 	rpos = ppos = mlen = 0; */
	/* } */
	PosInfo(int p1, int p2, int l)
	{
		rpos = p1; ppos = p2; mlen = l;
	}
};

typedef std::tr1::unordered_map< PathId, PosInfo > PathPosInfoMap;
typedef std::tr1::unordered_map< ReadId, PosInfo > ReadPosInfoMap;

struct ReadMatchPos
{
	ReadId   read;
	int      rpos;
	int      ppos;
	int      mlen;
	ReadMatchPos(ReadId r, int p1, int p2, int l)
	{
		read = r; rpos = p1; ppos = p2; mlen = l;
	}
	void print()
	{
		printf("(i:%u,r:%d,p:%d,l:%d)", read, rpos, ppos, mlen);
	}
};

/* struct PathMatchPos */
/* { */
/* 	PathId   path; */
/* 	int      rpos; */
/* 	int      ppos; */
/* 	int      mlen; */
/* 	PathMatchPos(PathId &p, int p1, int p2, int l) */
/* 	{ */
/* 		path = p; rpos = p1; ppos = p2; mlen = l; */
/* 	} */
/* }; */

/* typedef std::list<PathMatchPos> PathMatchPosList; */
/* //typedef std::vector<PathMatchPosList> PathMatchPosL */
typedef std::list<ReadMatchPos> ReadMatchPosList;
typedef std::tr1::unordered_map<PathId, ReadMatchPosList> ReadMatchPossMap;	

typedef std::tr1::unordered_map< PathId, IntArray> PathPossMap;

//typedef std::tr1::unordered_map<ReadId, bool> ReadFlagMap;

//typedef std::vector<GsaType> GsaTypeArray;

/**
 * \class ReadJoinerSA
 * \brief Path joiner for bridging read support
 */
class ReadJoinerSA : public ReadJoiner
{
 private:
	/* EndStrings lend_strs; */
	/* EndStrings rend_strs; */
	char **lend_path_strs;
	char **rend_path_strs;
	GSA *lend_gsa; ///< GSA of left end str of paths
	GSA *rend_gsa; ///< GSA of right end str of paths

	size_t *lend_path_lens;
	size_t *rend_path_lens;

	PathIdMap gsa_ids;
	//PathIdMap rgsa_ids;

	ReadMatchPossMap lend_path_reads;
	ReadMatchPossMap rend_path_reads;

	ReadToPathsMap lend_read_paths;
	ReadToPathsMap rend_read_paths;
	int max_read_len;


	double t_suffix, t_trim, t_align, t_save;
	double t_trim_right, t_trim_left, t_trim_erase;
	double t_trim_insert;
	double t_align_check, t_align_insert, t_align_rlen, t_align_getsa, t_align_plen;
	size_t n_align_total, n_align_check, n_align_valid;
 private:

	// pure virutal function from Connectr class
	void initMaps();

	void getMaxReadLength();

	// pure virutal function from Connectr class
	bool latchOverlapPaths( PathId &spid, int direction );

	// pure virutal function from Connectr class
	void updateIndex( PathId spid,
					  std::string &old_pivot,
					  std::string &new_pivot,
					  int direction);
	

	/* // pure virutal function from ReadJoiner */
	/* void stitchPaths(); */

	// function overloading (Base: Connecter)
	void joinPaths();

	void extractEndStrings();
	void extractOneEndStrings(int direction);
	
	void makeSuffixArrays();
	void makeOneSuffixArray(int direction);

	/* void makeReadToPathMaps(); */

	/* void makeOneReadToPathMap( ReadToPathsMap    *end_read_paths, */
	/* 						   PathToGsaTypesMap *end_path_reads); */

	/* bool latchReadBridgingPaths( PathId spid, int direction ); */
	void stitchPath( PathId spid, int direction );

	bool latchReadBridgingPaths( PathId spid, int direction );

	void findCandidates( std::multimap<size_t, PathId> &count_map,
						 PathId spid,
						 int direction );
	
	bool tryLatch( PathId spid,
				   PathId mpid,
				   const ReadMatchPosList *pivot_poss,
				   std::tr1::unordered_map<ReadId, bool> &pivot_reads,
				   int direction );

	bool extractSharedReads( const ReadMatchPosList *pivot_poss,
							 const ReadMatchPosList *match_poss,
							 ReadPosInfoMap &pivot_rmap, 
							 ReadPosInfoMap &match_rmap, 
							 ReadFlagMap &comm_reads,
							 PathId spid,
							 PathId mpid,
							 int direction );
	
	bool determineBridgingEvidence( ReadIdSet &bridges,
									int &overlap,
									std::string &mid_str,
									std::string &pivot_str,
									std::string &match_str,
									ReadPosInfoMap &pivot_poss,
									ReadPosInfoMap &match_poss,
									int direction );

	void addReads( PathId spid,
				   PathId mpid,
				   ReadIdSet &bridges,
				   ReadPosInfoMap & pivot_poss,
				   ReadPosInfoMap & match_poss,
				   AlignSummary &summary,
				   int direction );
	


	/* bool goodRead( IntPair &poss,  */
	/* 			   GsaType &sa,  */
	/* 			   std::string &pstr,  */
	/* 			   int direction ); */

	void updateMap( PathId spid, 
					PathId mpid, 
					AlignSummary &summary,
					int direction );

	bool sameBases( const char *str1,
					const char *str2,
					int &s1,
					int &s2,
					const int &l1,
					const int &l2,
					int &len );

	int countUnrecruitedReads();
	void extractBridgeReads();
	void __extractBridgeReadsSP();
	void __extractBridgeReadsMP();

	void dropCommonPaths( PathPossMap &lsupports, 
						   PathPossMap &rsupports, 
						   const BoundType &lsuffix,
						   const BoundType &rsuffix );

	void extractValidPaths( int i,
							PathPosInfoMap &lgood_paths,
							PathPosInfoMap &rgood_paths,
						   const BoundType &lsuffix,
						   const BoundType &rsuffix );


 public:
	ReadJoinerSA();
	virtual ~ReadJoinerSA();

	/* void stitchAll(); */

	// pure virtual function from Connecter class
	void purgeTemp();

};

#endif
