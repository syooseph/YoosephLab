/** 
 * \file      rjoiner.h
 * \brief     Path extension by bridging reads
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \bug       None.
 * \warning   None.
 * \date      Tue 2014-02-04 05:22:50 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __READ_JOINER_H__
#define __READ_JOINER_H__

#include "connecter.h"

typedef std::tr1::unordered_map<ReadId, PathId> ReadPathMap;
typedef std::list<GsaType> GsaTypeList;
typedef std::tr1::unordered_map<PathId, GsaTypeList> PathToGsaTypesMap;
typedef std::tr1::unordered_map<ReadId, PathIdSet> ReadToPathsMap;
typedef std::unordered_map<ReadId, IntPair> ReadPathPoss;

/**
 * \class ReadJoiner
 * \brief Path joiner for bridging read support
 */
class ReadJoiner : public Connecter
{
 private:
	/* PathToGsaTypesMap   lend_path_reads; */
	/* PathToGsaTypesMap   rend_path_reads; */

	/* ReadToPathsMap      lend_read_paths; */
	/* ReadToPathsMap      rend_read_paths; */

 protected:
	ReadPathMap     recruited;
	BridgeLog       bridge_log;
	
 protected:
	void initMaps();

	/* void getGsaTypes( GsaTypeList &gsas, */
	/* 				  std::string &end_str ); */

	/* virtual bool tryLatch( PathId spid, */
	/* 					   PathId mpid, */
	/* 					   const GsaTypeList *pivot_gsas, */
	/* 					   int direction ) = 0; */
	
	/* virtual bool extractSharedReads( const GsaTypeList *pivot_gsas, */
	/* 								 const GsaTypeList *match_gsas, */
	/* 								 std::tr1::unordered_map<ReadId,GsaType> &pivot_rmap,  */
	/* 								 std::tr1::unordered_map<ReadId,GsaType> &match_rmap,  */
	/* 								 std::tr1::unordered_map<ReadId, bool> &comm_reads, */
	/* 								 PathId spid, */
	/* 								 PathId mpid, */
	/* 								 int direction ); */

	/* bool dropPartialReads( std::tr1::unordered_map<ReadId,GsaType> &pivot_rmap,  */
	/* 					   std::tr1::unordered_map<ReadId,GsaType> &match_rmap,  */
	/* 					   std::tr1::unordered_map<ReadId, bool> &comm_reads, */
	/* 					   ReadPathPoss &pivot_poss,  */
	/* 					   ReadPathPoss &match_poss, */
	/* 					   std::string &pivot_str, */
	/* 					   std::string &match_str, */
	/* 					   int direction ); */
	
	/* bool goodRead( IntPair &poss, */
	/* 			   GsaType &sa, */
	/* 			   std::string &pstr, */
	/* 			   int direction ); */

	/* bool determineBridgingEvidence( ReadIdSet &bridges, */
	/* 								int &overlap, */
	/* 							  std::string &mid_str, */
	/* 							  std::string &pivot_str, */
	/* 							  std::string &match_str, */
	/* 							  ReadPathPoss &pivot_poss,  */
	/* 							  ReadPathPoss &match_poss, */
	/* 							  int direction ); */
	
	/**
	 * Base class function.
	 * Find overlapping paths and join them.
	 */
	bool latchOverlapPaths( PathId &spid, int direction );

	/* void addReads( PathId spid, */
	/* 			   PathId mpid, */
	/* 			   ReadIdSet &bridges, */
	/* 			   ReadPathPoss &pivot_poss, */
	/* 			   ReadPathPoss &match_poss, */
	/* 			   AlignSummary &summary, */
	/* 			   int direction ); */
	
	/* void updateMap( PathId spid, PathId mpid, int direction ); */

	/**
	 * Base class function 
	 */
	void updateIndex( PathId spid,
					  std::string &old_pivot,
					  std::string &new_pivot,
					  int direction);
	
	void updateSequence( PathId lpid,
						 PathId rpid,
						 int overlap,
						 std::string &mid_str,
						 int direction );
	
	void adjustAlignPositions( AlignSummary &summary, 
							   const int &overlap, 
							   const std::string &mid_str, 
							   const PathId &spid, 
							   const PathId &mpid, 
							   const int &direction );
	

	size_t getAddedReadsCount();

 public:
	ReadJoiner();
	virtual ~ReadJoiner();

	// function overloading
	void run();

	/** Release temporary storages */
	void purgeTemp();

	void printAddedReads();

	
};

#endif
