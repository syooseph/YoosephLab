/** 
 * \file      sjoiner.h
 * \brief     Short overlapping path extension
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \bug       None.
 * \warning   None.
 * \date      Fri 2014-01-043
 * \copyright J. Craig Venter Institute.
 */

#ifndef __SHORT_JOINER_H__
#define __SHORT_JOINER_H__

#include "connecter.h"

/**
 * \class TinyJoiner
 * \brief Path joiner for short overlapping paths.
 */
class TinyJoiner : public Connecter
{
 private:
	EndStringMap  lend_map; ///< path to head strings map
	EndStringMap  rend_map; ///< path to tail strings map

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
	 * Find overlapping paths by checking sequence tips.
	 */
	void findPathsByEndString( std::vector<PathId> &matches,
							   std::vector<size_t> &lengths,
							   PathId spid, 
							   std::string &pivot,
							   std::tr1::unordered_map<PathId, bool> checked,
							   int direction );
	
	
	/** 
	 * Join a path if possible.
	 */
	bool tryConnectOverlap( size_t length,
							PathId &mpid,
							PathId &spid, 
							std::string &pivot, 
							std::tr1::unordered_map<ReadId, bool> &preads,
							int direction );
	
	/** 
	 * Elongate of overlapping sequence in both direction.
	 */
	bool getLadderString( std::string &ladder,
						  size_t length,
						  std::string &pivot,
						  std::string &match,
						  int direction );

	/**
	 * Update path-maps of extended path & merged path
	 */
	void updateIndex( PathId spid,
					  std::string &old_pivot,
					  std::string &new_pivot,
					  int direction);

	/**
	 * Update entries of a path with given string.
	 */
	void updateEndString( std::string &str,
						  PathId pid );


	/** 
	 * Drop entries of a path with given string.
	 */
	void deleteEndString( std::string &str,
						  PathId pid );
	
 public:
	TinyJoiner();
	virtual ~TinyJoiner();

	/** Release temporary storages */
	void purgeTemp();
};

#endif
