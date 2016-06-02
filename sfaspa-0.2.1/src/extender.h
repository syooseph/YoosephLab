/** 
 * \file      extender.h
 * \brief     Extended path object
 * \details   This keeps track of aligned position of joined paths. 
 *            The extend object updates representative sequence, alignment summaries 
 *            as well as path end criteria.
 *            
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \bug       None.
 * \warning   None.
 * \date      Modified on Fri 2013-12-13 02:10:26 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __EXTENDER_H__
#define __EXTENDER_H__

#include "extracter.h"

/** PathId to aligned position map */
typedef std::tr1::unordered_map<PathId, int> StartPosMap;
/** PathId to alignment summary map */
typedef std::tr1::unordered_map<PathId, AlignSummary> PathAlignMap;

/**
 * \class Extender
 * \brief An extended entry during path extension
 */
struct Extender
{
	std::string sequence;     ///< extended path sequence
	PathIdSet members;        ///< newly joined paths
	StartPosMap positions;    ///< aligned position
	PathAlignMap aligns;      ///< alignments between pivot and joined paths
	int lstop,rstop;          ///< stop condition
	int ltrim,rtrim;          ///< trim condition

	Extender();
	Extender(std::string s);
	~Extender();

	/** intialization */
	void init();

	/** binaray dump */
	void dump(std::fstream &out);

	/** read binary object */
	void load(std::fstream &in);
};

/** Extened path ID */
typedef size_t LatchId; 
/** Extended path ID and extended path object map */
typedef std::tr1::unordered_map<LatchId, Extender> ExtenderMap;

#endif
