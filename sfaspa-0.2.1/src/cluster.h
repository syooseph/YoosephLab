/** 
 * @file       cluster.h
 * @brief      Cluster declaration
 * @date       Modified on Fri 2014-04-11 12:59:30 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include "extracter.h"

/** 
 * a list of alignment summaries 
 */
typedef std::list<AlignSummary> AlignSummaryList;

/**
 * \class    Cluster
 * \brief    Cluster contents
 * \deftails This struct defines cluster center sequence, member paths, 
 *           alignment summaries, etc.
 */
struct Cluster
{
	std::string sequence;       ///< cluster center sequence
	PathIdList  members;        ///< a list of member paths
	AlignSummaryList summarys;  ///< alignment summaries to members
	int lstop;                  ///< left stop condition
	int rstop;                  ///< right stop condition
	int ltrim;                  ///< left trimmed bases
	int rtrim;                  ///< right trimmed bases
	int lbase;                  ///< left new extra bases
	int rbase;                  ///< right new extra bases

	/** add new member */
	void add(PathId &p, AlignSummary &a );

	/** write to a file */
	void dump( std::fstream &out );

	/** load from a file */
	void load( std::fstream &in );

	/** 
	 * Find maximum end (head or tail) gaps from members.
	 * @param g maximum gap to be reported
	 * @param p maximum pid to be reported
	 * @param l left extension?
	 */
	void findMaxGap( size_t &g, PathId &p, bool l);

	/** 
	 * Update leading gap and indel positions 
	 */
	void updatePositions( size_t &max_lgap );
};


#endif
