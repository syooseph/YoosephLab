/** 
 * @file rplace.h
 * @date Tue 2013-12-17 02:43:37 PM
 * @author Youngik Yang
 * @version 0.2
 * @brief Read placement.
 * @details 
 * Single read placement to a path.
 * @copyright J. Craig Venter Institute.
 */

#ifndef __READ_PLACEMENT_H__
#define __READ_PLACEMENT_H__

#include "path.h"
#include <list>
#include <fstream>

/**
 * \brief Single read placement to a path
 */
struct ReadPlacement
{
	ReadId rid;                    ///< Read ID
	int read_pos;                  ///< Read alignment start positioin
	int path_pos;                  ///< Alignment position in a path
	int length;
	std::list<Mismatch> ilist;     ///< Insertions
	std::list<Mismatch> dlist;     ///< Deletions

	void print(std::ostream &out); ///< Print placement
	void dump(std::fstream &out);  ///< Dump binary object
	void load(std::fstream &in);   ///< Load binary object
};

typedef std::list<ReadPlacement> ReadPlacementList;

#endif
