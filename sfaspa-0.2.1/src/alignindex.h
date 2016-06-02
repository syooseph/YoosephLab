/**
 * @file alignindex.h
 * @brief Pair-wise alignment position.
 * @author    Youngik Yang
 * @version   0.2
 * @date      2010-2013
 * @date      Modfied on Thu 2014-01-16 01:58:56 PM
 * @copyright J. Craig Venter Institute.
 */

#ifndef __ALIGN_INDEX_H__
#define __ALIGN_INDEX_H__

#include "file.h"

/** 
 * \class AlignIndex
 * \brief Pair-wise alignment position
 */
struct AlignIndex
{
	int aln_pos; ///< alignment position
	int qry_pos; ///< query position
	int ref_pos; ///< sbjct position

	/**
	 * Empty default constructor
	 */
	AlignIndex() {}

	/**
	 * Constructor
	 * \param a alignment position
	 * \param q query position
	 * \param s sbjct position 
	 */
	AlignIndex( int a, int q, int s )
	{ 
		aln_pos = a; qry_pos = q; ref_pos = s; 
	}

	/** 
	 * Dump a binary object to file 
	 * \param out fstream writer
	 */
	void dump( std::fstream &out ) 
	{
		out.write( (char*)&aln_pos, sizeof(int) );
		out.write( (char*)&qry_pos, sizeof(int) );
		out.write( (char*)&ref_pos, sizeof(int) );
	}

	/** 
	 * Load a binary object
	 * \param out fstream reader
	 */
	void load( std::fstream &in ) 
	{
		in.read( (char*)&aln_pos, sizeof(int) );
		in.read( (char*)&qry_pos, sizeof(int) );
		in.read( (char*)&ref_pos, sizeof(int) );
	}

	/** 
	 * Print positions 
	 * \param out fstream writer
	 */
	void print( std::ostream &out ) 
	{
		out << "align:" << aln_pos << " query:" << qry_pos << " sbjct:" << ref_pos << "\n";
	}
};

#endif
