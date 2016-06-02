/** 
 * @file       IntervalArray.h
 * @brief      Quick interval lookup
 * @details    Internal LCP implementation
 * @date       Modified on Tue 2013-12-17 06:46:57 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __INTERVAL_ARRAY_H__
#define __INTERVAL_ARRAY_H__

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <stack>
#include <cassert>
//#include <boost/unordered_map.hpp>
#include <tr1/unordered_map>
#include "file.h"

//===========================================
// NOTE:
// ValueType must be compatible with LCP type
//===========================================
typedef unsigned char ValueType;

//================================================
// Overflow map: LCP is bigger than char.
// This feature will be implemented in the future.
//================================================
typedef std::tr1::unordered_map<size_t, size_t> BigValueMap;

const ValueType MAX_INTERVAL_VALUE = 255;

/**
 * \class IntervalArray
 * \brief Internal LCP value look up.
 */
class IntervalArray
{
 private:
	size_t size; ///< size of internal LCPs
	size_t leaf; ///< size of lead LCPs

	ValueType *leaf_ints;  ///< pointer of leaf LCP values
	ValueType *intervals;  ///< internal LCP values

 private:
	/** 
	 * Determine index of given range 
	 * \param l left
	 * \param r right
	 */
	long getIndex( size_t l , size_t r );

	/** 
	 * Recurively fill the interval array of given range.
	 * Fill the minimum of left and right sub-trees.
	 * \param l left
	 * \param r right
	 */
	ValueType fill( size_t l, size_t r );

 public:
	IntervalArray();
	/**
	 * Constructor
	 * \param ints Leaf LCPs
	 * \param n size of leaf LCPs
	 */
	IntervalArray(ValueType *ints, size_t n);
	~IntervalArray();

	/**
	 * Build interval LCPs
	 * \param ints Leaf LCPs
	 * \param n size of leaf LCPs
	 */
	void build(ValueType *ints, size_t n);

	/**
	 * Release memory
	 */
	void clear();

	/**
	 * Print all internal LCPs
	 */
	void print( std::ostream &out );

	/** 
	 * \return size of internal LCPs
	 */
	size_t   getSize() { return size; }

	/**
	 * \return interval arrays
	 */
	ValueType* getIntervals() { return intervals; }

	/**
	 * \return LCP value in range
	 * \param i Leaf LCPs
	 * \param j size of leaf LCPs
	 */
	ValueType  getValue(size_t i, size_t j);

	/** Dump the array to a file */
	void dump( const char *filename );

	/** Load an array from a file */
	void load( const char *filename, size_t filesize );
};


#endif
