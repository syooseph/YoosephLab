/**
 * \file      rcmat.h
 * \brief     Read count matrix
 * \details   Read counts of from minimum length sub-path (k+2 mer) to 
 *            maximum length sub-path (n-back mer) are stored/updated 
 *            in a matrix.
 * \author    Youngik Yang
 * \version   0.001
 * \date      Tue 2013-12-17 06:25:06 PM
 * \date      Modified on Mon 2014-02-24 12:39:04 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __RANGE_COUNT_MATRIX__
#define __RANGE_COUNT_MATRIX__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <vector>
#include <iostream>
#include "core.h"


typedef std::vector<std::vector<int> > Matrix;

/**
 * \brief Count of supporting reads of k+2-mer to nback-mer in a given path
 */
class RangeCountMatrix
{
 private:
	int nrow;
	int ncol;
	Matrix matrix;

 private:


 public:
	RangeCountMatrix();
	RangeCountMatrix(int r, int c);
	RangeCountMatrix(int len, int min, int max);
	~RangeCountMatrix();

	void clear();
	void init(int len, int min, int max) ;
	void init(int r, int c);
	void increment(int r, int b, int e);
	void increment(int pos, int len, int min, int max);

	int get(int row, int col);
	int get(int left, int right, int min, int max);

	int getRow() { return nrow; }
	int getCol() { return ncol; }

	void trim(int s, int e);
	void print();
};
#endif
