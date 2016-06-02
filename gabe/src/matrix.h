#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <random>
#include <vector>
#include <unordered_map>
#include <string.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>


using namespace boost::numeric::ublas;

#if SPARSE == 1
typedef compressed_matrix<double> MatrixType;
#else
typedef matrix<double> MatrixType;
#endif

class Matrix
{
 private:
	MatrixType matrix;
	size_t   nrow;   /** N */
	size_t   ncol;   /** M */
	
	std::vector<std::string> colnames;
	
 public:
	Matrix();
		
	Matrix( const char *csv, char delim);

	~Matrix();
	
	void init();

	void load( const char *data, char delim );
	
	/** Parse first line and get column names */
	void parseHeader( std::string &, char );
	
	size_t getRowSize() { return nrow; }
	size_t getColSize() { return ncol; }
	std::vector<std::string> getColNames() { return colnames; }
	
	/** Return matrix */
	MatrixType *getMatrix() { return &matrix; }
	
	/** Sampling with replacement */
	void resample( MatrixType *m, size_t n, double seed=0 );

	/** Subsampling (no replacement) */
	void subsample( MatrixType *m, size_t n, double seed=0 );

};

#endif
