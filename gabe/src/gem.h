/*
 Youngik Yang 
 Fri 2015-06-19 09:05:50 PM
*/

#ifndef __GEM_H__
#define __GEM_H__

#include <stdint.h>
#include <stdio.h>
#include <cmath>
#include <numeric>
#include <cassert>
#include <vector>
#include <list>
#include <iostream>
#include <thread>
#include <iomanip>
#include <algorithm>
#include "timer.h"
#include <unordered_map>
#include <cfloat>
#include "matrix.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif


const size_t  MAX_ITER  = 1.0e+5;
const double  TOLERANCE = 1.0e-6;
const double  PSEUDO    = 1.0e-12;


/* using EntryMap = std::unordered_map<uint32_t, double>; */
/* using SparseEntryMap = std::vector<EntryMap>; */

class GEM
{
 private:
	/* Inputs */
	MatrixType *matrix;
	size_t   nrow;      // no. of row
	size_t   ncol;      // no. of different column
	size_t   max_iter;  // max iteration
	size_t   ncpu;      // no. of CPUs
	double   tolerance; // tolerance of likelihood change
	double   pseudo;    // small number added to matrix to avoid zero possibility
	bool     even;      // even mixing coefficient during initialization 
	bool     verbose;   // verbose option

	/* EM stuffs */
	size_t iter;
	double likelihood; // record complete log likelihood
	double likelidiff;
	double   merr;      // mean absolue errorerror between two consecutive iterations
	double   rmse;
	double *mixing;    // mixing coefficients
	double *mixold;    // mixing coefficients of previous iteration
	bool converge;               // convergency
	MatrixType Q; // joint probability of mixing cofficient and count
	MatrixType R; // posterior probability (responsibilites)

	/* Temporary stuffs */
	//double *LogLiks; // loglikelihood of each read (only for multi-threaded execution)

	std::vector<std::vector<bool>> empties; //(ncpu, std::vector<bool>(ncol, true));

	double t_exp, t_max, t_con;
	
	/* Internal functions */
	void allocate();
	void performEM();
	void initMixingCoefficient();
	void doExpectation();
	void doMaximization();
	void checkConvergency();
	void showIterationSummary( double t0 );
	double computeLikelihood();
	//double computeError();
	double computeMeanAbsError();
	double computeRMSE();	

	void getEmptyColumns( size_t w, size_t r );
		
 public: 
	GEM();
	GEM( MatrixType *m, size_t nd, size_t nc, size_t cpu=1,size_t it=MAX_ITER, double to=TOLERANCE, double psu=PSEUDO,  bool eq=false, bool v=false);
	~GEM();
	void run();

	/* Getters */
	bool   converged()            { return converge; }
	double *getCoefficient()      { return mixing; } 
	double getLogLik()            { return likelihood; }
	double getError()             { return merr; } 
	size_t getIteration()         { return iter; }
	MatrixType *getPosterior()  { return &R; }

	/* void computeAbundance( double *a, uint32_t *l);	 */
};

#endif
