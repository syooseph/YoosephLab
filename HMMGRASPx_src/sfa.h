/** 
 * @file    sfa.h
 * @date    2013
 * @date    Modified: Fri 2013-12-13 06:44:49 PM
 * @author  Youngik Yang
 * @version 0.002
 * @brief   Suffix array
 * @pre     divsufsort needs to be pre-installed.
 * @warning Data type is 32 bit signed integer.
 * @copyright J. Craig Venter Institute.
 */

#ifndef __SUFFIX_ARRAY_H__
#define __SUFFIX_ARRAY_H__

#include <stack>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <divsufsort64.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <cassert>
#include <sstream>
#include <limits>
#include "timer.h"
#include "file.h"
#include "interval_array.h"


typedef sauchar_t SfaChar;
typedef sauchar_t LcpType;
typedef saidx64_t IdxType; 

//#if BIGINT == 0
//typedef saidx_t   SfaType; 
//#else
//typedef saidx64_t   SfaType;       
//#endif

typedef std::pair<IdxType, IdxType> BoundType;

const char DELIMITER = '$';

class SFA
{
 protected:
	SfaChar       *concat;  ///< concatenated string
	IdxType       *SA;      ///< suffix array
	LcpType       *LCP;     ///< longest common prefix array
	LcpType       *mLCP;    ///< internal LCPs
	IdxType        size;     ///< size of suffix array
	char          **seqs;   ///< sequence reads
	int           nreads;   ///< # of reads
	//IntervalArray *ia;      ///< interval array
	IntervalArray ia;
	bool          verbose;  ///< verbosity

//=========================
// private member functions 
//=========================
 private:
	/** initialization */
	void init();               
 
	/** 
	 * Determine the size of the suffix array.
	 * Raise an error and terminate the program if the size is too big to load.
	 */
	void setStorageSize();     
	
	/** Concatenate all sequence reads */
	void concatenateAllSeqs();

	/** Clear all */
	void clear();

	/** Range search of given string using two precomputed LCP arrays. */
	BoundType searchWithLCPs( const SfaChar *srch, int len );

	/** Range search of a string without LCP arrays. */
	BoundType searchOnSFA( const SfaChar *srch, int len );

	/** Find left match boundary of query sequence. */
	IdxType getLeftBoundWithLCPs(const SfaChar* pat, int len);

	/** Find right match boundary of query sequence. */
	IdxType getRightBoundWithLCPs(const SfaChar* pat, int len);


	/** Find left match boundary of query sequence. */
	IdxType getLeftBoundOnSFA(const SfaChar* pat, int len);

	/** Find right match boundary of query sequence. */
	IdxType getRightBoundOnSFA(const SfaChar* pat, int len);
	
//===========================
// protected member functions 
//===========================
 protected:
	//==========
	// Clearners
	//========== 
	

	//========
	// Writers
	//========
	/** Write suffix array to a binary file */
	void writeSFA( const char *filename );

	/** Write concatenated string to a binary file */
	void writeDoc( const char *filename );

	/** Write LCP to a binary file */
	void writeLCP( const char *filename );

	/** Write internal LCP to a binary file */
	void writeMLCP( const char *filename );

	//========
	// Readers
	//========	
	/** Read suffix array from a file */
	void readSFA( const char *filename );

	/** Read concatenated string from a file */
	void readDoc( const char *filename );

	/** Read LCP from a file */
	void readLCP( const char *filename );

	/** Read internal LCP from a file */
	void readMLCP( const char *filename );

	//==================
	// Utility functions
	//==================
	/** Find the longest commont prefix of two strings */
	int lcp(const SfaChar* a, const SfaChar* b);

	/** Get internal LCP of given range */
	// Return value changed to LCPType
	LcpType getMcp(IdxType left, IdxType right );


//========================
// Public member functions
//========================
 public:
	/** Empty default constructor */
	SFA(); 

	/** 
	 * Constructor
	 * \param seqs sequences
	 * \param nseq number of sequences
	 */
	SFA(char **seqs, int nseq);

	/** Empty destructor */
	~SFA();

	//========
	// Getters
	//========
	/** Get suffix array */
	IdxType*  getSFA() { return SA; }

	/** Get concatenated string */
	SfaChar*  getDoc() { return concat; }

	/** Get the size of suffix array */
	IdxType   getSize() { return size; }

	/** Get the suffix starting at a given position */
	char*     getSuffix(IdxType p);

	/** Get the suffix array value of given positioin */
	IdxType   getAt(IdxType p);

	/** Get LCP value at position i */
	LcpType getLcpAt(IdxType i) { return LCP[i]; }

	/** Get internal LCP value at position i */
	LcpType getMcpAt(IdxType i) { return mLCP[i]; }


	//========
	// Setters
	//========
	/** Set a set of sequences */
	void setSequences(char **s) { seqs    = s; }

	/** Set read count */
	void setReadCount(int    n) { nreads  = n; }

	/** Set verbose option */
	void setVerbosity(bool   v) { verbose = v; }


	//=========
	// Builders
	//=========
	/** Construct SFA */
	void buildSFA();

	/** Construct LCP */
	void buildLCP();

	/** Construct internal LCP to make fast query */
	void buildMLCP();


  /** Delete suffix array */
	void purgeSA();

	/** Delete concatenated string */
	void purgeDoc();

	/** Delete LCP */
	void purgeLCP();

	/** Delete internal LCP */
	void purgeMLCP();
	//=============
	// Print values
	//=============
	/** Print suffix array */
	void printSFA();    

	/** Pint concatenated string */
	void printDoc();

	/** Print longest common strings of consecutive suffixes */
	void printLCP(); 

	/** Print longest common strings of internal nodes */
	void printMLCP();  

	/** Print all suffixes */
	void printSuffix();


	//=================
	// Search functions
	//=================
	/** 
	 * Simple binary search of string match.
	 * All matching positions are reported.
	 */
	void search( std::vector<IdxType> &pos,
                  const SfaChar *srch, 
                  int len );

	/** 
	 * Wrapper for string matching range search.
	 * A range in SFA is reported.
	 */
	BoundType search( const SfaChar *srch, int len );

	//=============
	// File options
	//=============
	/** Dump objects as binaries */
	void dump( const char *sfa_file, 
			   const char *doc_file, 
			   const char *lcp_file, 
			   const char *mcp_file );

	/** Load binary objects */
	void load( const char *sfa_file, 
			   const char *doc_file, 
			   const char *lcp_file, 
			   const char *mcp_file );

	
};

#endif
