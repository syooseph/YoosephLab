/** 
 * @file    gsa.h
 * @date    2013
 * @date    Modified: Fri 2013-12-20 02:15:59 PM
 * @author  Youngik Yang
 * @version 0.002
 * @brief   Generalized suffix array
 * @details SFA wrapper for multiple sequences.
 * @copyright J. Craig Venter Institute.
 */

#ifndef __GENERALIZED_SUFFIX_ARRAY_H__
#define __GENERALIZED_SUFFIX_ARRAY_H__

#include <algorithm>
#include <divsufsort64.h>
#include "sfa.h"

typedef saint_t RidType;

/** 
 * \brief Generalized suffix array type definition
 */
struct GsaType
{
	RidType doc; ///< document ID
	LcpType pos; ///< position
	GsaType(){}
	GsaType(RidType d, LcpType p) 
	{
		doc = d; pos = p;
	}
};


/**
 * \brief Generalized suffix array
 */
class GSA : public SFA
{
 private:
	RidType *Ids; ///< Array of read IDs
	LcpType *Pos; ///< Array of positions in concatenated string

//=========================
// Private member functions	
//=========================
 private:
  
  
	/**
	 * Initialize GSA
	 * \param Strings a set of sequence reads
	 * \param Integer a number of reads
	 */
	void init(char **s, RidType n);

	/** 
	 * Build GSA
	 * First, build SFA. 
	 * Then, find read ID and positions in concatenated string by
	 * calling __convertWithArrays().
	 */
	void buildSFA();

	/**
	 * Build GSA and write GSA to file. 
	 */
	void buildSFA(const char*);

	/** 
	 * Alternative read ID and position search using binary search.
	 * It is more memory efficient but slow.
	 */
	void __convertWithBinarySearch();

	/**
	 * Read ID and position search from concatenated string.
	 * It is fast but requires more RAMs 
	 * because of two additional arrays to keep track of reads and IDs.
	 */
	void __convertWithArrays();

	/** Build LCP */
	void buildLCPs();
	
	/** Build LCP and internal LCP 	*/
	void buildLCPs(const char*, const char*);


	/** Range search of given string without LCP arrays */
	BoundType searchOnGSA( const SfaChar *srch, LcpType len );
	
	/** Get left boundary of query sequence match using LCP arrays */
	IdxType getLeftBoundWithLCPs(const SfaChar* pat, LcpType len);

	/** Get right boundary of query sequence match using LCP arrays */
	IdxType getRightBoundWithLCPs(const SfaChar* pat, LcpType len);

	/** Get left boundary of query sequence match without LCP arrays */
	IdxType getLeftBoundOnGSA(const SfaChar* pat, LcpType len);

	/** Get right boundary of query sequence match without LCP arrays */
	IdxType getRightBoundOnGSA(const SfaChar* pat, LcpType len);
	
	/** Get a sequence of given position */
	char* getSeq( IdxType p );

	/** Write generalized suffix array */
	void writeSFA( const char *filename );

	/** Read generalized suffix array from a file */
	void readSFA( const char *filename );


//=========================
// Public member functions	
//=========================
 public:
	/** Default constructor */
	GSA();

	/** 
	 * Constructor
	 * \param s a set of sequence reads
	 * \param n Integer a number of sequences
	 */
	GSA( char **s, RidType n );

	/** 
	 * Constructor
	 * \param s a set of sequence reads
	 * \param n a number of sequences
	 * \param f create LCPs?
	 */
	GSA( char **s, RidType n, bool f);

	/** 
	 * Constructor
	 * \param seqs a set of sequence reads
	 * \param nseq a number of sequences
	 * \param flag create LCPs?
	 * \param gsa_file GSA output file name
	 * \param lcp_file LCP output file name
	 * \param mcp_file internal LCP output file name
	 */
	GSA( char **seqs, 
		 RidType    nseq, 
		 bool   flag, 
		 const char *gsa_file, 
		 const char *lcp_file, 
		 const char *mcp_file);

	/** Default destructor */
	~GSA();

	/** Clear all objects */
	void clear();

	/** Get a suffix of read starting from given position */ 
	char*    getSuffix(IdxType);

	/** Get a suffix of entire string from given position */ 
	char*    getEntireSuffix(IdxType);

	/** Get a record at given position */
	GsaType  getAt(IdxType);

	/** Get a number of reads */
	RidType  getReadCount() { return nreads; }

	/** Set sequence reads */
	void setSequences( char **s ) { seqs   = s; }

	/** Set read count */
	void setReadCount( RidType n )    { nreads = n; }

	/** Set verbosity */
	void setVerbosity( bool v)    { verbose = v; }

	/* /\** */
	/*  * Build generalized suffix array */
	/*  * \param  */
	/*  *\/ */
	/* void buildSFA(bool); */

	/** Print suffix array */
	void printSFA();

	/** Print suffix */
	void printSuffix();
	
	/** 
	 * Search a matching range of a pattern.
	 * \param srch pattern to search
	 * \param len  a length of the pattern
	 */
	BoundType search( const SfaChar *srch, LcpType len );

	/**
	 * Load binary objects
	 * \param lcp_file LCP file
	 * \param mcp_file internal LCP file
	 * \param gsa_file GSA file
	 */
	void load( const char* lcp_file, 
			   const char* mcp_file, 
			   const char* gsa_file);

	/**
	 * Load binary objects
	 * \param sfa_file SFA file
	 * \param con_file concatenated string file
	 * \param lcp_file LCP file
	 * \param mcp_file internal LCP file
	 * \param gsa_file GSA file
	 */
	void load( const char* sfa_file, 
			   const char* con_file, 
			   const char* lcp_file, 
			   const char* mcp_file, 
			   const char* gsa_file);

	/**
	 * Dump binary objects
	 * \param sfa_file SFA file
	 * \param con_file concatenated string file
	 * \param lcp_file LCP file
	 * \param mcp_file internal LCP file
	 * \param gsa_file GSA file
	 */
	void dump( const char* sfa_file, 
			   const char* con_file, 
			   const char* lcp_file, 
			   const char* mcp_file, 
			   const char* gsa_file);
			   
  /********************Newly added or redefined************************************/
  /** Range search of given string using LCP arrays */
	BoundType searchWithLCPs( const SfaChar *srch, LcpType len );
	char* getSuffix_explicit(RidType rid, LcpType pos);
	char* getSequence_explicit(RidType rid);
	LcpType  getSuffixLength(IdxType p);
	LcpType  getFullSequenceLength(IdxType p);

	RidType getId(IdxType p);
  LcpType getPos(IdxType p);	
  LcpType getLcpArbitrary(IdxType p, IdxType q);
  
  BoundType searchWithLCPs_bounded(
      const SfaChar *srch, LcpType len, 
      IdxType left_bound, IdxType right_bound
  );
	IdxType getLeftBoundWithLCPs_bounded(
	    const SfaChar* pat, LcpType len, 
	    IdxType left_bound, IdxType right_bound
	);
	IdxType getRightBoundWithLCPs_bounded(
	    const SfaChar* pat, LcpType len, 
	    IdxType left_bound, IdxType right_bound
	);
	
	LcpType getLcp(IdxType p);
	LcpType getSeqLength_RID(RidType rid);

};

#endif
