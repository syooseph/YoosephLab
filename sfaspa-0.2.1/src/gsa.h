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
#include "sfa.h"


/** 
 * \brief Generalized suffix array type definition
 */
struct GsaType
{
	SfaType doc; ///< document ID
	LcpType pos; ///< position
	GsaType(){}
	GsaType(SfaType d, LcpType p) 
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
	SfaType *Ids; ///< Array of read IDs
	LcpType *Pos; ///< Array of positions in concatenated string

//=========================
// Private member functions	
//=========================
 private:

	/**
	 * Initialize GSA
	 * \param s a set of sequence reads
	 * \param n a number of reads
	 */
	void init(char **s, int n);

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




	/** Range search of given string using LCP arrays */
	BoundType searchWithLCPs( const SfaChar *srch, int len );

	/** Range search of given string without LCP arrays */
	BoundType searchOnGSA( const SfaChar *srch, int len, SfaType lmin, SfaType rmax );

	/** Get left boundary of query sequence match using LCP arrays */
	SfaType getLeftBoundWithLCPs(const SfaChar* pat, int len);

	/** Get right boundary of query sequence match using LCP arrays */
	SfaType getRightBoundWithLCPs(const SfaChar* pat, int len);

	/** Get left boundary of query sequence match without LCP arrays */
	SfaType getLeftBoundOnGSA(const SfaChar* pat, int len, SfaType lmin, SfaType rmax);

	/** Get right boundary of query sequence match without LCP arrays */
	SfaType getRightBoundOnGSA(const SfaChar* pat, int len, SfaType lmin, SfaType rmax);


	SfaType refineLeftBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax);
	SfaType refineRightBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax);

	SfaType getLeftEndBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax);	
	SfaType getRightEndBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax);	

	/** Get a sequence of given position */
	char* getSeq( SfaType p );

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
	GSA( char **s, int n );

	/** 
	 * Constructor
	 * \param s a set of sequence reads
	 * \param n a number of sequences
	 * \param f create LCPs?
	 */
	GSA( char **s, int n, bool f);

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
		 int    nseq, 
		 bool   flag, 
		 const char *gsa_file, 
		 const char *lcp_file, 
		 const char *mcp_file);

	/** Default destructor */
	~GSA();

	/** Clear all objects */
	void clear();

	/** Get a suffix of read starting from given position */ 
	char*    getSuffix(size_t) const;

	/** Get a suffix of entire string from given position */ 
	char*    getEntireSuffix(size_t);

	/** Get a record at given position */
	GsaType  getAt(size_t) const;

	/** Get a number of reads */
	int  getReadCount() { return nreads; }

	/** Set sequence reads */
	void setSequences( char **s ) { seqs   = s; }

	/** Set read count */
	void setReadCount( int n )    { nreads = n; }

	/** Set verbosity */
	void setVerbosity( bool v)    { verbose = v; }

	/** Print suffix array */
	void printSFA();

	/** Print suffix */
	void printSuffix();
	
	/** 
	 * Search a matching range of a pattern.
	 * \param srch pattern to search
	 * \param len  a length of the pattern
	 */
	BoundType search( const SfaChar *srch, int len );
	//BoundType search( const SfaChar *srch, int len, SfaType left, SfaType right );
	
	BoundType refine( const SfaChar *srch, int len, SfaType lmin, SfaType rmax );

	BoundType getEndBound( const SfaChar *srch, int len, SfaType lmin, SfaType rmax );


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

};

#endif
