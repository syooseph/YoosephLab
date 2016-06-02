/** 
 * @file       loader.h
 * @brief      Objects loader.
 * @details    This loads sequences and suffix arrays prior to assembly.
 * @date       Modified on Tue 2013-12-17 06:46:57 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __LOADER_H__
#define __LOADER_H__

#include "core.h"

const ReadId NOT_PAIR = 4294967295; // 2^32-1
static const double INIT_TIME = mytime();

typedef std::vector<int> IntArray;

/* Strands */
enum strands { POSITIVE, /**< positive strand. */ 
			   NEGATIVE, /**< negative strand. */ 
			   NONSENSE  /**< otherwise        */ 
};                                     

/**
 * \brief Object loaders prior to assembly.
 */
class Loader
{
 private:
	int             nreads;        // Total number of reads
	char            **seqs;        // sequence reads
	char            **tags;        // sequence names
	GSA             *gsa;          // Generalized suffix array
	int             nparts;        // No. of suffix array partitions
    ReadId          *pairs;        // Read pairing
    char            *strands;      // Read strands
	IntArray        bad_reads;     // Bad sequences
	bool            status;        // Ojbect loaded?
	int             max_len;       // max read len
	int             min_len;       // min read len
	
 private:
	/** Load sequence reads */
	void loadReads();

	void getLengths();

	/** Mark bad sequences awary */
	void trimReads();

	/** Load suffix arrays */
	void loadSuffixArray();

	/** Set sequence strands */
	void setStrands(char**, int);

	/** Set pair-end reads */
	void setReadPairs(char **, int, bool);

	/** Remove pair-info after consuming pair info */
	void dropPairInfo(char **, int);

 public:
	Loader();
	~Loader();

	/** Load all objects */
	void   load();

	/** Clear all objects */
	void   clear();

	//========
	// Getters
	//======== 
	int       getCount()    { return nreads; }
	char**    getReads()    { return seqs;   }
	char**    getNames()    { return tags;   }
	GSA*      getGSA()      { return gsa;    }
	int       getParts()    { return nparts; }
	ReadId*   getPairs()    { return pairs; }
	char*     getStrands()  { return strands; }
	IntArray* getBadReads() { return &bad_reads; }
	int       getMaxLength(){ return max_len; }
	int       getMinLength(){ return min_len; }

	/** Release generalize suffix array */
	void purgeGSA();

	bool getStatus() { return status; }

};
#endif

