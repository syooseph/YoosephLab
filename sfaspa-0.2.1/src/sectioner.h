/** 
 * @file       sectioner.h
 * @brief      Range finder of two sequences
 * @details    Matching (merging/extension) region is searched based on kmer matches.
 * @date       Thu 2014-01-16 01:21:25 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __SECTIONER_H__
#define __SECTIONER_H__

#include "core.h"
#include "PathKmerPosIndex.h"

//typedef std::tr1::unordered_map<KmerId, std::vector<int> > KmerPossMap;
typedef std::multimap<int,int> PosPairMap;

/**
 * @brief Match ranges of a pair of sequences
 */
struct Section
{
	int qbeg, qend; ///< query start, end
	int sbeg, send; ///< sbjct start, end
	int count;      ///< match count

	Section()
	{
		count = 0;
	}

	Section( int qs, int qe, int ss, int se, int c ) 
	{
		qbeg = qs; qend = qe;
		sbeg = ss; send = se;
		count = c;
	}
	
	bool operator<( const Section & other ) const
	{
		if ( count < other.count ) return true;
		return false;
	}

	void clear()
	{
		qbeg = qend = sbeg = send = count = 0;
	}
};


/**
 * @brief Similar region detection of a pair of sequences.
 */
class Sectioner
{
 protected:
	std::string *query;
	std::string *sbjct;
	/* KmerArray  query_kmers; ///< kmers from query sequence */
	/* KmerArray  sbjct_kmers; ///< kmers from sbjct sequence */
	Section    section;     ///< match  range
	bool verbose;

 protected:
	/** Initalization */
	//void init( KmerArray &qkmers, KmerArray &skmers );
	void init( std::string *query, std::string *sbjct );

 public:
	//============================
	// Constructors/Initialization
	//============================
	Sectioner();
	Sectioner( std::string *query, std::string *sbjct );

	//=======
	// Setter
	//======= 
	void setVerbosity( bool v ) { verbose = v; }

	//=======
	// Getter
	//=======
	Section getSection() { return section; }

	/** Find sub-area */
	virtual bool find() = 0;
};

#endif

