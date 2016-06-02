/** 
 * \file      anchor.h
 * \brief     Find anchor kmers of two sequences and calculate simiarity.
 * \details   Two anchor kmers are scanned from each end of a sequence. 
 *            Distance between anchor points are must be same in query and sbjct. 
 *            Then, the anchor region is extended until either sequence is ended.
 *			  Sequence similarity score is computed in the region only by 
 *			  comparing base by base. 
 *			  NO gap is allowed.
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \date      Mondified on Thu 2013-10-24 01:03:19 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __ANCHOR_ALIGN_H__
#define __ANCHOR_ALIGN_H__

#include "core.h"
#include "path.h"
#include "sectioner.h"
#include "log.h"
#include "PathKmerPosIndex.h"

typedef std::map<IntPair, bool> IntPairMap;
typedef std::tr1::unordered_map<size_t, std::vector<size_t> > PosPosMap;

/**
 * \class Anchor
 * \brief Find anchor kmers of two sequences and calculate base-to-base simiarity
 */
class Anchor : public Sectioner
{
 private:
	PathId sbjct_pid;
	ReadId query_rid;

	//std::string strim;    ///< sbjct in anchor region
	//std::string qtrim;    ///< query in anchor region

	//KmerArray *skmers;    ///< sbjct kmers
	const KmerPosMap *skposs;
	const KmerArray  *qkmers;    ///< query kmers
	
	PosPairMap poss_pair; ///< kmer match position pairs
	//PosPosMap poss_pair;
	IntPairMap checked;   ///< previously checked region in sbjct 

	int    kmer;          ///< size of kmer
	int    min_filter;    ///< minimum kmer match
	size_t length;        ///< extended anchor region length
	double score;         ///< alignment score
	bool   found;         ///< anchor search status
	bool   verbose;       ///< verbosity

	AnchorLog log;

 private:
	
	/* /\**  */
	/*  * Initialize */
	/*  *\/ */
	/* void init( std::string *s,  */
	/* 		   std::string *q,  */
	/* 		   //KmerArray *sk,  */
	/* 		    KmerPosMap *sk, */
	/* 		   KmerArray *qk, */
	/* 		   int m, */
	/* 		   int k,  */
	/* 		   bool v ); */


	/**
	 * Fill the position of kmer matches between query and sbjct
	 */
	void initPosMap();

	/**
	 * Find extended anchor region
	 */
	bool findAnchor();

	/** 
	 * Find a range from positions in query and sbjct
	 */
	bool  findRange( int qpos, int spos );

	/**
	 * Is the current range was already checked?
	 */
	bool exist();


 public:
	/** 
	 * Constructor
	 */
	Anchor( std::string *s, 
			std::string *q,
			//KmerArray *sk,
			const KmerPosMap *sk,
			const KmerArray *qk,
			int m,
			int k, bool v );

	/** 
	 * Constructor
	 */
	Anchor( ReadId &r,
			PathId &p,
			std::string *s, 
			std::string *q,
			//KmerArray *sk,
			const KmerPosMap *sk,
			const KmerArray *qk,
			int m,
			int k, bool v );

	/** 
	 * Destructor
	 */
	~Anchor();

	/**
	 * Runner function.
	 */ 
	bool find();

	/**
	 * Setter: verbosity
	 */
	void   setVerbose(bool v) { verbose = v; }

	/**
	 * Getter: anchor search status
	 */
	bool   getStatus() { return found; }

	/**
	 * Getter: extended anchor length
	 */
	size_t getLength() { return length; }

	/**
	 * Getter: similarity score of extended anchor region
	 */
	double getScore()  { return score; }

	AnchorLog getLog() { return log; }
};


#endif
