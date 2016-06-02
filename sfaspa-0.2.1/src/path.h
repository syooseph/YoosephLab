/** 
 * @file       path.h
 * @brief      Graph path declaration.
 * @date       Wed 2011-06-08 08:04:21 PM
 * @date       Modified on Tue 2013-12-17 06:41:48 PM
 * @author     Youngik Yang
 * @version    0.001
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __SPAPATH_H__
#define __SPAPATH_H__

#include <iostream>
#include <list>
#include <set>
#include <boost/unordered_map.hpp>
#include "kmer.h"
#include "biostr.h"
#include "semiglobal.h"
#include "param.h"
#include "gsa.h"

typedef int PathId;
typedef std::tr1::unordered_map<PathId, PathId> PathIdMap;
typedef std::vector<BoundType> BoundArray;
typedef std::list<BoundArray> SuffixBounds;
typedef std::list<size_t> TraceSizeList;
typedef std::pair<PathId, AlignSummary> PathAlignPair;
typedef std::list<PathAlignPair> PathAlignPairList; 
typedef std::tr1::unordered_map<PathId, PathAlignPairList> PathBins;

typedef std::tr1::unordered_map<PathId, std::string> PathSequences;
typedef std::vector<PathId> PathIdArray;

/* struct UsedInfo */
/* { */
/* 	unsigned gsa; */
/* 	unsigned len; */
/* 	unsigned num; */
/* 	UsedInfo( unsigned g, unsigned l, unsigned n ) { */
/* 		gsa = g; len = l; num = n; */
/* 	} */
/* }; */

struct UsedBound
{
	SfaType  min;
	SfaType  max;
	//unsigned gsa;
	LcpType  gsa;
	//unsigned len;

	//UsedBound(SfaType l, SfaType m, unsigned g, unsigned n) 
	UsedBound(SfaType l, SfaType m, LcpType g)
	{
		/* min = l; max = m; gsa = g; len = n; */
		min = l; max = m; gsa = g;
	}
	
	bool operator<( const UsedBound& other) const
	{
		if ( min == other.min ) {
			if ( max == other.max ) {
				/* if ( gsa == other.gsa )  */
				/* 	return len < other.len; */
				return gsa < other.gsa;				
			}
			return max < other.max;
		}
		return min < other.min;
	}

	void print() {
		//printf("l:%d, r:%d, gsa:%u, len:%u", min, max, gsa, len);
		printf("l:%d, r:%d, gsa:%u", min, max, gsa);
	}
};

//typedef std::multimap<BoundType, UsedInfo> BoundUsedMap;
//typedef std::pair<BoundType, UsedInfo> BoundUsedInfo;
//typedef std::set<BoundType, UsedInfo> BoundUsedMap;

//typedef std::pair<UsedBound, unsigned> BoundUsedInfo;
typedef std::pair<UsedBound, unsigned> UsedBoundInfo;
typedef std::list<UsedBoundInfo> UsedBoundInfoList;
typedef std::list<UsedBoundInfoList> UsedBoundInfoLists;

typedef std::map<UsedBound, size_t> UsedBoundMap;

const PathId NOT_PATH = -1;
const PathId BAD_READ = -2; 

/**
 * Path search stop criteria.
 */
enum STOPCONDITION { PATHEND,     /**< Path ends by either sink node or stop condon */
					 EXTENDFAIL,  /**< Extend fails because all nodes are weakly supported */
					 WEAKMAX,     /**< Weak depth of maximum neighbor node */
					 BADSTART,    /**< Stop codon at source */
					 CYCLICPATH,  /**< Cycle */
					 REDUNDANT,   /**< Redundant */
};


/** 
 * \brief Path object identified from a graph
 */
struct GraphPath
{
	KmerId        seed;    ///< seed kmer 
    PathType      path;    ///< a set of nodes
	int           spos;    ///< seed kmer position in a path
	KmerArray     kmers;   ///< kmers that corresponds to nodes
    int           lstop;   ///< left extension termination codition
	int           rstop;   ///< right extension termination condition
	SuffixBounds  bounds;  ///< suffix array search results for each node extension
	TraceSizeList traces;  ///< n-back trace size for each node extension
	UsedBoundInfoLists used_begs;
	UsedBoundInfoLists used_ends;

	int getLStop() { return lstop; }
	int getRStop() { return rstop; }
	std::string toString(int k);
	void dump(std::fstream &out);
	void load(std::fstream &in, size_t);
	void clear();
	void release();
	void trim( int n, bool l );
};


/**
 * \brief Mismatch position index
 */
struct Mismatch {
	ReadId read;    ///< read id
	int    qry_pos; ///< read pos
	int    ref_pos; ///< path pos
	
	Mismatch(){}

	Mismatch( ReadId rid, int sp, int rp ) 
	{
		read = rid; qry_pos = sp; ref_pos = rp;
	}

	Mismatch( const Mismatch &source )
	{
		read = source.read; qry_pos = source.qry_pos; ref_pos = source.ref_pos;
	}

	Mismatch& operator= ( const Mismatch &source )
	{
		read = source.read; qry_pos = source.qry_pos; ref_pos = source.ref_pos;
		return *this;
	}

	void dump(std::ostream &out) 
	{
		out.write((char*)&read, sizeof(ReadId));
		out.write((char*)&qry_pos, sizeof(int));
		out.write((char*)&ref_pos, sizeof(int));
	}

	void load(std::istream &in) 
	{
		in.read((char*)&read, sizeof(ReadId));
		in.read((char*)&qry_pos, sizeof(int));
		in.read((char*)&ref_pos, sizeof(int));
	}
};

typedef std::pair<int, int> IntPair;


#endif

