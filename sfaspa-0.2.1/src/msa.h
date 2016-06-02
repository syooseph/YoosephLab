/** 
 * @file       msa.h
 * @brief      Multiple sequence alignment
 * @date       Modified on Tue 2013-12-17 06:48:36 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __MSA_H__
#define __MSA_H__

#include "palign.h"


typedef int SeqPos;     
typedef int RefPos;
typedef unsigned Count;
typedef std::pair<SeqPos, RefPos> PosPair;  
typedef std::list<PosPair> PosPairList;   
typedef std::map<ReadId, PosPairList> ReadPosPairListMap;
typedef std::tr1::unordered_map<ReadId, Count> ReadCountsMap;
typedef std::map<RefPos, ReadCountsMap> PosReadCountsMap;
typedef std::tr1::unordered_map<ReadId, Count> ReadCountsMap;
typedef std::map<RefPos, Count> PosCountsMap;     
typedef std::vector<RefPos> RefPosArray;
typedef std::vector<std::string> Sequences;

/**
 * \brief Multiple sequence alignment
 */
class MSA
{
 private:
	std::string       pivot;
	ReadIdArray       rids;
	Sequences         nreads;
	RefPosArray       nstart;
	char              **reads;
	ReadPlacementList *places;


 public:
	MSA();
	MSA( PathAligner *a, char **r);//, Param &p );
	~MSA();
	void run();
	void print(std::ostream &out, int csize);
	
	std::string getPivot()   { return pivot; }
	Sequences   getSequences() { return nreads; }
	void setPivot(std::string &p) { pivot = p; }

 private:
	void insertGapsToReads();
	void insertGapsToPivot();

	ReadPosPairListMap loadInsertions();
	PosReadCountsMap countReadPositions(ReadPosPairListMap &pmap);
	PosCountsMap getMaxCounts(PosReadCountsMap &prmap);
	void insertGapsToPivot(PosCountsMap &cmap);
	void stretch();
	void adjustInsertions(ReadPosPairListMap &rpmap,
							   PosReadCountsMap &prmap,
							   PosCountsMap &cmap );
	void equalize();

};

#endif
