/** 
 * @file       kmer.h
 * @brief      K-mer and read type declaration
 * @date       2010-10-11 10:59:00 AM
 * @date       Modified on Tue 2013-12-17 06:52:49 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

/*
 * To do:
 * Test with a big kmer (>6)
 */

#ifndef __KMER_H__
#define __KMER_H__

#include <cstdlib>
//#include <boost/unordered_map.hpp>
#include <tr1/unordered_map>
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include "alpha.h"
#include "read.h"

typedef std::set<ReadId> ReadIdSet;
typedef std::vector<ReadId> ReadIdArray;
typedef std::list<ReadId> ReadIdList;
typedef std::string KmerType;

#if LONGKMER == 1
typedef unsigned long long KmerId;
#else
typedef unsigned KmerId;
#endif

typedef std::set<KmerId> KmerSet;
typedef std::vector<KmerId> KmerArray;
typedef std::tr1::unordered_map<KmerId, CoverageType> CoverageMap;
typedef std::tr1::unordered_map<KmerId, CoverageType*> LeftAAsCoverageMap;
typedef std::list<KmerId> KmerList;       
typedef std::tr1::unordered_map<KmerId, int> KmerScoreMap;
typedef std::tr1::unordered_map<KmerId, KmerId> PredKmerMap;
	
#endif
