//=============================================================================
// $Id: kmer.h 2010-10-11 10:59:00 AM yyang@jcvi.org $
//=============================================================================

#ifndef __KMER_H__
#define __KMER_H__

#include <cstdlib>
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

//const static size_t BRANCH = 23;

typedef std::string KmerType;
//typedef unsigned KmerId;

#if LONGKMER == 1
typedef unsigned long long KmerId;
#else
typedef unsigned KmerId;
#endif


/*
 * To do:
 * unsigned long long when kmer size > 6
 */

//typedef std::string TagType;

typedef std::vector<KmerId> KmerArray;
//typedef std::set<ReadId> ReadArray;
typedef std::tr1::unordered_map<KmerId, CoverageType> CoverageMap;
typedef std::tr1::unordered_map<KmerId, CoverageType*> LeftAAsCoverageMap;
//typedef std::tr1::unordered_map<KmerId, CoverageType*> PrevAAsCoverageMap;
//typedef std::map<TagType, ReadId> ReadIdMap;
//typedef std::tr1::unordered_map<KmerId, ReadArray> KmerToReadMap;
//typedef std::tr1::unordered_map<KmerId, KmerArray> KmerMemberMap;


#endif
