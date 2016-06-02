#ifndef __RELATION_H__
#define __RELATION_H__

#include <queue>
#include "kmer.h"
#include "debruijn.h"
#include "graph.h"


//using namespace boost;

//typedef void* NodeType;
typedef std::list<KmerId> KmerList;
typedef Vertex NodeType;
typedef std::tr1::unordered_map<NodeType, KmerId> VertexToKmerMap;
//typedef std::tr1::unordered_map<NodeType, KmerList> VertexToKmerListMap;

typedef EdgeWeight WeightType;
typedef Vertex NodeType;
//typedef std::vector<NodeType> PathType;
typedef std::vector<NodeType> NodeArray;
typedef std::list<NodeType> PathType;
typedef PathType Path;
//typedef std::list<NodeType> PathType;
//typedef std::vector<PathType> Paths;
//typedef std::vector<PathType> PathArray;
//typedef std::list<PathType> Paths;
typedef std::list<NodeList> PathList;
typedef PathList Paths;

//typedef std::tr1::unordered_map<NodeType, KmerId> VertexMap;
//typedef std::tr1::unordered_map<NodeType, std::vector<KmerId> > VertexToKmerMap;
typedef std::tr1::unordered_map<KmerId, NodeType> KmerToVertexMap;
//typedef std::tr1::unordered_map<KmerId, KmerList> KmerListMap;
//typedef std::tr1::unordered_map<NodeType, NodeType> NodeToNodeMap;

// removed edge because of indel
//typedef std::tr1::unordered_map<NodeType, std::pair<NodeType, EdgeWeight> > LostEdgeMap;

//typedef std::set<PathType> PathSet;
typedef std::queue<PathType> PathQueue;

//typedef std::pair<NodeType, NodeType> NodePair;
//typedef std::set<NodePair> NodePairSet;

//typedef std::pair<KmerId, KmerId> KmerPair;
//typedef std::set<KmerId> KmerIdPairSet;

//typedef std::vector<ReadId> ReadIdArray;

//typedef std::list<ReadIdList> ReadsList;

#endif
