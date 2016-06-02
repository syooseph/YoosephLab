#ifndef __TRAVERSAL_H__
#define __TRAVERSAL_H__

#include <iostream>
#include <queue>
#include <vector>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "debruijn.h"
#include "utils.h"
#include "kmer.h"
#include "alpha.h"
//#include "biostr.h"
//#include "relation.h"
#include "timer.h"
#include "param.h"
#include "file.h"

/* Shortcut of graph properties */
typedef boost::graph_traits<DeBruijnGraph>::vertex_descriptor  Vertex;       ///< Graph vertex type
typedef boost::graph_traits<DeBruijnGraph>::edge_descriptor    Edge;         ///< Graph edge type     
//typedef boost::graph_traits<DeBruijnGraph>::adjacency_iterator Adj_iter;     ///< Adjacent 
typedef boost::graph_traits<DeBruijnGraph>::vertex_iterator    Vertex_iter;  ///< Vertex iterator
typedef boost::graph_traits<DeBruijnGraph>::edge_iterator      Edge_iter;    ///< Edge iterator
typedef boost::graph_traits<DeBruijnGraph>::out_edge_iterator  OutEdge_iter; ///< Out-edge iterator

/* Shortcut for vertex containers */
typedef std::vector<Vertex> NodeArray; ///< A vector of vertices
typedef std::set<Vertex>    NodeSet;   ///< A set of vertices
typedef std::list<Vertex>   NodeList;  ///< A list of vertices
 
/* Hash table shortcuts */
typedef std::tr1::unordered_map<Vertex, int>        InDegreeMap;  ///< Map: no. of indegrees of a vertex
typedef std::tr1::unordered_map<Vertex, NodeArray > AdjacencyMap; ///< Map: neighboring nodes of a vertex


typedef std::list<KmerId> KmerList;
typedef Vertex NodeType;
typedef std::tr1::unordered_map<NodeType, KmerId> VertexToKmerMap;

typedef EdgeWeight WeightType;
typedef Vertex NodeType;

typedef std::vector<NodeType> NodeArray;
typedef std::list<NodeType> PathType;
typedef PathType Path;

typedef std::list<NodeList> PathList;
typedef PathList Paths;

typedef std::tr1::unordered_map<KmerId, NodeType> KmerToVertexMap;

typedef std::queue<PathType> PathQueue;


//

//#include "container.h"

enum searchdirection { LEFT, RIGHT };
enum nodetype { SOURCE, SINK, INTERNAL };
enum searchmode { STRICT, AGGRESSIVE };


/* enum { WHITE, GRAY, BLACK }; */
enum { IN, OUT }; // indegree, outdegree

namespace graph
{
	void build( DeBruijnGraph &graph,
				CoverageMap &kmer_coverage,
				VertexToKmerMap &vertex_map,
				int kmer_size, 
				std::string &graph_input,
				bool verbose );

	void buildGraph( DeBruijnGraph &graph, 
					 CoverageMap &kmer_coverage,
					 KmerToVertexMap &vertex_map,
					 int kmer_size,
					 std::string &graph_input,
					 bool verbose );
	
	VertexToKmerMap revertMap( KmerToVertexMap & );

	/* void updateTopology(DeBruijnGraph &,  */
	/* 					KmerToVertexMap &,  */
	/* 					KmerId,  */
	/* 					unsigned,  */
	/* 					Param & ); */
	void updateTopology( DeBruijnGraph &graph, 
						 KmerToVertexMap &vertex_map,
						 KmerId kmer_id,
						 CoverageType left[],
						 int kmer_size );
	
	void dropAncestor( AdjacencyMap &ancestors, 
					   Vertex key,
					   Vertex val);

	void dropAncestors( AdjacencyMap &ancestors, 
						PathType &minor,
						PathType &major );
		
	
	NodeSet trimGraph( DeBruijnGraph &graph,
					   CoverageMap &kmer_coverage,
					   VertexToKmerMap &vertex_map,
					   int min_depth,
					   bool verbose);
	
	NodeList trimVertices( DeBruijnGraph &graph,
					   CoverageMap &kmer_coverage,
					   VertexToKmerMap &vertex_map,
					   size_t min_coverage );
	
	void trimEdges( DeBruijnGraph &graph,
					size_t min_coverage );

	NodeList dropIslands( DeBruijnGraph &graph,
						  VertexToKmerMap &vertex_map);
	
	
	bool formCycle( NodeArray &path, 
					NodeType  target );
	
	bool formCycle( NodeList &path, 
					NodeType target );
	bool formCycle( NodeArray &path, 
					NodeArray &nodes );

	bool hasNode( PathType &path,
				  NodeType target );
		
	
	void BFS( PathQueue &queue,
			  DeBruijnGraph &graph,
			  int max_depth );

	void reverseBFS( PathQueue &queue,
					 DeBruijnGraph &graph,
					 AdjacencyMap &ancestors,
					 int max_depth );

/* 	void DFS( PathQueue &queue, */
/* 			  DeBruijnGraph &graph, */
/* 			  int max_depth ); */

/* 	void reverseDFS( PathQueue &queue, */
/* 					 DeBruijnGraph &graph, */
/* 					 AdjacencyMap &ancestors, */
/* 					 int max_depth ); */





	
	InDegreeMap indegrees(DeBruijnGraph &g);

	NodeArray predecessors(DeBruijnGraph &g, 
								  Vertex v);

	AdjacencyMap predecessorMap(DeBruijnGraph &g);

	NodeArray successors(DeBruijnGraph &g, Vertex v);

	AdjacencyMap successorMap(DeBruijnGraph &g);

	NodeList getSources(DeBruijnGraph &g);

	NodeList getSinks(DeBruijnGraph &g);

	NodeList getJunctions(DeBruijnGraph &g);

	bool hasVertex(DeBruijnGraph &graph, Vertex vertex);

	bool isSource(Vertex v,
				  DeBruijnGraph &g, 
				  AdjacencyMap &ancestors);
	bool isSink(Vertex v, 
				DeBruijnGraph &g, 
				AdjacencyMap &ancestors) ;

	void graphSummary( DeBruijnGraph & graph, std::ostream &out );

	void __extractIndegrees(DeBruijnGraph &graph, InDegreeMap &in_degree, double indeg[]);

	void __extractOutdegrees(DeBruijnGraph &graph, double outdeg[]);
	
	void degreeStats( DeBruijnGraph &graph );
	
}

#endif
