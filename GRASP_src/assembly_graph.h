#include "index_sample.h"
#include "sfa.h"
#include "gsa.h"
#include "scoring_function.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/thread.hpp>
#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/smart_ptr.hpp>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <cassert>
#include <map>
#include <unordered_map>
#include <stack>

#ifndef _MAX_VALUE
#define _MAX_VALUE 30000
#endif

#ifndef _ASSEMBLY_GRAPH_H_
#define _ASSEMBLY_GRAPH_H_

typedef int AlignScoreType;
typedef int AlignmentPositionType;
typedef unsigned int QueryPositionType;
typedef unsigned int AlignmentIDType;
typedef enum{EXT_LEFT, EXT_RIGHT} DirectionType;

typedef boost::shared_mutex LockType;
typedef boost::shared_lock<LockType> SLock;
typedef boost::unique_lock<LockType> ULock;

// TODO: obsolete definition
/******************************************************
// a data structure denoting a bridging read, containing the read_ID and position of the seed,
// as well as the tree nodes that are connected by the bridging read
struct BridgingReadType {
  RIDType read_ID;
  POSType position; // for left extension, such position indicates the beginning of the suffix,
                    // for right extension, such position indicates the end of the prefix
  DirectionType extend_direction;
};
*****************************************************/

struct ReadAssignmentType {
  // remember we are recording the 3'end poistion of the read and see where it maps to
  RIDType read_ID;
  POSType position;
  AlignmentIDType alignment_ID;
  AlignmentPositionType alignment_position;
  DirectionType extend_direction;
  SfaType rank_in_sfa;
  int len_search_seed;
  int offset_to_start;  
  //int full_seq_len;
};

typedef ReadAssignmentType BridgingReadType;

struct ReadVertexType {
  AlignScoreType max_score, min_score;
  // expected to be "BoostVertex" or "boost::graph_traits<BoostGraph>::vertex_descriptor" type
  // "BoostGraph" is defined as "boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, VertexProperty, EdgeProperty>"
  void *max_vertex_source, *max_vertex_target, *min_vertex_source, *min_vertex_target; 
};

struct AssembledPathType  {
  int mapped_begin, mapped_end;
  std::string assembled_sequence;
  std::list<std::pair<RIDType, int> > mapped_locations;
  AlignScoreType align_score;
};

struct EdgeProperty {
  // the outgoing sequence that corresponds to the current extension
  //std::string outgoing_seq;
  unsigned int outgoing_seq_ID;
  unsigned int outgoing_seq_Pos;
  // the end of the query of the current alignment (the start of the query is defined within the class)
  QueryPositionType query_end;
  // the reads that have been included in this extension
  std::list<ReadAssignmentType> included_reads;
  // the alignment score achieved by this extension
  AlignScoreType current_score;
  AlignScoreType increased_score;
  // note that the extend_direction indicates when the path is discovered, i.e. from left or right extension
  // such information is used to traceback the assembled sequence (basically we have flip the sequence if the direction is left)
  // The direction is different from the graph direction, which always goes from 5' to 3'
  DirectionType extend_direction;
  std::list<BridgingReadType> candidate_bridging_reads;
};

struct VertexProperty {
  // the sequence that have been assembled so far
  //std::string assembled_seq;
  // the highest (global) alignment score achieved so far, given that the alignment edge ends at this vertex
  AlignmentPositionType query_seed_begin;
  AlignmentPositionType query_end;
  bool touched; // indicates whether the vertex is touched (with all its children put in the stack)
  bool visited; // indicates whether the vertex is visited (touched and max_score_decendent resolved)
  AlignScoreType max_score_below;  // the best (highest) score that can be achieved with the vertex as root
  AlignScoreType max_score_above; 
  int num_updates_below;  // number of updates received from its children
  void *max_vertex_above, *max_vertex_below;
  std::unordered_map<RIDType, ReadVertexType> reads_in_children;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, VertexProperty, EdgeProperty> BoostGraph;
typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;      
typedef boost::graph_traits<BoostGraph>::edge_descriptor BoostEdge;             
typedef boost::graph_traits<BoostGraph>::vertex_iterator BoostVertexIter; 
typedef boost::graph_traits<BoostGraph>::edge_iterator BoostEdgeIter;


class AssemblyGraph {
 public:
  AssemblyGraph();
  ~AssemblyGraph();
  // adding a vertex, returning the vertex descriptor
  BoostVertex AddVertex(const VertexProperty& data);
  // adding an edge between vertices u and v, returnnig the edge descriptor
  BoostEdge AddEdge(const BoostVertex& u, const BoostVertex& v, const EdgeProperty& data);
  VertexProperty AccessVertex(const BoostVertex& v);
  EdgeProperty AccessEdge(const BoostEdge& e);
  BoostEdge GetEdge(const BoostVertex& u, const BoostVertex& v);

  // update the information stored in a vertex
  void UpdateVertex(const BoostVertex& v, const VertexProperty& n);
  // update the information stored in an edge
  void UpdateEdge(const BoostEdge& e, const EdgeProperty& n);
  // get the number of vertices in the graph
  int GetNumVertices(void);
  // get the number of edges in the graph
  int GetNumEdges(void);
  // get the range of iterators
  std::pair<BoostVertexIter, BoostVertexIter> GetVertices(void);
  std::pair<BoostEdgeIter, BoostEdgeIter> GetEdges(void);
  // get the in/out edges
  void GetInEdges(const BoostVertex& v, std::list<BoostEdge>& ie_list);
  void GetOutEdges(const BoostVertex& v, std::list<BoostEdge>& oe_list);
  // get the in/out/total degree
  int GetInDegree(const BoostVertex& v);
  int GetOutDegree(const BoostVertex& v);
  int GetDegree(const BoostVertex& v);
  // get the source and target of an edge
  BoostVertex GetSource(const BoostEdge& e);
  BoostVertex GetTarget(const BoostEdge& e);
  // delete edge or vertex
  void RemoveEdge(const BoostEdge& e);
  void RemoveVertex(const BoostVertex& v);
  // take the source node of a path and spell all paths that follow it
  void SpellMinimumPaths(
      GSA& suffix_array, GSA& reverse_suffix_array,
      const BoostVertex& source, std::set<RIDType> reads_to_visit,
      std::list<AssembledPathType>& assembled_seqs,
      AssembledPathType& best_seq
  );
  void SpellAllPaths(
    GSA& suffix_array, GSA& reverse_suffix_array, 
    const BoostVertex& source, std::list<AssembledPathType>& assembled_seqs
  );
  void TraverseGraphUnder(
      const BoostVertex& source, 
      std::unordered_map<RIDType, ReadVertexType>& recorded_reads
  );
  
 protected:
  LockType mutex_;  // the mutex for the entire graph, will be locked only when adding edges/vertices
  BoostVertex source_;  // the vertex that corresponds to the seed extension, beginning of the graph traversal
  BoostGraph graph_;    // the graph that we work on 
  void UpdateAbove(
      BoostVertex& source, BoostVertex& target, 
      AlignScoreType path_score_max
  );
  void UpdateBelow(
      BoostVertex& source, BoostVertex& target, 
      AlignScoreType path_score_max
  );
  void ExtendEdge(EdgeProperty& edge_current, EdgeProperty& edge_extend, EdgeProperty& edge_next);
  /*
  void TrackDownBestPath(
      const BoostVertex& pivot,
      EdgeProperty& e_current,
      EdgeProperty& e_complete
  );
  */
  void MergeReadsWithParent(BoostVertex& v_parent, BoostVertex& v_child);
  void AddEdgeReadsToParent(
      BoostVertex& v_parent, BoostVertex& v_child,
      AlignScoreType max_score
  );
};

#endif
