#ifndef _DEBRUIJN_GRAPH_H_
#define _DEBRUIJN_GRAPH_H_

#include <string>
#include <list>
#include <stack>
#include <unordered_map>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>

#include "frequency_table.h"
#include "kmer.h"
//#include "hash_bucket.h"

class VertexType{
 public:
  explicit VertexType() :
    multiplicity_(0),
    color_(0),
    traversed_(false)
  {}
  explicit VertexType(std::string &seq) :
    seq_(seq),
    multiplicity_(1),
    color_(0),
    traversed_(false)
  {}
  
  inline VertexType& operator++ (){++ this->multiplicity_; return *this;}
  VertexType& operator= (const VertexType &vertex)  {
    this->seq_ = vertex.seq_;
    this->multiplicity_ = vertex.multiplicity_;
    this->color_ = vertex.color_;
    this->traversed_ = vertex.traversed_;
    return *this;
  }
  std::string seq_;
  int multiplicity_;
  int color_;
  bool traversed_;
};

class EdgeType{
 public:
  explicit EdgeType() :
    multiplicity_(1)
  {}
  explicit EdgeType(int multi) : 
    multiplicity_(multi)
  {}
  inline EdgeType& operator++ (){++ this->multiplicity_; return *this;}
  inline EdgeType& operator-- (){-- this->multiplicity_; return *this;}
  inline EdgeType& operator= (const EdgeType &edge) {
      this->seq_ = edge.seq_;
      this->multiplicity_ = edge.multiplicity_; 
      return *this;
  }
  std::string seq_;
  int multiplicity_;
};

// obsolete definition
/*
class KmerInfo{
 public:
  // constructor
  explicit KmerInfo() {}
  explicit KmerInfo(const Kmer &kmer_obj, const int count)  {
    kmer_obj_ = kmer_obj;
    count_ = count;
  }
  // Destructor
  ~KmerInfo(){}
  
  Kmer kmer_obj_;
  int count_;
  inline bool operator==(const KmerInfo &kinfo)  {
    return (this->kmer_obj_ == kinfo.kmer_obj_);
  }
};
*/

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, VertexType, EdgeType> BoostGraph;
typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;
typedef boost::graph_traits<BoostGraph>::edge_descriptor BoostEdge;
typedef boost::graph_traits<BoostGraph>::vertex_iterator BoostVertexIter;
typedef boost::graph_traits<BoostGraph>::edge_iterator BoostEdgeIter;


struct ShortBridgeType  {
  BoostVertex source, target;
  int multiplicity, len;
};

class DeBruijnGraph {
 public:
  explicit DeBruijnGraph(KmerUnitcoder &encoder, const int kmer_size) : 
    initialized_(false), 
    encoder_(encoder),
    kmer_size_(kmer_size)
    {}
  ~DeBruijnGraph();
  
  // check if all unit k-mer in seq have higher frequency than freq_cutoff
  bool PassFreqFilter(
      const std::string &seq, 
      const FrequencyTable<KmerUnitType> &unit_freq, const int freq_cutoff
  );
  // construct the de bruijn graph from the sequence
  void Construct(
      const int num_seq, char ** const seq, std::vector<bool> &read_mark,
      const FrequencyTable<KmerUnitType> &unit_freq, const int freq_cutoff
  );
  // delete k-mers that are less frequent than the threshold and remove corresponding edges
  void RefineGraph(const int freq_cutoff);
  void RemoveEdges(std::list<BoostEdge> &to_delete);
  void RemoveOrphantVertices();
  // compute a set of low-multiplicity edges to delete 
  // to transform the graph into a minimum spanning tree
  // min_cov is the minimum coverage for a specific edge
  void FormMinSpanningTree(const int min_cov);
  // setting the color for all vertices in vertex_list to the input color
  void ReplaceColor(std::list<BoostVertex> &vertex_list, const int color);
  void TraverseUnbranched(std::list<std::string> &paths, const int min_length = 60);
  
  // remove head/tail tips. for all paths shorter than the length, remove them
  void RemoveTipsHead(const int length = 10);
  void RemoveTipsTail(const int length = 10);
  
  // the function delete all edges whose coverage is less than 
  // fraction*(coverage of one of its neighouring vertices)
  void BreakBadCovBranch(const float fraction = 0.05);
  
  void CondenseGraph();
  void Condense(const BoostEdge source_edge);
  
  void BridgeThreading(const int max_len);
  
  void PrintGraphSizes()  {
    std::cout << "The graph contains " << num_vertices(*p_graph_) << " vertices." << std::endl;
    std::cout << "The graph contains " << num_edges(*p_graph_) << " edges." << std::endl;
  }
  
  void PrintBranchingVertices();
  
  // obsolete function
  //void TraverseUnbranched2(std::list<std::string> &paths);
  
  /*
  void Construct2(
      const int num_seq, char ** const seq, 
      const FrequencyTable<KmerUnitType> &unit_freq, const int freq_cutoff
  );
  */
  
 protected:
  KmerUnitcoder encoder_;
  int kmer_size_;
  bool initialized_;
  BoostGraph *p_graph_;
  std::unordered_map<std::string, BoostVertex> *p_info_map_;
  // obsolete variable
  //HashBucket<KmerUnitType, KmerInfo> *p_info_map2_;
};

#endif
