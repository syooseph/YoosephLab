#ifndef _STRING_GRAPH_H_
#define _STRING_GRAPH_H_

#include "bwt.h"
#include "bwt_search.h"
#include "frequency_table.h"
#include "kmer.h"
#include "align_batch.h"
#include "bio_alphabet.h"
#include "minimizer_sort.h"

#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <stack>
#include <unordered_map>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>



class STRVertexType {
 public:
  explicit STRVertexType(void) {rid_ = -1; traversed_ = false;}
  explicit STRVertexType(const int i) {rid_ = i; traversed_ = false;}
  inline void SetID(const int i) {if(rid_ == -1)  rid_ = i;}
  inline void SetTraversedTag(const bool i) {traversed_ = i;}
  inline bool IsTraversed(void) {return traversed_;}
  int rid_;   // the read ID for the current vertex
  int len_;   // the length of the read
  bool traversed_;  // tag to check whether the vertex has been visited
};

class STREdgeType {
 public:
  explicit STREdgeType(void) {len_ = -1; seq_ = "";}
  explicit STREdgeType(const int l) {len_ = l; seq_ = "";}
  inline void SetLen(const int i) {if(i > len_)  len_ = i;}
  inline void SetSeq(const std::string &s) {seq_ = s;}
  inline void SetTraversedTag(const bool i) {traversed_ = i;}
  inline bool IsTraversed(void) {return traversed_;}
  inline std::string GetSeq(void) {return seq_;}
  int sid_;   // the sequence ID for the current edge
  int len_;   // the length of the overlap between the two vertices 
              // (only make sense before condensing graphs)
  std::string seq_; // the sequence of the unitig
  bool traversed_;  // tag to check whether the edge has been traversed
  std::vector<int> path_info_;  // a vector contanins the path information of the edge
                                // with the following format: 
                                // read_ID:overlap_len:read_ID:overlap_len:read_ID:...:read_ID 
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, STRVertexType, STREdgeType> BoostSTRGraph;
typedef boost::graph_traits<BoostSTRGraph>::vertex_descriptor BoostSTRVertex;
typedef boost::graph_traits<BoostSTRGraph>::edge_descriptor BoostSTREdge;
typedef boost::graph_traits<BoostSTRGraph>::vertex_iterator BoostSTRVertexIter;
typedef boost::graph_traits<BoostSTRGraph>::edge_iterator BoostSTREdgeIter;

class ExtType {
 public:
  BWTIDX source_;  // the BWT index in the reverse BWT
  int rid_;        // the read ID of the source read
  // the set of irreducible reads' reverse BWT positions 
  // and the corresponding overlap lengths with the source
  std::vector<BWTIDX> ir_position_;
  std::vector<int>  ir_overlap_;
};

struct SeqAlnInfoType  {
  int q1, q2; // the interval for the query sequence
  int id;     // the ID of the sequence
  int score;  // the alignment score of the sequence
};

class StringGraph {
 public:
  explicit StringGraph(void) {initialized_ = false;}
  ~StringGraph(void) {}

  // import all information in the extension into the string graph
  // notice that the function can be called multiple times 
  // for differnt sets of extension information
  void ImportExtension(std::vector<ExtType> &extension);
  // batch execution of ComputeExtension
  void MultiComputeExtension(
      const int threads, const int min_overlap, 
      const int offset, const int n, char **seq, 
      BWT &bwt, BWT &rev_bwt, std::vector<ExtType> &extension
  );
  // computing the extension for a set of reads 
  // NOTE: bwt and rev_bwt are expected to be generated from the entire data set
  // output a set of "extension"s; each of them record the source and targets
  void ComputeExtension(
      const int min_overlap, const int offset, const int n, char **seq, 
      BWT &bwt, BWT &rev_bwt, std::vector<ExtType> &extension
  ); 
  // check if there exists some un-labeled vertex or un-labeld path
  // delete those un-labeld vertex or paths and clear the graph
  void CheckGraph(void);
  // check if self-cycle exists; if yes, remove them
  void CheckSelfCycle(void);
  // remove all orphant vertices (in_degree == 0 and out_degree == 0) in the graph 
  void RemoveOrphantVertices();
  // remove tips in the graph; return the number of vertices deleted
  int RemoveTipsBeforeCondense();
  // connect the unipaths and condense the graph
  void CondenseGraph(char** seq);
  // explore the condened graph and record all sequences that are longer than min_length
  void TraverseUnbranched(std::list<std::string> &paths, const int min_length = 60);
  // dump the graph to file for future re-loading
  void WriteGraph(BioAlphabet &alphabet, char** seq, const std::string &file_name);
  // loading the index file and construct the graph
  // "orphan_rid" and "orphan_reads" correspond to the sequences that cannot be overlapped
  // to any other reads; direct alignment of these sequences are sufficient
  void LoadGraph(
      const std::string &file_name, 
      std::vector<int> &orphan_rid, std::vector<std::string> &orphan_seq
  );
  // assign each edge sequence an ID, record that in the edge property "sid_"
  // also copy the sequences into "seqs" with an order consistent to the assigned sequence ID
  // the resulting sequence can be used to facilitate batch alignment
  // return the number of sequences
  int RecordEdgeSeqs(std::vector<std::string> &seqs);
  
  // computes the mapping between edge ID and edge pointer
  // "graph_edge[i]" stores the edge pointer which points to an edge with ID i
  void ComputeGraphEdgeMapping(std::vector<BoostSTREdge> &graph_edge);
  
  // use depth-first search to traverse the graph and extract high-scoring sequences
  // any edge whose alignment score is not found will receive a score 0; and all paths
  // with score-sum higher than the "cutoff" and with length close to "query_len" will 
  // be stored in "high_scoring_seqs"; "edge_seqs" are the set of edge sequences;
  // "graph_edge[i]" stores the edge index in the graph for the ith sequence in "edge_seqs";
  // "edge_ID" records the edge ID for the set of alignments, which corresponds to the scores
  // in "score"; in other words, the alignment between the query and the edge sequence with ID 
  // edge[i] has an alignment score of score[i]; "q_interval" are the intervals in the query
  // used to generate the alignment
  void GetHighScoringPaths(
      const int cutoff, const int query_len,
      std::vector<std::string> &edge_seqs,
      std::vector<BoostSTREdge> &graph_edge,
      std::vector<int> &edge_ID, std::vector<int> &score, 
      std::vector<std::pair<int, int> > &q_interval,
      const double sim_cutoff,
      std::vector<std::string> &high_scoring_seqs
  );
  
  void GetHighScoringPaths(
      const int cutoff, const int query_len,
      std::vector<std::string> &edge_seqs,
      std::vector<BoostSTREdge> &graph_edge,
      const int num_seqs, int *edge_ID, int *score, 
      int *q_interval, const double sim_cutoff,
      std::vector<std::string> &high_scoring_seqs
  );
  
  // given the "source_edge", progressively searching for incoming edges that are
  // found in hash table "aligned_seqs", which has the key as the target sequence ID
  // and the value as the alignment score. "expected_len" is the expected length for extension,
  // "seqs" is the set of target sequences, "expanded_seqs" stores the results 
  // and "visited_seqs" records the target sequences that have been traversed during expansion
  void ProgressiveExtendLeft(
      BoostSTREdge &source_edge, 
      const int expected_len, std::vector<std::string> &seqs,
      std::unordered_map<int, int> &aligned_seqs, 
      std::vector<std::pair<std::string, int> > &expanded_seqs,
      std::unordered_map<int, int> &visited_seqs
  );
  
  // refer to ProgressiveExtendLeft
  void ProgressiveExtendRight(
      BoostSTREdge &source_edge, 
      const int expected_len, std::vector<std::string> &seqs,
      std::unordered_map<int, int> &aligned_seqs, 
      std::vector<std::pair<std::string, int> > &expanded_seqs,
      std::unordered_map<int, int> &visited_seqs
  );
  
  // merging the left and right extended sequences
  void MergeExtensions(
      std::string &source_seq, 
      std::vector<std::pair<std::string, int> > &left_seqs, 
      std::vector<std::pair<std::string, int> > &right_seqs,
      std::vector<std::string> &merged_seqs
  );
  
  double EstSeqSimilarity(
      const std::string &seq1, const std::string &seq2, const bool is_left
  );
   
  /************ obsolete function *****************
  // clustering the sequences that correspond to branches from the same source
  // "identity" is the percetage of identity for sequence to be clustered together
  //void ClusteringSeqs(const int identity);
  *****************************/
  
  void PrintGraphSizes()  {
    std::cout << "The graph contains " << num_vertices(*p_graph_) << " vertices." << std::endl;
    std::cout << "The graph contains " << num_edges(*p_graph_) << " edges." << std::endl;
  }
  void Purge()  {
    if(initialized_) delete p_graph_;
  }
  
 private:
  bool initialized_;
  // the string graph
  BoostSTRGraph *p_graph_;
  // the hash table recording the rev_BWT index 
  // (which is expected to be corresponding to a full read with delimitor, i.e. $r$)
  // and the corresponding vertex label in the graph "p_graph_"
  std::unordered_map<BWTIDX, BoostSTRVertex> node_hash_;
  
  // condense the graph (for all reachable components from the source_edge)
  void Condense(char** seq, const BoostSTREdge source_edge);
};

#endif
