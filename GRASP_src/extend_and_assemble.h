#include "scoring_function.h"
#include "seq_align.h"
#include "seq_align_extend.h"
#include "index_sample.h"
#include "tree_basic.h"
#include "assembly_graph.h"
#include "gsa.h"

#include <cstdlib>
#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <unordered_map>
#include <tuple>
#include <algorithm>

#ifndef _MAX_VALUE
#define _MAX_VALUE 30000
#endif

#define DEBUG 0

#ifndef _EXTEND_AND_ASSEMBLE_H_
#define _EXTEND_AND_ASSEMBLE_H_

typedef int AlignScoreType;
typedef unsigned int PriorityType;  

struct ReadTerminateType  {
  // records where the read finally being aligned with and whats the current alignment
  // score of the entire assembled sequence
  //RIDType read_ID;
  AlignmentPositionType alignment_position;
  AlignScoreType assembled_score_left, assembled_score_right;
  // the edge and vertex where the read goes through and goes into
  //BoostEdge end_edge_left, end_edge_right;
  //BoostVertex end_vertex_left, end_vertex_right;
  //SfaType rank_in_sfa_left, rank_in_sfa_right;
  int phase_counter_left, phase_counter_right;
  bool left_extended, right_extended;  
};

struct OutgoingPathType {
  SfaType key_sequence_position;          // the representative sequence position in the suffix array
  unsigned int len_search_seed;           // the length of the search seed sequence that should be trimed
      // one can spells out the sequence with sfa::getSuffix_explicit
  std::list<ReadAssignmentType> included_reads;      // a list of included reads for this outgoing path
  std::list<BridgingReadType> candidate_bridging_reads; // a list of candidate bridging reads at this outgoing path
  DirectionType extend_direction;
};


// a data structure to be stored as a node in the assembled tree
struct AlignmentNodeType  {
  std::string outgoing_seq;                 // the outgoing sequence between the current node and its parental node
  AlignScoreType max_score;               // the (accumulative) alignment score up to this node 
  std::list<ReadAssignmentType> recruited_reads;
                                          // a list of recruited reads
};


// the structure that stores the seeding vertices (the ones that connect the left and right extensions)
struct SeedVertexType {
  BoostVertex seed_vertex_left, seed_vertex_right;
  bool has_initialized_left, has_initialized_right;
  AlignScoreType seed_match_score;
  std::string seed_sequence;
};

struct VertexPairType {
  BoostVertex vertex_left, vertex_right;
  std::string gap_sequence;
  //std::string seed_sequence;
  int seed_in_gap_seq;
  std::list<ReadAssignmentType> bridging_reads;
  int seed_start;
  AlignScoreType seed_match_score, left_score, right_score;
  bool has_initialized_left, has_initialized_right;
};

template<typename SCORETYPE>
struct AlignmentInfoType  {
  AlignmentIDType alignment_ID;             // the ID of the alignment
  unsigned int current_query_length;
  //std::string current_query_seq;               // the first sequence of the alignment
  std::string current_assembled_seq;          // the second sequence of the alignment
  std::string seed_extend_seq;                // the seed sequence to be used to extend the alignment (with n_back_check concaternated)
  DirectionType extend_direction;             // the direction of the alignment (EXT_LEFT or EXT_RIGHT)
  bool is_seed_complete;                    // indicate whether the seed contains complete extension sequence and n_back_check
  unsigned int seed_coverage;               // the coverage of the seed sequence
  EdgeScoreBundle<SCORETYPE> edge_score;  // the resulting edge scores of the alignment
  SCORETYPE max_score;                    // the maximum alignment score of the alignment
  double best_bit_score;
  PriorityType priority_score;                 // the prioirty of the alignment for extend 
                                           // (a combination of seed_coverage and max_score)
  AlignScoreType prev_score;              // the alignment score of the previous extension
  BoostVertex source_vertex;              // the vertex in the assembly graph for the current alignment
  bool has_source_initialized;              // tag indicating whether the source vertex has been initialized
  // temporary buffer for the extensions that have not been submitted. the sizes of the following lists should be 
  // the same, and the size of each of them (i.e. non_increasing_count) should not exceeds a threshold, which will
  // be passed in as a parameter
  unsigned int non_increasing_count;
  unsigned int non_inclusion_count;
  std::list<VertexProperty> vertex_holder;   // the list holds the vertices that have not been submitted due to 
                                        // the non-increasing alignment score.
  std::list<EdgeProperty> edge_holder;    // the list holds the edges that associate with the vertices in vertex holder
                                        // the sizes of these two lists should be the same
  //std::list<SeedVertexType>::iterator seed_iterator;   // the iterator of the list where the seed of the current alignment path resides 
};


class ExtendAndAssemble {
 public:
  ExtendAndAssemble();
  ExtendAndAssemble(
      // the following three arguments locate the seed in the query sequence
      const unsigned int& in_num_threads,
      const std::string& in_query_seq, 
      const unsigned int in_start,
      const unsigned int in_seed_length,
      // a list of matched seed positions in the target/read set
      const std::string& in_seed_target,
      // the forward and reverse suffix arrays
      GSA& in_suffix_array,
      GSA& in_reverse_suffix_array,
      // the length of overlap that a read has to have with the assembled sequence
      const unsigned int in_n_back_check,
      const unsigned int in_dropoff,
      // the parameter controlling how close two reads map to the query sequence should be considered as one mapping
      std::map<int, int>& in_clump_map,
      // the scoring function for alignment
      ScoringFunction<AlignScoreType>* in_score_scheme,
      // parameters for the banded sequence alignment
      const unsigned int in_right_band, const unsigned int in_down_band,
      // the cutoff for a path to be recorded
      const double in_bit_score_cutoff,
      // outputs: the read IDs and their roughly aligned locations
      std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* in_reads_to_record,
      std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* in_assigned_reads,
      LockType* in_mutex_assigned_reads,
      // outputs: the assembly graph that encodes all the assembled sequences
      AssemblyGraph& in_assembly_graph,
      std::list<VertexPairType>& in_source_vertices
  );
  ~ExtendAndAssemble();
  //void LoadSeedSequences(std::list<std::string>& loaded_seeds);  // This functions
      // loads all seed sequences from the sample (short reads), which will then 
      // be used as extension seeds. Note that the current seed in the query may 
      // not have exact match with the seed in the sample due to the reduced
      // alphabet that has been used.
  void ExtendToBothDirections();  // This function 
      // initialize the alignments that extend the seed to both left and right. 
  /*****************************obsolete definition**************************************************
  //void ExtendToRight(const std::list<std::string>& loaded_seeds); // This function 
      // initialize the alignments that extend the seed to the right.
  //void ExtendToLeft(const std::list<std::string>& loaded_seeds);  // This function
      // initialize the alignments that extend the seed to the left.
  *******************************************************************************/
  void RightExtensionWorker(void* alignment_to_extend, void* out_path, BoostEdge& ext_e, BoostVertex& ext_v);
  void LeftExtensionWorker(void* alignment_to_extend, void* out_path, BoostEdge& ext_e, BoostVertex& ext_v);
  // increase/decrease the thread counter
  void IncreaseThreadCounter(void);
  void DecreaseThreadCounter(void);
  
 protected:
  // mutex for shared variables
  LockType* mutex_assigned_reads_;
  // variable definition
  AlignmentIDType alignment_counter_;
  //std::unordered_map<AlignmentIDType, FinishedAlignmentType> traversed_alignments_;
  unsigned int num_threads_;
  unsigned int num_concurrent_threads_;
  const std::string* query_seq_;
  unsigned int start_;
  unsigned int seed_length_;
  std::string seed_target_;
  unsigned int right_band_; // intersection between the band the the first row
  unsigned int down_band_;  // intersection between the band and the first column
  GSA* suffix_array_;
  GSA* reverse_suffix_array_;
  ScoringFunction<AlignScoreType>* score_scheme_;
  unsigned int n_back_check_;
  unsigned int dropoff_;
  std::map<int, int>* clump_map_;
  double bit_score_cutoff_;
  //std::unordered_map<unsigned int, bool>* covered_region_;
  std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* reads_to_record_;
  std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* assigned_reads_;
  AssemblyGraph* assembly_graph_;
  SeedVertexType seed_vertices_;
  // a list holding all high-score
  // alignments that are qualified to be extend further. The extended alignment
  // will be removed from the list, and newly qualified alignments will be added
  // to the end of the list.
  std::list<AlignmentInfoType<AlignScoreType> > high_score_alignments_right_;
  std::list<AlignmentInfoType<AlignScoreType> > high_score_alignments_left_; 
  // a map of left-most vertices that serve as the source node for traversal
  std::list<VertexPairType>* source_vertices_; 
  // function definition
  int GetKmerCoverageLeft(std::string kmer);
  int GetKmerCoverageRight(std::string kmer);
  void ExtendToLeftPhase(void);
  void ExtendToRightPhase(void);    // This function takes one alignment from the
      // high_score_alignments_ and extend it. It will also add new alignments
      // that are qualified for extension into high_score_alignments_.
  void ComputeAlignmentPriority(AlignmentInfoType<AlignScoreType>* a);
  void InsertToAlignmentHolderRight(const AlignmentInfoType<AlignScoreType>* a);
  void InsertToAlignmentHolderLeft(const AlignmentInfoType<AlignScoreType>* a);
  void ComputeOutgoingPathsLeft(
      AlignmentInfoType<AlignScoreType>& source_alignment,
      std::list<OutgoingPathType>& out_paths
  );
  void ComputeOutgoingPathsRight(
      AlignmentInfoType<AlignScoreType>& source_alignment,
      std::list<OutgoingPathType>& out_paths
  );
  std::string GetQueryExtendSeqLeft(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  // This function computes the possible outgoing paths to the right direction
  std::string GetQueryExtendSeqRight(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  void SelectValidSuffixLeft(
      const std::string& assembled_seq, // the assembled sequence which should be the prefix of the valid suffix
      const std::string& seed_seq,      // the seeding sequence that defines the range in the prefix/suffix array
      BoundType seed_range,             // the range returned by searching the suffix array using seed_seq
      std::unordered_map<unsigned int, bool>& valid_indicator // a list indicating whether the corresponding suffix
          // is valid. The first field correspond the suffix's position in the suffix array
  );
  void SelectValidSuffixRight(
      const std::string& assembled_seq, // the assembled sequence which should be the prefix of the valid suffix
      const std::string& seed_seq,      // the seeding sequence that defines the range in the prefix/suffix array
      BoundType seed_range,             // the range returned by searching the suffix array using seed_seq
      std::unordered_map<unsigned int, bool>& valid_indicator // a list indicating whether the corresponding suffix
          // is valid. The first field correspond the suffix's position in the suffix array
  );
  void IntersectBoundRanges(
      GSA* ref_array,             // the reference array whose elements are being tested
      const BoundType ref_array_bound,  // the range in the reference array
      GSA* rev_ref_array,         // the reverse of the ref_array to check the prefix
      const BoundType rev_ref_array_bound,  // the range in the reverse reference array
      int seed_len,     // the length of seed that is used to search both suffix arrays
      int n_back_check, // size of the read that overlaps with the assembled sequence, in additional to the seed   
      int max_rev_length, // the maximum length of the reverse reference sequence, longer than such value will be deemed invalid
      std::unordered_map<unsigned int, bool>& valid_indicator
  );
  /*********************************************************************************
  void AddRecruitedReadsLeft(
      const OutgoingPathType& current_path, 
      const AlignmentInfoType<AlignScoreType>& source_alignment
  );
  void AddRecruitedReadsRight(
      const OutgoingPathType& current_path, 
      const AlignmentInfoType<AlignScoreType>& source_alignment
  );
  *********************************************************************************/
  //void SubmitRecruitedReads(const AlignmentInfoType<AlignScoreType>& current_alignment);
  // an updated version of SubmitRecruitedReads (taking the reads from the recorded edges)
  int SubmitEdgeReads(
      AlignmentInfoType<AlignScoreType>& current_alignment, 
      BoostEdge& current_edge, BoostVertex& current_vertex
  );
  void ClearRecruitedReads(AlignmentInfoType<AlignScoreType>& current_alignment);
  void RecordRecruitedReadsLeft(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      const std::list<SfaType>& reads_positions,
      std::list<ReadAssignmentType>& reads_assignments
  );
  void RecordCandidateReadsLeft(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      const std::list<SfaType>& candidate_positions,
      std::list<BridgingReadType>& candidate_bridging_reads
  );
  void RecordRecruitedReadsRight(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      const std::list<SfaType>& reads_positions,
      std::list<ReadAssignmentType>& reads_assignments
  );
  void RecordCandidateReadsRight(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      const std::list<SfaType>& candidate_positions,
      std::list<BridgingReadType>& candidate_bridging_reads
  );
  bool HasAssignedLeft(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  bool HasAssignedRight(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  // check whether the read has been assigned (recorded into the "assigned_reads_" hash table)
  // the function LookupAssignedVertex() also record the ending vertex of the read, if it has been assigned
  bool HasAssigned(const RIDType read_ID, const AlignmentPositionType position);
  bool LookupAssignedVertex(const ReadAssignmentType& lookup_read);
  // computing the Right/Left boundaries in the query for the read that corresponds
  // to the extend_seq (in case of right extension the extend_seq_len is the length
  // of the suffix, and in case of left extension the extend_seq_len is the length
  // of the prefix)
  int ComputeBridgingReadsRightBound(const SfaType& rank, const DirectionType& direction);
  int ComputeQueryRightBoundRight(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  int ComputeQueryLeftBoundRight(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  int ComputeQueryRightBoundLeft(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  int ComputeQueryLeftBoundLeft(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      SfaType extend_position
  );
  void FillExtendSeq(
      const std::string& out_path,
      AlignmentInfoType<AlignScoreType>& current_alignment
  );
  // double checks the bridging reads to see if their sequences are compatible with 
  // the currently assembled sequences (the left and the right tree)
  //void SubmitBridgingReads(TreePairType* tree_pair);
  
  // the function based on the bridging reads to link the extended paths together, and recruit those valid bridging reads
  // note that this function is expected to alter the structure of the graph
  void ProcessBridgingReads(const SeedVertexType& sv);
  void ProcessIsolatedReads(const SeedVertexType& sv);
  //void RecordBridgingRead(const BridgingReadType& to_record);
  
  // function for doing alignment at each extension phase
  void DoExtendAlignment(
      // the source alignment information
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      // the query extension and assembled sequences (not the cummulative, but only the part gets extended)
      const std::string& query_seq_phase, const std::string& out_seq_phase,
      // a list of included reads that are associated with the out_seq_phase
      const unsigned int& read_ID, const unsigned int& position, 
      const std::list<ReadAssignmentType>& included_reads,
      // the candidate bridging reads
      const std::list<BridgingReadType>& candidate_reads, 
      // the alignment information recorder that serves as the output media
      AlignmentInfoType<AlignScoreType>& current_alignment
  ); 
  // the function tries to look up whether the current "extend_path" has reached a vertex that have been previously extended
  // if not, return false, if yes return true and record such vertex into "end_vertex"
  bool LookupExtendedPaths(
      const OutgoingPathType& extend_path, 
      int& num_new_reads, BoostEdge& end_edge, BoostVertex& end_vertex
  );
  // this function is used to enforce that the extend sequence in the merged path is actually consistent to ensure that
  // a single traversal of the path will spell out the correct assembled sequence
  /*
  bool ConsolidateMergedPaths(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      OutgoingPathType& extend_path, const BoostEdge& end_edge, const BoostVertex& end_vertex, 
      EdgeProperty& consolidated_edge
  );
  */
  int ComputeOutSeqDifference(const ReadAssignmentType& extend_read, const ReadAssignmentType& edge_read);
  /*
  void BuildNewEdge(
      const AlignmentInfoType<AlignScoreType>& source_alignment,
      const OutgoingPathType& extend_path, const EdgeProperty& edge_content, 
      const DirectionType& extend_direction, const DirectionType& read_direction, int distance,
      EdgeProperty& new_edge
  );
  */
  bool UpdateAssignedReads(const RIDType& read_ID, const ReadTerminateType& mapped_info);
  bool IsSeqEndMatch(std::string seqA, std::string seqB);
  
  // function to initialize the seeds in the assembly graph
  void InitSeedPair(
      AlignmentInfoType<AlignScoreType>& left_alignment, 
      AlignmentInfoType<AlignScoreType>& right_alignment
  );
  
  // handling the merging of the paths
  bool DetectPathMerge(
    AlignmentInfoType<AlignScoreType>& source_alignment,
    OutgoingPathType& out_path
  );
  // get how many threads have been issued
  unsigned int GetNumConcurrentThreads(void);
  
  
};

#endif
