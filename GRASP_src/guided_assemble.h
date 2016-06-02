#include "extend_and_assemble.h"
#include "index_sample.h"
#include "gsa.h"
#include "sequence.h"
#include "assembly_graph.h"

#include "timer.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <cstdlib>
#include <set>

#ifndef _GUIDED_ASSEMBLE_H_
#define _GUIDED_ASSEMBLE_H_

typedef int AlignScoreType;
typedef std::pair<int, int> RegionType;
struct KmerInfoType {
  int kmer_begin;
  int kmer_coverage;
};

struct ReadMapType  {
  //QIDType query_ID;
  //std::vector<int> mapped_region;
  // the mapped region vector contains 4 entries, which are the start and end locations
  // of the query, and start and end locations of the target
  //AlignScoreType alignment_score;
  //AlignScoreType assembled_score;
  int phase_counter;
  int query_seq_ID;
  //BoostVertex max_vertex_source, max_vertex_target;
  //bool is_bridging_read;
  double best_evalue;
};

template<typename AlignScoreType>
struct AssembleParameterType  {
  unsigned int num_threads;
  unsigned int seed_len;
  enum Alphabet alphabet;
  enum MatrixName scoring_matrix;
  AlignScoreType gap_open;
  AlignScoreType gap_extension;
  unsigned int band_size;
  double evalue_cutoff;
  unsigned int n_back_check;
  unsigned int dropoff;
  bool write_certificate;
};

template<typename AlignScoreType>
struct CertificateType  {
  std::vector<int> alignment_regions;
  std::string alignment_query, alignment_target, alignment_symbol;
  AlignScoreType alignment_score;
  std::list<std::pair<RIDType, int> > mapped_locations;
};

class GuidedAssemble  {
 public:
  GuidedAssemble();
  GuidedAssemble(
      const std::string& in_working_directory,
      const std::string& in_result_directory,
      const std::string& in_query_sequence_file,
      const std::string& in_read_mfasta_file
  );
  GuidedAssemble(
      const std::string& in_working_directory,
      const std::string& in_result_directory,
      const std::string& in_query_sequence_file,
      const std::string& in_read_mfasta_file,
      const AssembleParameterType<AlignScoreType>& set_parameters
  );
  ~GuidedAssemble();
  void PrepareRunEnvironment();
  void AssembleAndAlign(std::unordered_map<RIDType, ReadMapType>& mapped_reads);
  void FetchSeeds(const std::string& sequence, std::set<std::string>& matched_seeds);
  void PreFilterAll(void);
  //void GetReadTags(std::vector<std::string>& tags);
  std::string GetFileStem(const std::string& path);  
  friend class ExtendFunctor;
  friend class CertificateFunctor;
  friend class RealignmentFunctor;
  
 protected:
  bool info_loaded_;
  unsigned int num_threads_;
  unsigned int num_running_threads_;
  LockType mutex_threads_;
  std::string working_directory_;
  std::string result_directory_;
  std::string query_sequence_file_;
  std::string read_mfasta_file_;
  std::vector<std::string> query_sequence_info_;
  std::vector<std::string> query_sequence_;
  unsigned int seed_len_;
  enum Alphabet alphabet_;
  AlignScoreType gap_open_;
  AlignScoreType gap_extension_;
  enum MatrixName scoring_matrix_;
  ScoringFunction<AlignScoreType> *score_scheme_;
  unsigned int right_band_, down_band_;
  double evalue_cutoff_;
  unsigned int n_back_check_;
  unsigned int dropoff_;
  unsigned int query_num_seqs_;
  //char** sample_tags_;
  char** sample_seqs_;
  double sample_size_;
  char** reversed_sample_seqs_;
  unsigned int sample_num_reads_;
  GSA *suffix_array_;
  GSA *reverse_suffix_array_;
  std::unordered_map<KmerType, std::list<PositionType> > kmer_positions_;
  bool write_certificate_;
  /*
  void SelectNonRedundantSeeds(
      IndexSample& index_caller,
      const std::string& query_sequence,
      const KmerInfoType& kmer, 
      const std::unordered_map<RIDType, std::list<ReadTerminateType> >& mapped_locations, 
      std::list<PositionType>& seed_positions
  );
  */
  void PrepareScoringFunction(void);
  void PrepareWorkSpace(void);
  void PrepareQuerySequence(void);
  void PrepareSampleSequence(void);
  void PrepareIndex(void);
  void PrepareSuffixArray(void);
  void PrepareReversedSuffixArray(void);
  void ExtendAllSeeds(void);
  void FinalRealignment(void);
  void DefineNonRepeatRegions(const std::string sequence, std::list<RegionType>& regions);
  void DefineKmersToSearch(
      const std::string sequence, 
      const std::list<RegionType>& regions, 
      std::list<KmerInfoType>& kmers
  );
  
  void ClearForwardSeq();
  void ClearReverseSeq();
  //void ClearTag();
  void ClearSuffixArrays();
  void TraceBack(
      const std::string& query_sequence,
      AssemblyGraph& graph, std::list<VertexPairType>& source_vertices,
      std::unordered_map<RIDType, ReadMapType>& mapped_reads
  );
  
  int ReconcileReads(
    const std::string& query_sequence,
    const AlignScoreType& cutoff,
    const AlignScoreType& opposite_max_score,
    const AlignScoreType& opposite_min_score,
    std::unordered_map<RIDType, ReadVertexType>& recorded_reads,
    std::unordered_map<RIDType, ReadMapType>& mapped_reads 
  );
  
  int ReconcileBridgingReads(
    const std::string& query_sequence,
    const AlignScoreType& cutoff,
    const AlignScoreType& left_max_score,
    const AlignScoreType& left_min_score,
    const AlignScoreType& right_max_score,
    const AlignScoreType& right_min_score,
    VertexPairType seed_pair,
    std::unordered_map<RIDType, ReadMapType>& mapped_reads 
  );
  
  unsigned int GetNumConcurrentThreads(void);
  void IncreaseNumThreads(void);
  void DecreaseNumThreads(void);
  void PartitionQuery(
      const std::string& sequence, 
      std::map<int, int>& clump_map
  );   
  void WriteAlignments(
      const int& query_ID, const std::string& query_sequence, 
      const std::list<CertificateType<AlignScoreType> >& certificates
  );
  void SelectRecruitedReads(
      const int& query_ID, const std::string& query_sequence, const std::list<CertificateType<AlignScoreType> >& certificates, 
      std::unordered_map<RIDType, ReadMapType>& mapped_reads
  );
  void WriteLog(void);
  void WriteRecruitedReads(const int& query_ID, std::unordered_map<RIDType, ReadMapType> mapped_reads);
  void WriteAlignment(const std::string& query_sequence, const CertificateType<AlignScoreType>& certificate);
  int NumGaps(const std::string& seq);
  AlignScoreType EvalHamming(const std::string seq_a, const std::string seq_b);
  void ProgressiveSearch(const std::string& seed, int step_forward, int branch_cutoff, std::set<RIDType>& reached_reads);
  void SearchSingleFw(std::set<std::string>& search_seed, std::set<RIDType>& containing_reads);
  void SearchSingleRe(std::set<std::string>& search_seed, std::set<RIDType>& containing_reads);
  void DefineNewSearchSeqFw(std::set<RIDType>& reads, std::set<std::string>& end_sequence);
  void DefineNewSearchSeqRe(std::set<RIDType>& reads, std::set<std::string>& end_sequence);
};

#endif
