#include "reachable_reads.h"
#include "database_index.h"
#include "scoring_function.h"
#include "seq_align_extend.h"
#include "sequence_build.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <deque>
#include <vector>
#include <tuple>
#include <string>
#include <list>

#ifndef _ALIGN_EXTEND_
#define _ALIGN_EXTEND_

struct AlignInfoType  {
  std::string contig;
  EdgeScoreBundle<int> edge_score;
  int opt_score, opt_bit_score;
  int score, bit_score;
  int prev_bit_score;
  int opt_contig_len;
  bool is_fw;
  int q_begin, q_end;
  int opt_q_begin, opt_q_end;
  RIDType last_rid;
  int num_extensions;
};

struct AlignmentPrintType {
  std::string seq1, seq2, symbol, header;
  std::unordered_map<int, int> nuc_match;
};

struct ContigType {
  std::string sequence;
  int score; 
  double bit_score;
  double e_value;
  int q_begin, q_end;
  bool valid;
  AlignmentPrintType al_print;
};

struct RangeType  {
  int begin, end;
};

struct ReadPosType  {
  RIDType rid;
  int begin, end;
};

class AssembleExtend  {
 public:
  AssembleExtend(void);
  ~AssembleExtend(void);
  void AssembleAllContigs(
      std::string& query, SequenceBuild& seq_obj,
      ReachableReads& link_obj, 
      ScoringFunction<int>& score_obj, int band_size,
      std::map<int, std::list<ReadPairType> >& seed_reads,
      int assembly_depth,
      double e_value_cutoff,
      std::list<ContigType>& contigs
  );
  std::string EncodeHeader(std::string& barcode, ContigType& contig);
  bool DecodeHeader(std::string& header, std::string& barcode, ContigType& contig);
 protected:
  inline bool IsReadRedundant(
      std::unordered_map<RIDType, bool>& alned_reads, 
      RIDType read_id
  );
  inline void UpdateRedundancyRecord(
      std::unordered_map<RIDType, bool>& alned_reads, 
      RIDType read_id
  );
  void ExtendSeedReadFW(
      std::string& query, SequenceBuild& seq_obj,
      ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
      std::unordered_map<RIDType, bool>& alned_reads,
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext, 
      std::stack<AlignInfoType>& fw_aln, std::list<ContigType>& fw_contigs
  );
  void ExtendSeedReadRE(
      std::string& query, SequenceBuild& seq_obj,
      ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
      std::unordered_map<RIDType, bool>& alned_reads,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_rext, 
      std::stack<AlignInfoType>& re_aln, std::list<ContigType>& re_contigs
  );
  void MergeContigs(
      std::string& query, SequenceBuild& seq_obj, ReachableReads& link_obj,
      ScoringFunction<int>& score_obj, 
      double e_value_cutoff, int seed_score,
      std::list<ContigType>& fw_contigs, std::list<ContigType>& re_contigs,
      std::list<ContigType>& merged_contigs
  );
  double ComputeJaccardIndex(int i1, int j1, int i2, int j2);
  std::string FetchQuerySeqFW(std::string& query, int q_begin, int len);
  std::string FetchQuerySeqRE(std::string& query, int q_end, int len);
  void InitFW(
      std::string& query, SequenceBuild& seq_obj, ReadPairType& seed_read, 
      ScoringFunction<int>& score_obj, int band_size,
      ReachableReads& link_obj,
      std::stack<AlignInfoType>& fw_aln
  );
  void InitRE(
      std::string& query, SequenceBuild& seq_obj, ReadPairType& seed_read, 
      ScoringFunction<int>& score_obj, int band_size,
      ReachableReads& link_obj,
      std::stack<AlignInfoType>& re_aln
  );
  void ExtendFWPhase(
      std::string& query, SequenceBuild& seq_obj,
      ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
      std::unordered_map<RIDType, bool>& alned_reads,
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext, 
      AlignInfoType& src_aln, std::list<AlignInfoType>& emitted_alns,
      std::list<AlignInfoType>& terminus
  );
  void ExtendREPhase(
      std::string& query, SequenceBuild& seq_obj,
      ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
      std::unordered_map<RIDType, bool>& alned_reads,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_rext, 
      AlignInfoType& src_aln, std::list<AlignInfoType>& emitted_alns,
      std::list<AlignInfoType>& terminus
  );
  bool DoExtendAlnFW(
      std::string& query, RIDType ext_rid,
      std::string& seed_query, std::string& seed_contig,
      std::string& extd_query, std::string& extd_contig, 
      ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
      AlignInfoType& src_aln, AlignInfoType& emit_aln
  );  
  bool DoExtendAlnRE(
      std::string& query, RIDType ext_rid,
      std::string& seed_query, std::string& seed_contig,
      std::string& extd_query, std::string& extd_contig, 
      ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
      AlignInfoType& src_aln, AlignInfoType& emit_aln
  );
  
};

inline void AssembleExtend::UpdateRedundancyRecord(
    std::unordered_map<RIDType, bool>& alned_reads, 
    RIDType read_id
) {
  alned_reads[read_id] = true;
  return;
}

inline bool AssembleExtend::IsReadRedundant(
    std::unordered_map<RIDType, bool>& alned_reads, 
    RIDType read_id
) {
  auto it = alned_reads.find(read_id);
  if(it == alned_reads.end())  {
    // if the read is not contained then there is no redundancy
    return false;
  } else  {
    return true;
  }
  return true;
}

#endif
