#include "reachable_reads.h"
#include "database_index.h"
#include "sequence_build.h"
#include "hmm_profile.h"
#include "seq_align_extend_p.h"

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
#include <bitset>

#ifndef _ASSEMBLE_EXTEND_P_
#define _ASSEMBLE_EXTEND_P_

struct AlignInfoType  {
  bool is_fw;                     // indicate whether the extension is FW or RE
  int score;                      // score achieved in current extension
  int opt_score;                  // optimal score has been achieved so far
  int opt_contig_len;             // length of the contig that corresponds to optimal score
  RIDType opt_rid;                // the rid corresponds to the optimal terminus
  int mbegin;                     // begin(FW)/end(RE) of the suffix(FW)/prefix(RE) of the model
  int m_filled, s_filled;         // regions that have been filled in the DP-matrix
  RIDType last_rid;               // the read corresponds to the suffix of the contig
  int ext_len;                    // length of the extension (length of suffix of last_rid)
  int ext_id;                     // pointer to the next extension of the last_rid
  int acc_len;                    // accumulated length
  int num_extensions;             // number of extension have been performed so far
  bool recorded;                  // indicate whether the node has been recorded
                                  // it is a soft-mask, recorded alignment cannot be used
                                  // as a leaf node for back-tracking
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
  RIDType fw_rid, re_rid;
};

class AssembleExtendP  {
 public:
  AssembleExtendP(void) {;}
  ~AssembleExtendP(void) {;}
  void AssembleAllContigs(
      HMMProfile &profile_obj, SequenceBuild& seq_obj,
      ReachableReads& link_obj, int band_size,
      std::map<int, std::list<SeedMatchType> >& seed_reads,
      bool dup_seeds, int assembly_depth, double p_value_cutoff,
      std::list<ContigType>& contigs
  );
  // performs progressive extension until a branching path is found
  void ProgressiveExtension(
      SequenceBuild& seq_obj,
      ReachableReads& link_obj,
      int n,
      std::list<ContigType>& contigs
  );
 protected:
  inline void UpdateRedundancyRecord(
      uint8_t *alned_reads, 
      RIDType read_id
  );
  inline bool IsReadRedundant(
      uint8_t *alned_reads, 
      RIDType read_id
  );
  void InitFW(
      SequenceBuild& seq_obj, ReadPairType& seed_read,
      ReachableReads& link_obj,
      SeqAlignExtendP &align_obj, AlignInfoType *fw_aln, int &fw_aln_index
  );
  void InitRE(
      SequenceBuild& seq_obj, ReadPairType& seed_read,
      ReachableReads& link_obj,
      SeqAlignExtendP &align_obj, AlignInfoType *re_aln, int &re_aln_index
  );
  void ExtendSeedReadFW(
      SequenceBuild& seq_obj, int assembly_depth,
      uint8_t *alned_reads,
      std::unordered_map<RIDType, std::vector<OverlapType> >& fw_rext,
      SeqAlignExtendP &fw_align_obj,
      bool dup_seeds, int score_cutoff,
      AlignInfoType *fw_aln, int &fw_aln_index, 
      std::list<ContigType>& fw_contigs
  );
  void ExtendSeedReadRE(
      SequenceBuild& seq_obj, int assembly_depth,
      uint8_t *alned_reads,
      std::unordered_map<RIDType, std::vector<OverlapType> >& re_rext,
      SeqAlignExtendP &re_align_obj,
      bool dup_seeds, int score_cutoff,
      AlignInfoType *re_aln, int &re_aln_index,
      std::list<ContigType>& re_contigs
  );
  void MergeContigs(
      HMMProfile &q_profile, ReadPairType &seed_pair,
      SequenceBuild& seq_obj, ReachableReads& link_obj,
      int score_cutoff, int seed_score,
      std::list<ContigType>& fw_contigs, std::list<ContigType>& re_contigs,
      std::list<ContigType>& merged_contigs
  );
  void FilterSeedReads(
      SequenceBuild& seq_obj,
      ReachableReads& link_obj,
      std::list<ContigType>& contigs,
      std::unordered_map<std::string, bool> &filtered_seeds
  );
};

inline void AssembleExtendP::UpdateRedundancyRecord(
    uint8_t *alned_reads,
    RIDType read_id
) {
  int q = read_id / 8;
  int r = read_id % 8;
  //uint8_t c = 1;
  alned_reads[q] |= 1 << (8 - r - 1);
  //std::cout << "!!! " << (unsigned int) read_id << std::endl;
  return;
}

inline bool AssembleExtendP::IsReadRedundant(
    uint8_t *alned_reads, 
    RIDType read_id
) {
  int q = read_id / 8;
  int r = read_id % 8;
  uint8_t k = alned_reads[q] & (1 << (8 - r - 1));
  //std::cout << "??? " << (unsigned int) read_id << "  " << (bool) k <<  std::endl;
  return (bool) k;
}


#endif
