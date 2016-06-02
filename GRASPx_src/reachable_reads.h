#include "sequence_build.h"
#include "database_index.h"
#include "scoring_function.h"

#include <cstdlib>
#include <cstring>
#include <cmath>
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

#ifndef _REACHABLE_READS_
#define _REACHABLE_READS_

struct SeedType {
  std::string seed_seq;
  int q_pos;
};

struct ReadType {
  RIDType rid;
  int q_begin, q_end; // the approximate begin positions of the read in query
  int score;
  double aln_e_value;
  int q_pos, r_pos; // the position of the match in the query and the read, respectively
};

struct ReachabilityType {
  RIDType rid;
  int q_pos;  // the index of the end of the extension in query sequence
  int reached_step;
};

struct ReadConnectType  {
  RIDType rid;
  int q_begin, q_end;
  double aln_e_value;
  int overlap_len_prev, overlap_len_aft;
  bool connect_terminate;
};

class ReachableReads  {
 public:
  ReachableReads(
      std::string& reduc_map_file, std::string& seed_ext_file, 
      std::string& read_ext_file, std::string& hsm_file
  );
  ReachableReads();
  ~ReachableReads(void);
  void SelectSeeds(
      double seed_score_scale, ScoringFunction<int>& score_scheme, std::string& query, 
      std::map<int, std::list<SeedType> >& candidate_seeds
  );
  void CollectNeighbourReads(
      int num_steps,
      std::map<int, std::list<ReadPairType> >& seed_read_pairs, 
      std::unordered_map<RIDType, bool>& candidate_reads
  );
  void GetSeedReads(
      SequenceBuild& seq_obj,
      std::map<int, std::list<SeedType> >& candidate_seeds,
      std::map<int, std::list<ReadPairType> >& seed_reads
  );
  int GetOverlapLen(void);
  void ComputeHighScoreMatch(
      std::vector<char> alphabet, double def_high_frac,
      int mer_size, ScoringFunction<int>& score_scheme,
      ReducedAlphabet& reduc_alph
  );
  friend class GreedyAssembly;
  friend class AssembleExtend;
 protected:
  // variables
  int alph_id_;
  int seed_len_;
  int overlap_len_;
  std::unordered_map<std::string, std::unordered_map<std::string, bool> > reduc_alph_map_;
  std::unordered_map<std::string, std::list<ReadPairType> > seed_ext_;
  std::unordered_map<RIDType, std::list<OverlapType> > fw_read_ext_;
  std::unordered_map<RIDType, std::list<OverlapType> > re_read_ext_;
  std::unordered_map<std::string, std::list<std::string> > high_score_match_;
  // functions
  void RecruitFWReads(
      int num_steps, 
      RIDType seed_rid, std::unordered_map<RIDType, bool>& fw_reads
  );
  void RecruitREReads(
      int num_steps, 
      RIDType seed_rid, std::unordered_map<RIDType, bool>& re_reads
  );
  void InitContig(
      RIDType init_read,
      std::map<RIDType, bool>& read_table, 
      std::vector<std::deque<ReadConnectType> >& pre_contigs
  );
  void ConnectFW(
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_rext, 
      std::vector<std::deque<ReadConnectType> >& pre_contigs, int index,
      std::map<RIDType, bool>& read_table,
      std::unordered_map<int, std::list<int> >& fw_link,
      std::unordered_map<int, std::list<int> >& re_link
  );
  void ConnectRE(
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_rext, 
      std::vector<std::deque<ReadConnectType> >& pre_contigs, int index,
      std::map<RIDType, bool>& read_table,
      std::unordered_map<int, std::list<int> >& fw_link,
      std::unordered_map<int, std::list<int> >& re_link
  );
  void RefineLinks(
      std::unordered_map<RIDType, bool>& read_candidates,
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_link,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_link
  );
  bool CheckFWTerminate(
      std::list<OverlapType> ext_fw, 
      std::unordered_map<RIDType, std::list<OverlapType> >& re_rext
  );
  bool CheckRETerminate(
      std::list<OverlapType> ext_re, 
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext
  );
  void RefineSeedSameScore(
      int min_overlap, int query_len, 
      std::map<int, std::list<SeedType> >& candidate_seeds
  );
  void RefineSeedDiffScore(
      int min_overlap, 
      std::map<int, std::list<SeedType> >& candidate_seeds
  );
  bool IsSeqCompatible(
      SequenceBuild& seq_obj, int seed_len,
      RIDType fw_rid, POSType fw_pos,
      RIDType re_rid, POSType re_pos
  );
  void GetHighScoreSeq(
      int seed_len, int mer_len,
      std::string& seq, std::unordered_map<std::string, bool>& hs_seq
  );
  
};

#endif
