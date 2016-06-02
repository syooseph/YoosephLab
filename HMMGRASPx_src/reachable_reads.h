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
  int num_extension;  // the number of extension that this seed contains
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
      SequenceBuild &seq_obj,
      std::string& reduc_map_file, std::string& seed_ext_file, 
      std::string& read_ext_file, std::string& hsm_file
  );
  ReachableReads();
  ~ReachableReads(void);
  /*
  void SelectSeeds(
      double seed_score_scale, ScoringFunction<int>& score_scheme, std::string& query, 
      std::map<int, std::list<SeedType> >& candidate_seeds
  );
  */
  /*
  void CollectNeighbourReads(
      int num_steps,
      std::map<int, std::list<ReadPairType> >& seed_read_pairs, 
      std::unordered_map<RIDType, bool>& candidate_reads
  );
  */
  void GetSeedReads(
      SequenceBuild& seq_obj,
      std::map<int, std::list<SeedType> >& candidate_seeds,
      std::map<int, std::list<ReadPairType> >& seed_reads
  );
  int GetOverlapLen(void);
  int GetSeedLen(void);
  /*
  void ComputeHighScoreMatch(
      std::vector<char> alphabet, double def_high_frac,
      int mer_size, ScoringFunction<int>& score_scheme,
      ReducedAlphabet& reduc_alph
  );
  */
  std::list<OverlapType>::iterator GetExtIterBeginFW(RIDType rid);
  std::list<OverlapType>::iterator GetExtIterBeginRE(RIDType rid);
  void ReadExtToVectors(void)  {
    for(auto it = fw_read_ext_.begin(); it != fw_read_ext_.end(); ++ it)  {
      fw_read_ext_vt_[it->first].resize(it->second.size());
      unsigned int id = 0;
      for(auto it_o = it->second.begin(); it_o != it->second.end(); ++ it_o)  {
        fw_read_ext_vt_[it->first][id] = *it_o;
        ++ id;
      }
    }
    for(auto it = re_read_ext_.begin(); it != re_read_ext_.end(); ++ it)  {
      re_read_ext_vt_[it->first].resize(it->second.size());
      unsigned int id = 0;
      for(auto it_o = it->second.begin(); it_o != it->second.end(); ++ it_o)  {
        re_read_ext_vt_[it->first][id] = *it_o;
        ++ id;
      }
    } 
    return;
  }
  friend class GreedyAssembly;
  friend class HMMProfile;
  friend class AssembleExtend;
  friend class AssembleExtendP;
 protected:
  // variables
  //int alph_id_;
  int seed_len_;
  int overlap_len_;
  std::unordered_map<std::string, std::unordered_map<std::string, bool> > reduc_alph_map_;
  std::unordered_map<std::string, std::list<ReadPairType> > seed_ext_;
  std::unordered_map<RIDType, std::list<OverlapType> > fw_read_ext_;
  std::unordered_map<RIDType, std::list<OverlapType> > re_read_ext_;
  std::unordered_map<RIDType, std::vector<OverlapType> > fw_read_ext_vt_;
  std::unordered_map<RIDType, std::vector<OverlapType> > re_read_ext_vt_;
  std::map<uint16_t, std::list<uint16_t> > high_score_match_;
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
  void RefineSeed(
      int min_overlap, int query_len, 
      std::map<int, std::list<SeedType> > &candidate_seeds
  );
  void RefineSeedSameScore(
      int min_overlap, int query_len, 
      std::map<int, std::list<SeedType> >& candidate_seeds
  );
  void RefineSeedDiffScoreHigh(
      int min_overlap, 
      std::map<int, std::list<SeedType> >& candidate_seeds
  );
  void RefineSeedDiffScoreLow(
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
