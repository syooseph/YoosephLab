#include "sequence_build.h"
#include "assemble_extend_p.h"

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

#ifndef _CONTIG_REFINEMENT_P_
#define _CONTIG_REFINEMENT_P_

struct MerPosType {
  int cid;
  int pos;
  int fw_len;
  int re_len;
};

struct ContigOverlapType  {
  int cid;
  int count;
  int ref_begin, ref_end;
  int target_begin, target_end;
};

class ContigRefinementP  {
 public:
  ContigRefinementP();
  ~ContigRefinementP();
  void RefineContigsNoAln(
      SequenceBuild &seq_obj, int mer_len, 
      std::list<ContigType> &contigs, std::list<ContigType> &refined_contigs
  );
 private:
  bool init_contig_holder_;
  std::vector<ContigType> contig_holder_;
  void IndexContigs(
      int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& mer_contig 
  );
  void IncorporateContig(
      bool handle_index, int mer_len, 
      std::unordered_map<std::string, std::list<MerPosType> >& kmer_map,
      int contig_index, ContigType& current_contig, std::list<int>& incorporated_contigs
  );
  bool TryMergeOverlapedContigs(
      bool handle_index, int mer_len, ContigOverlapType& overlap_record, 
      ContigType& ref_contig, int& ref_index
  );
  bool FindConsensusTail(
      std::string& seq1, std::string& seq2, std::string& consensus
  );
  bool FindConsensusHead(
      std::string& seq1, std::string& seq2, std::string& consensus
  );
  std::string FixedWidthString(int len, int num);
};

#endif
