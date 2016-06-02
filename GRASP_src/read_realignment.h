#include "scoring_function.h"
#include "seq_align.h"
#include "sfa.h"
#include "gsa.h"
#include "index_sample.h"
#include "extend_and_assemble.h"
#include "assembly_graph.h"

#include <iostream>
#include <assert.h>
#include <vector>
#include <list>
#include <string>
#include <cstdlib>
#include <time.h>
#include <map>
#include <unordered_map>

#ifndef _READ_REALIGNMENT_H_
#define _READ_REALIGNMENT_H_

/*
  this class does the alignment between the recruited reads and the query sequence
  TODO: the output should be a bed-like file that tells where the reads are aligned
*/

typedef int AlignScoreType;
typedef unsigned int QIDType;

struct ReadMapType  {
  //QIDType query_ID;
  //std::vector<int> mapped_region;
  // the mapped region vector contains 4 entries, which are the start and end locations
  // of the query, and start and end locations of the target
  //AlignScoreType alignment_score;
  //AlignScoreType assembled_score;
  int phase_counter;
  BoostVertex max_vertex_source, max_vertex_target;
  bool is_bridging_read;
  double best_evalue;
};

class ReadRealignment {
 public:
  ReadRealignment();
  ReadRealignment(
    const std::string& query_seq,
    const int& query_begin, const int& query_end,
    const std::string& assembled_seq,
    ScoringFunction<AlignScoreType>& score_scheme,
    const unsigned int& down_band, const unsigned int& right_band,
    std::list<ReadAssignmentType>& reads_in_path
  ); 
  ~ReadRealignment();
  void Align();
 protected:
  QIDType query_ID_;
  std::string query_seq_;
  std::string assembled_seq_;
  int query_begin_, query_end_;
  ScoringFunction<AlignScoreType>* score_scheme_;
  unsigned int right_band_, down_band_;
  std::list<ReadAssignmentType> path_reads_;
  
  void ReAlignSequences(
      AlignScoreType& opt_score, double& opt_bit_score,
      std::unordered_map<int, int>& nuc_match, 
      std::string& aligned_query, std::string& aligned_assembled, std::string& aligned_symbols
  );
  /*
  void FormatOutput(
      AlignScoreType& opt_score, double& opt_bit_score,
      std::unordered_map<int, int>& nuc_match, 
      std::string& aligned_query, std::string& aligned_assembled, std::string& aligned_symbols
  );
  */
};

#endif
