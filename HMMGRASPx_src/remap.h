#include "gsa.h"
#include "sequence.h"
#include "sequence_build.h"
#include "assemble_extend.h"

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

#ifndef _REMAP_
#define _REMAP_

struct MapReadType  {
  int rid;  // the read ID in the database
  int contig_id;  // the contig ID from the assembled ones
  int num_errors;
  double aligned_portion;
};

struct MapInfoType  {
  int contig_id;  // the ID of the aligned contig
  double e_value;
  int num_errors;
  double aligned_portion;
};

class ReMap {
 public:
  ReMap();
  ~ReMap();
  void LoadContigs(
      std::string& contig_file, std::list<ContigType>& loaded_contigs
  );
  void MapToContig(
      SequenceBuild& seq_obj, std::list<ContigType>& contigs, 
      int seed_len, double portion, int num_errors, bool gap_allowed,
      std::list<MapReadType>& mapped_reads
  );
  std::string GetContigInfo(int contig_id);
 protected:
  int num_seqs_;
  bool contigs_loaded_;
  char **header_;
  char **sequence_;
  inline int CountUngappedError(std::string& s1, std::string& s2);
  inline int Get5PSeq(
      int num_buffer, std::string& s1, int s1_index, std::string& s2, int s2_index,
      std::string& s1_substr, std::string& s2_substr
  );
  inline int Get3PSeq(
      int num_buffer, std::string& s1, int s1_index, std::string& s2, int s2_index,
      std::string& s1_substr, std::string& s2_substr
  ); 
  inline double CalErrorScore(int num_errors, double aligned_portion);
  inline bool CompMapQuality(bool is_qual, MapInfoType &m1, MapInfoType &m2);
  void MapUngapped(
      std::string &r_seq, std::string &t_seq,   // the input sequences
      int r_anchor, int t_anchor,               // the beginning of the anchor
      int seed_len,                             // the length of the seed match
      int num_sub_errors,                       // the number of allowed substitution errors
      int &detected_errors,                     // Outputs: number of detected errors
      int &mapped_begin, int &mapped_end        // Outputs: the indexes in r_seq
  ); 
};


inline double ReMap::CalErrorScore(int num_errors, double aligned_portion)  {
  // 1 is added as pseudo-count
  return (double) (num_errors + 1) / aligned_portion;
}

inline bool ReMap::CompMapQuality(bool is_qual, MapInfoType &m1, MapInfoType &m2) {
  // is_qual set to 1 means quality first (i.e. error_score is compared with priority)
  // otherwise e_value is compared with priority
  // the return value indicates: 1: m2 is better than m1; 0: m2 is not better than m1
  double m1_error_score = CalErrorScore(m1.num_errors, m1.aligned_portion);
  double m2_error_score = CalErrorScore(m2.num_errors, m2.aligned_portion);
  if(is_qual && 
      (m2_error_score < m1_error_score || 
      (m2_error_score == m1_error_score && m2.e_value < m1.e_value))
  )  {
    return true;
  } else if(!is_qual &&
      (m2.e_value < m1.e_value ||
      (m2.e_value == m1.e_value && m2_error_score < m1_error_score))
  )  {
    return true;
  }
  return false;
}


#endif
