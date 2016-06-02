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
};

struct MapInfoType  {
  int contig_id;  // the ID of the aligned contig
  double e_value; // the e_value of the aligned contig
  double error_score; // defined by the ((num_errors + 1) / aligned_portion)
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
};

inline int CountUngappedError(std::string& s1, std::string& s2) {
  int ne = 0;
  if(s1.length() != s1.length())  {
    std::cout << "Error: ReMap::CountUnGappedError:  ungapped matching applied on sequences with different lengths" << std::endl;
    return 9999;
  }
  for(unsigned int i = 0; i < s1.length(); ++ i) {
    if(s1[i] != s2[i])  {
      ++ ne;
    }
  }
  return ne;
}

inline int Get5PSeq(
    int num_buffer, std::string& s1, int s1_index, std::string& s2, int s2_index,
    std::string& s1_substr, std::string& s2_substr
) {
  // gets the 5' end sequence
  if(s1_index <= s2_index)  {
    s1_substr = s1.substr(0, s1_index);
    int s2_len = s2_index >= s1_index + num_buffer ? s1_index + num_buffer : s2_index;
    s2_substr = s2.substr(s2_index - s2_len, s2_len); 
    return s1_index;
  } else  {
    s2_substr = s2.substr(0, s2_index);
    int s1_len = s1_index >= s2_index + num_buffer ? s2_index + num_buffer : s1_index;
    s1_substr = s1.substr(s1_index - s1_len, s1_len); 
    return s2_index;
  }
  return 0;
}

inline int Get3PSeq(
    int num_buffer, std::string& s1, int s1_index, std::string& s2, int s2_index,
    std::string& s1_substr, std::string& s2_substr
) {
  // gets the 3' end sequence
  int s1_suf_len = s1.length() - s1_index;
  int s2_suf_len = s2.length() - s2_index;
  if(s1_suf_len <= s2_suf_len)  {
    s1_substr = s1.substr(s1_index, s1_suf_len);
    int s2_len = s2_suf_len >= s1_suf_len + num_buffer ? s1_suf_len + num_buffer : s2_suf_len;
    s2_substr = s2.substr(s2_index, s2_len); 
    return s1_suf_len;
  } else  {
    s2_substr = s2.substr(s2_index, s2_suf_len);
    int s1_len = s1_suf_len >= s2_suf_len + num_buffer ? s2_suf_len + num_buffer : s1_suf_len;
    s1_substr = s1.substr(s1_index, s1_len); 
    return s2_suf_len;
  }
  return 0;
}

inline double CalErrorScore(int num_errors, double aligned_portion)  {
  // 1 is added as pseudo-count
  return (double) (num_errors + 1) / aligned_portion;
}

inline bool CompMapQuality(bool is_qual, MapInfoType &m1, MapInfoType &m2) {
  // is_qual set to 1 means quality first (i.e. error_score is compared with priority)
  // otherwise e_value is compared with priority
  // the return value indicates: 1: m2 is better than m1; 0: m2 is not better than m1
  if(is_qual && 
      (m2.error_score < m1.error_score || 
      (m2.error_score == m1.error_score && m2.e_value < m1.e_value))
  )  {
    return true;
  } else if(!is_qual &&
      (m2.e_value < m1.e_value ||
      (m2.e_value == m1.e_value && m2.error_score < m1.error_score))
  )  {
    return true;
  }
  return false;
}


#endif
