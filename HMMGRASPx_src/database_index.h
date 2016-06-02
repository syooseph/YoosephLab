#include "sequence_build.h"
#include "reduced_alphabet.h"
#include "scoring_function.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>
#include <list>
#include <bitset>

#ifndef _DATABASE_INDEX_
#define _DATABASE_INDEX_

struct OverlapType  {
  RIDType rid;
  POSType len;  // the length of the overlap
};

struct ReadPairType {
  RIDType rid;
  POSType r_pos; 
  int q_pos; 
  int score;  // matching score for the seeds in the query and target
  double aln_e_value;
};

struct MinReadPairType  {
  RIDType rid;
  POSType pos;
};

class DatabaseIndex {
 public:
  DatabaseIndex(int in_alph_id, int in_seed_len, int in_overlap_len);
  DatabaseIndex();
  ~DatabaseIndex();
  uint16_t Convert3mer(const char *s);
  void Interpret3merIndex(const uint16_t c, char *s);
  void BuildSeedmerMap(SequenceBuild& seq_obj); // map from seed-mer to sequence position
  // the key is the reduced alphabet string, the value is
  // a list of RID (read ID) - POS (position) pairs of the 
  // corresponding strings in the original alphabet
  void CreateReducedMap(
      std::unordered_map<std::string, std::list<std::string> >& reduc_alph_map
  );
  void DumpReducedMap(
      std::unordered_map<std::string, std::list<std::string> >& reduc_alph_map,
      std::string& out_file
  );
  void LoadReducedMap(
      std::string& in_file, 
      std::unordered_map<std::string, std::unordered_map<std::string, bool> >& reduc_alph_map
  );
  // the key is a pair of RID-POS that specify a given seed k-mer
  // the value is a list of RIDs which can be used to extend the seed
  /*
  void CreateSeedExt(
      int min_seed_coverage, SequenceBuild& seq_obj, SequenceBuild& rev_seq_obj, 
      std::unordered_map<std::string, std::list<ReadPairType> >& ext_seed
  );
  */
  void DumpSeedExt(
      std::list<MinReadPairType> &ext_seed,
      std::string& out_file
  );
  void LoadSeedExt(
      std::string& in_file,
      SequenceBuild& seq_obj,
      std::unordered_map<std::string, std::list<ReadPairType> >& ext_seed
  );
  // the key is a RID, and the value is a list of RIDs that have
  // significant overlap with the key read
  void CreateReadExt(
      int min_ext_coverage,
      SequenceBuild& seq_obj, 
      SequenceBuild& rev_seq_obj,
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read
  );
  void DumpReadExt(
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read,
      std::string& out_file
  );
  void LoadReadExt(
      std::string& in_file,
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read
  );
  void CreateHighScore3mer(std::map<uint16_t, std::list<uint16_t> > &high_score_match);
  void DumpHighScore3mer(std::map<uint16_t, std::list<uint16_t> > &high_score_match, std::string& out_file);
  void LoadHighScore3mer(
      std::string& in_file,
      std::map<uint16_t, std::list<uint16_t> > &high_score_match
  );
  int GetSeedLen(void);
  int GetOverlapLen(void);
  int GetAlphID(void);
  void GetSeedExtNew(
    SequenceBuild &seq_obj, int min_seed_coverage, std::list<PositionType> &ext
  );  
  void GetSeedExt(
      SequenceBuild &seq_obj, int min_seed_coverage,
      std::unordered_map<std::string, std::list<PositionType> > &ext
  );
  void GetSeedExtRevNew(
    SequenceBuild &seq_obj, int min_seed_coverage, std::list<PositionType> &rev_ext
  );  
  void GetSeedExtRev(
      SequenceBuild &rev_seq_obj, int min_seed_coverage,
      std::unordered_map<std::string, std::list<PositionType> > &rev_ext
  );
  void MatchSeedPairNew(
      SequenceBuild &seq_obj,
      std::list<PositionType> &ext, std::list<PositionType> &rev_ext,
      std::list<MinReadPairType> &ext_pair
  );
  /*
  void MatchSeedPair(
      SequenceBuild &seq_obj,
      std::unordered_map<std::string, std::list<PositionType> > &ext,
      std::unordered_map<std::string, std::list<PositionType> > &rev_ext,
      std::unordered_map<std::string, std::list<ReadPairType> >& ext_pair
  );
  */
  void CreateReadExtWorker(
      int min_ext_coverage, SequenceBuild& seq_obj, 
      std::unordered_map<RIDType, std::list<OverlapType> >& ext_read
  );
  
 protected:
  bool is_smm_built_; // check if the seed_mer_map is built
  bool is_alph_set_;
  bool is_seed_len_set_;
  bool is_overlap_len_set_;
  int alph_id_;
  int seed_len_;
  int overlap_len_;
  std::unordered_map<std::string, PositionType> seed_mer_map_;
  bool IsExtRedundant(
      std::vector<std::list<PositionType> >& fw_ext, 
      std::vector<std::list<PositionType> >& re_ext, 
      std::unordered_map<RIDType, std::set<int> >& read_map_table, 
      std::list<PositionType>& fw_phase, std::list<PositionType>& re_phase, 
      int& map_ID
  );
  bool IsPositionListMatch(
      std::list<PositionType>& l1, std::list<PositionType>& l2
  );
  bool IsSeqCompatible(
      SequenceBuild& seq_obj, int seed_len,
      RIDType fw_rid, POSType fw_pos,
      RIDType re_rid, POSType re_pos
  );
  /*
  void MatchSeedPairSingle(
      SequenceBuild& seq_obj,
      std::list<PositionType>& fw_sp, std::list<PositionType>& re_sp,
      std::list<MinReadPairType>& read_pair
  );
  */
};

inline uint16_t DatabaseIndex::Convert3mer(const char *s)  {
  //std::cout << "<<<<<<<<<<<<<<<<<<<" << std::endl;
  uint16_t c = 0B0000000000000000;
  for(int i = 0; i < 3; ++ i)  {
    uint16_t d;
    // '*' denotes the stop codon
    if(s[i] == '*') {
      d = 27;
    } else  {
      uint16_t o = (uint16_t) s[i];
      d = o - 64;  // 'A' = 1, 'B '= 2, ...
    }
    // bit-wise OR
    //std::bitset<16> y(d);
    //char r = s[i];
    //std::bitset<8> z(r);
    //std::cout << "char: " << s[i] << " " << y << "  " << z << std::endl;
    c = c | d;
    if(i < 2)  {
      c = c << 5;
    }
    //std::bitset<16> x(c);
    //std::cout << x << std::endl;
  }
  return c;
}

inline void DatabaseIndex::Interpret3merIndex(const uint16_t c, char *s)  {
  
  uint16_t t = c;
  //std::bitset<16> y(t);
  //std::cout << "Input:  " << y << std::endl;
  for(int i = 0; i < 3; ++ i)  {
    uint16_t d = 31; // 0000000000000111
    d = t & d;
    // '*' denotes the stop codon
    if(d == 27) {
      s[2 - i] = '*';
    } else  {
      uint16_t o = d + 64;
      s[2 - i] = (char) o;
    }
    
    //uint16_t r = d + 64;
    //char a = (char) r;
    //std::bitset<16> z1(r);
    //std::bitset<8> z2(a);
    //std::cout << "char: " << z1 << "  " << z2 << std::endl;
    t = t >> 5;
    //std::bitset<16> x(t);
    //std::cout << x << std::endl;
  }
  s[3] = '\0';
  //std::cout << ">>>>>>>>>>>>>" << std::endl;
  return;
}

#endif
