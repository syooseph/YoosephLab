#include "database_index.h"

using namespace std;

DatabaseIndex::DatabaseIndex(int in_alph_id, int in_seed_len, int in_overlap_len)  {
  is_smm_built_ = false;
  alph_id_ = in_alph_id;
  is_alph_set_ = true;
  seed_len_ = in_seed_len;
  is_seed_len_set_ = true;
  overlap_len_ = in_overlap_len;
  is_overlap_len_set_ = true;
  return;
}

DatabaseIndex::DatabaseIndex(void)  {
  is_smm_built_ = false;
  is_alph_set_ = false;
  is_seed_len_set_ = false;
  is_overlap_len_set_ = false;
  return;
}

DatabaseIndex::~DatabaseIndex(void) {
  return;
}

int DatabaseIndex::GetSeedLen(void) {
  if(!is_seed_len_set_)  {
    cout << "Error: DatabaseIndex::GetSeedLen: seed length is not set" << endl;
    exit(0);
  }
  return seed_len_;
}
  
int DatabaseIndex::GetOverlapLen(void)  {
  if(!is_overlap_len_set_)  {
    cout << "Error: DatabaseIndex::GetOverlapLen: overlap length is not set" << endl;
    exit(0);
  }
  return overlap_len_;
}

int DatabaseIndex::GetAlphID(void)  {
  if(!is_alph_set_)  {
    cout << "Error: DatabaseIndex::GetAlphID: alphabet is not set" << endl;
    exit(0);
  }
  return alph_id_;
}

void DatabaseIndex::BuildSeedmerMap(SequenceBuild& seq_obj) {
  // go over the suffix array and find all the kmers and their position
  if(!seq_obj.is_sfa_built_)  {
    cout << "Warning: DatabaseIndex::BuildSeedmerMap suffix array not built for input object" << endl;
  }
  if(!seq_obj.is_k_array_built_)  {
    cout << "Warning: DatabaseIndex::BuildSeedmerMap key array not built for input object" << endl;
  }
  if(!is_seed_len_set_ || !is_overlap_len_set_ || !is_alph_set_)  {
    cout << "Error: DatabaseIndex::BuildSeedmerMap: object not initialized with seed length/overlap length, only load functions are available" << endl;
    exit(0);
  }
  IdxType n = seq_obj.suffix_array_->getSize();
  IdxType i = 0;
  while(i < n) {
    if(seq_obj.suffix_array_->getSuffixLength(i) >= seed_len_)  {
      // try to see if the seed_mer_is present in the map or not
      string mer = seq_obj.suffix_array_->getSuffix(i);
      mer = mer.substr(0, seed_len_);
      //cout << "BuildSeedmerMap: seed-mer: " << mer << endl;
      auto it = seed_mer_map_.find(mer);
      if(it == seed_mer_map_.end())  { // when the seed is not presented in the map
        PositionType mer_pos;
        mer_pos.rid = (RIDType) seq_obj.suffix_array_->getId(i);
        mer_pos.pos = (POSType) seq_obj.suffix_array_->getPos(i);
        seed_mer_map_.insert({mer, mer_pos});
      }
      // jump to the next chunck
      i = seq_obj.key_array_[i];
    }
    ++ i;
  }
  is_smm_built_ = true;
  return;
}

void DatabaseIndex::CreateReducedMap(
    std::unordered_map<std::string, std::list<std::string> >& reduc_alph_map
) {
  // create the reudce alphabet map, where for each reduced alphabet list 
  // all sequences in the original alphabet
  if(!is_smm_built_)  {
    cout << "Error: DatabaseIndex::CreateReducedMap seed-mer map is not built in this object" << endl;
    exit(0);
  }
  if(!is_seed_len_set_ || !is_overlap_len_set_ || !is_alph_set_)  {
    cout << "Error: DatabaseIndex::CreateReducedMap: object not initialized with reduced alphabet/seed length/overlap length, only load functions are available" << endl;
    exit(0);
  }
  ReducedAlphabet reduc_alph_obj((Alphabet) alph_id_);
  // convert each seed-mer in reduced alphabet and record their positions
  for(auto it = seed_mer_map_.begin(); it != seed_mer_map_.end(); ++ it) {
    bool is_standard;
    //cout << "CreateReducedMap: seed-mer:  " << it->first << endl;
    string reduc_mer = reduc_alph_obj.Convert(it->first, is_standard);
    if(is_standard)  {
      reduc_alph_map[reduc_mer].push_back(it->first);
    }
  }
  // debug printing
  //for(auto it = reduc_alph_map.begin(); it != reduc_alph_map.end(); ++ it) {
  //  cout << it->first << endl;
  //  for(auto it_p = it->second.begin(); it_p != it->second.end(); ++ it_p) {
  //    cout << *it_p << ";";
  //  }
  //  cout << endl;
  //}
  return;
}

void DatabaseIndex::DumpReducedMap(
    std::unordered_map<std::string, std::list<std::string> >& reduc_alph_map,
    std::string& out_file
) {
  // create and check validity of the output file
  ofstream out_fh(out_file.c_str(), ios_base::out | ios_base::binary);
  if(!out_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::DumpReducedMap cannot write index file " << out_file << endl;
    exit(1);
  }
  // write config information: seed length, size_of_RIDType, size_of_POSType, total entries
  int total_entries = 0;
  out_fh.write((char*) &seed_len_, sizeof(int));
  out_fh.write((char*) &alph_id_, sizeof(int));
  streampos total_size_pos = out_fh.tellp();  // record the position for size
  out_fh.write((char*) &total_entries, sizeof(int));
  // write stored information from the input map
  for(auto it = reduc_alph_map.begin(); it != reduc_alph_map.end(); ++ it) {
    // write the seed-mer in reduced alphabet
    out_fh.write((char*) it->first.c_str(), sizeof(char) * seed_len_);
    // write the number of entries for this seed-mer
    RIDType num_entries = (RIDType) it->second.size();
    out_fh.write((char*) &num_entries, sizeof(RIDType));
    // write each position
    for(auto it_p = it->second.begin(); it_p != it->second.end(); ++ it_p) {
      out_fh.write((char*) it_p->c_str(), sizeof(char) * seed_len_);
      ++ total_entries;
    }
  }
  // update the number of entries in the file
  out_fh.seekp(total_size_pos);
  out_fh.write((char*) &total_entries, sizeof(int));
  out_fh.close();
  return;
}

void DatabaseIndex::LoadReducedMap(
    std::string& in_file, 
    std::unordered_map<std::string, std::unordered_map<std::string, bool> >& reduc_alph_map
) {
  // try open and check validity of the input file
  ifstream in_fh(in_file.c_str(), ios_base::in);
  if(!in_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::LoadReducedMap cannot read index file " << in_file << endl;
    exit(1);
  }
  // load basic configurations: seed length, size_of_RIDType, size_of_POSType
  int loaded_seed_len, loaded_alph_id, loaded_num_entries;
  in_fh.read((char*) &loaded_seed_len, sizeof(int));
  in_fh.read((char*) &loaded_alph_id, sizeof(int));
  alph_id_ = loaded_alph_id;
  is_alph_set_ = true;
  in_fh.read((char*) &loaded_num_entries, sizeof(int));
  if(loaded_num_entries <= 0)  {
    in_fh.close();
    return;
  }
  //cout << "loaded into: " << loaded_seed_len << " " << loaded_alph_id << "  " << loaded_num_entries << endl;
  // load information per seed-mer
  char* mer_cstr = new char[loaded_seed_len + 1];
  while(!in_fh.eof())  {
    // read seed_mer  
    in_fh.read((char*) mer_cstr, sizeof(char) * loaded_seed_len);
    mer_cstr[loaded_seed_len] = '\0';
    string mer = mer_cstr;
    //cout << mer << endl;
    // read num of entries
    int num_entries;
    in_fh.read((char*) &num_entries, sizeof(RIDType));
    for(int i = 0; i < num_entries; ++ i) {
      in_fh.read((char*) mer_cstr, sizeof(char) * loaded_seed_len);
      mer_cstr[loaded_seed_len] = '\0';
      string mer_ori = mer_cstr;
      reduc_alph_map[mer][mer_ori] = true;
      //cout << mer_ori << ",";
      -- loaded_num_entries;
    }
    //cout << endl;
    //cout << loaded_num_entries << endl;
    if(loaded_num_entries <= 0)  {
      break;
    }
  }
  delete [] mer_cstr;
  if(loaded_num_entries < 0)  {
    cout << "Warning: DatabaseIndex::LoadReducedMap corrupted file, please rebuild index " << endl;
  }
  in_fh.close();
  return;
}

bool _cmp_position(const PositionType& a, const PositionType& b) {
  if(a.rid < b.rid)  {
    return true;
  }
  return false;
}

bool DatabaseIndex::IsPositionListMatch(
    std::list<PositionType>& l1, std::list<PositionType>& l2
) {
  if(l1.size() != l2.size())  {
    return false;
  }
  auto it_l2 = l2.begin();
  for(auto it_l1 = l1.begin(); it_l1 != l1.end(); ++ it_l1) {
    if(it_l1->rid != it_l2->rid)  {
      return false;
    }
    ++ it_l2;
  }
  return true;
}

bool DatabaseIndex::IsExtRedundant(
    std::vector<std::list<PositionType> >& fw_ext, 
    std::vector<std::list<PositionType> >& re_ext, 
    std::unordered_map<RIDType, std::set<int> >& read_map_table, 
    std::list<PositionType>& fw_phase, std::list<PositionType>& re_phase, 
    int& map_ID
) {
  set<int> fw_pos, re_pos;
  // find out where the reads are mapped to
  for(auto it = fw_phase.begin(); it != fw_phase.end(); ++ it) {
    auto it_e = read_map_table.find(it->rid);
    if(it_e != read_map_table.end())  {
      for(auto it_p = it_e->second.begin(); it_p != it_e->second.end(); ++ it_p) {
        fw_pos.insert(*it_p);
      }
    }
  }
  for(auto it = re_phase.begin(); it != re_phase.end(); ++ it) {
    auto it_e = read_map_table.find(it->rid);
    if(it_e != read_map_table.end())  {
      for(auto it_p = it_e->second.begin(); it_p != it_e->second.end(); ++ it_p) {
        re_pos.insert(*it_p);
      }
    }
  }
  // find intersectio of the positions
  set<int> pos_intersect;
  for(auto it = re_pos.begin(); it != re_pos.end(); ++ it) {
    if(fw_pos.find(*it) != fw_pos.end())  {
      pos_intersect.insert(*it);
    }
  }
  // try to find left and right extensions that are exactly the same 
  for(auto it = pos_intersect.begin(); it != pos_intersect.end(); ++ it) {
    if(IsPositionListMatch(fw_ext[*it], fw_phase) && 
        IsPositionListMatch(re_ext[*it], re_phase))  {
      map_ID = *it;
      return true;
    }
  }
  return false;
}

void DatabaseIndex::CreateSeedExt(
    int min_seed_coverage, SequenceBuild& seq_obj, SequenceBuild& rev_seq_obj,
    std::unordered_map<std::string, std::list<ReadPairType> >& ext_seed
) {
  // check if the seed-mer map is built
  if(!is_smm_built_)  {
    cout << "Error: DatabaseIndex::CreateReducedMap seed-mer map is not built in this object" << endl;
    exit(0);
  }
  if(!is_seed_len_set_ || !is_overlap_len_set_)  {
    cout << "Error: DatabaseIndex::CreateSeedExt: object not initialized with seed length/overlap length, only load functions are available" << endl;
    exit(0);
  }
  // foreach seed mer, search both suffix arrays and record the maximal extension sequences
  //cout << "number of different k-mers:  " << seed_mer_map_.size() << endl;
  for(auto it = seed_mer_map_.begin(); it != seed_mer_map_.end(); ++ it) {
    list<PositionType> fw_ext_phase, re_ext_phase;
    // search for the forward sequence
    string seed = it->first;
    //cout << "Seed" << seed << endl;
    pair<IdxType, IdxType> fw_range = seq_obj.SearchSFA(seed);
    if(fw_range.second - fw_range.first + 1 < min_seed_coverage)  {
    //  cout << "Not enough coverage, skipped fw" << endl;
      continue;
    }
    seq_obj.GetMaxExtInfoWithinRange(fw_range, fw_ext_phase);
    // search for the reverse sequence
    string rev_seed(seed.rbegin(), seed.rend());
    pair<IdxType, IdxType> re_range = rev_seq_obj.SearchSFA(rev_seed);
    if(re_range.second - re_range.first + 1 < min_seed_coverage)  {
    //  cout << "Not enough coverage, skipped re" << endl;
      continue;
    }
    rev_seq_obj.GetMaxExtInfoWithinRange(re_range, re_ext_phase);  
    // recompute the position of the reverse seeds
    for(auto it_s = re_ext_phase.begin(); it_s != re_ext_phase.end(); ++ it_s) {
      it_s->pos = (POSType) 
          (strlen(seq_obj.sequence_[it_s->rid]) - seed_len_ - (int) it_s->pos);
    }
    cout << "Match found" << endl;
    MatchSeedPairSingle(seq_obj, fw_ext_phase, re_ext_phase, ext_seed[seed]);
    //fw_ext_phase.clear();
    //re_ext_phase.clear();
    //cout << "sizes: " << fw_ext_phase.size() << " " << re_ext_phase.size() << " " << ext_seed[seed].size() << endl;
  }
  return;
}

void DatabaseIndex::GetSeedExt(
    SequenceBuild &seq_obj, int min_seed_coverage,
    std::unordered_map<std::string, std::list<PositionType> > &ext
) {
  // check if the seed-mer map is built
  if(!is_smm_built_)  {
    cout << "Error: DatabaseIndex::CreateReducedMap seed-mer map is not built in this object" << endl;
    exit(0);
  }
  if(!is_seed_len_set_ || !is_overlap_len_set_)  {
    cout << "Error: DatabaseIndex::CreateSeedExt: object not initialized with seed length/overlap length, only load functions are available" << endl;
    exit(0);
  }
  for(auto it = seed_mer_map_.begin(); it != seed_mer_map_.end(); ++ it) {
    list<PositionType> fw_ext_phase;
    // search for the forward sequence
    string seed = it->first;
    //cout << "Seed" << seed << endl;
    pair<IdxType, IdxType> fw_range = seq_obj.SearchSFA(seed);
    if(fw_range.second + 1 < (IdxType) min_seed_coverage + fw_range.first)  {
    //  cout << "Not enough coverage, skipped fw" << endl;
      continue;
    }
    seq_obj.GetMaxExtInfoWithinRange(fw_range, fw_ext_phase);
    ext[seed] = fw_ext_phase;
  }
  return;
}

void DatabaseIndex::GetSeedExtRev(
    SequenceBuild &rev_seq_obj, int min_seed_coverage,
    std::unordered_map<std::string, std::list<PositionType> > &rev_ext
) {
  // check if the seed-mer map is built
  if(!is_smm_built_)  {
    cout << "Error: DatabaseIndex::CreateReducedMap seed-mer map is not built in this object" << endl;
    exit(0);
  }
  if(!is_seed_len_set_ || !is_overlap_len_set_)  {
    cout << "Error: DatabaseIndex::CreateSeedExt: object not initialized with seed length/overlap length, only load functions are available" << endl;
    exit(0);
  }
  for(auto it = seed_mer_map_.begin(); it != seed_mer_map_.end(); ++ it) {
    list<PositionType> re_ext_phase;
    // search for the forward sequence
    string seed = it->first;
    string rev_seed(seed.rbegin(), seed.rend());
    pair<IdxType, IdxType> re_range = rev_seq_obj.SearchSFA(rev_seed);
    if(re_range.second + 1 < (IdxType) min_seed_coverage + re_range.first)  {
    //  cout << "Not enough coverage, skipped re" << endl;
      continue;
    }
    rev_seq_obj.GetMaxExtInfoWithinRange(re_range, re_ext_phase);  
    // recompute the position of the reverse seeds
    for(auto it_s = re_ext_phase.begin(); it_s != re_ext_phase.end(); ++ it_s) {
      it_s->pos = (POSType) 
          (strlen(rev_seq_obj.sequence_[it_s->rid]) - seed_len_ - (int) it_s->pos);
    }
    rev_ext[seed] = re_ext_phase;
  }
  return;
}

void DatabaseIndex::MatchSeedPair(
    SequenceBuild &seq_obj,
    std::unordered_map<std::string, std::list<PositionType> > &ext,
    std::unordered_map<std::string, std::list<PositionType> > &rev_ext,
    std::unordered_map<std::string, std::list<ReadPairType> >& ext_pair
) {
  if(!is_smm_built_)  {
    cout << "Error: DatabaseIndex::CreateReducedMap seed-mer map is not built in this object" << endl;
    exit(0);
  }
  for(auto it = seed_mer_map_.begin(); it != seed_mer_map_.end(); ++ it) {
    string seed = it->first;
    //cout << seed << endl;
    if(ext.find(seed) != ext.end() && rev_ext.find(seed) != rev_ext.end())  {
      //cout << " Match found" << endl;
      MatchSeedPairSingle(seq_obj, ext[seed], rev_ext[seed], ext_pair[seed]);
    }
  }
  return;
}

bool _cmp_pos_fw(PositionType& a, PositionType& b) {
  if(a.pos > b.pos)  {
    return true;
  }
  return false;
}

bool _cmp_pos_re(PositionType& a, PositionType& b) {
  if(a.pos < b.pos)  {
    return true;
  }
  return false;
}

bool DatabaseIndex::IsSeqCompatible(
    SequenceBuild& seq_obj, int seed_len,
    RIDType fw_rid, POSType fw_pos,
    RIDType re_rid, POSType re_pos
) {
  //cout << "Comparing sequences" << endl;
  //cout << seq_obj.sequence_[fw_rid] << "  " << (int) fw_pos << endl;
  //cout << seq_obj.sequence_[re_rid] << "  " << (int) re_pos <<endl;

  string fw_seq = seq_obj.sequence_[fw_rid];
  string re_seq = seq_obj.sequence_[re_rid];
  int ixf = (int) fw_pos - 1;
  int ixr = (int) re_pos - 1;
  while(ixf >= 0 && ixr >= 0) {
    if(fw_seq[ixf] != re_seq[ixr])  {
      //cout << "false" << endl;
      return false;
    }
    -- ixf;
    -- ixr;
  }
  ixf = (int) fw_pos + seed_len;
  ixr = (int) re_pos + seed_len;
  while(ixf < (int) fw_seq.length() && ixr < (int) re_seq.length()) {
    if(fw_seq[ixf] != re_seq[ixr])  {
      //cout << "false" << endl;
      return false;
    }
    ++ ixf;
    ++ ixr;
  }
  //cout << "true" << endl;
  return true;
}

void DatabaseIndex::MatchSeedPairSingle(
    SequenceBuild& seq_obj,
    std::list<PositionType>& fw_sp, std::list<PositionType>& re_sp,
    std::list<ReadPairType>& read_pair
) {
  fw_sp.sort(_cmp_pos_fw);
  re_sp.sort(_cmp_pos_re);
  unordered_map<RIDType, bool> fw_taken; 
  unordered_map<RIDType, bool> re_taken;
  // try finding seed read pairs
  
  //cout << "*****************************" << endl;
  //for(auto itp = fw_sp.begin(); itp != fw_sp.end(); ++ itp) {
  //  cout << "FW:  " << itp->rid << endl;
  //}
  //for(auto itp = re_sp.begin(); itp != re_sp.end(); ++ itp) {
  //  cout << "RE:  " << itp->rid << endl;
  //}
  
  auto it_fw = fw_sp.begin();
  while(it_fw != fw_sp.end()) {
    //cout << "in merge loop lv1" << endl;
    if(fw_taken.find(it_fw->rid) == fw_taken.end())  {
      auto it_re = re_sp.begin();
      while(it_re != re_sp.end()) {
        //cout << "in merge loop lv2" << endl;
        //cout << "//////////////////////////////" << endl;
        //cout << "$$$: " << seq_obj.sequence_[it_fw->rid] << " " << (int) it_fw->pos << "  " << seq_obj.GetHeader(it_fw->rid) << endl;
        //cout << "$$$: " << seq_obj.sequence_[it_re->rid] << " " << (int) it_re->pos << "  " << seq_obj.GetHeader(it_re->rid) << endl;
        if(re_taken.find(it_re->rid) == re_taken.end())  {
          bool is_merge_success = false;
          ReadPairType rp;
          rp.rid_fw = it_fw->rid;
          rp.r_pos_fw = it_fw->pos;
          rp.rid_re = it_re->rid;
          rp.r_pos_re = it_re->pos;
          rp.init_fw = rp.init_re = true;
          rp.overlap = strlen(seq_obj.sequence_[it_re->rid]) - it_re->pos + it_fw->pos;
          if(rp.rid_fw == rp.rid_re)  {
            // if the same read is taken
            is_merge_success = true;
          } else if(
              IsSeqCompatible(seq_obj, seed_len_, it_fw->rid, it_fw->pos, it_re->rid, it_re->pos)
            ) 
          {
            // if the sequence is compatible
            if(rp.overlap >= 12)  {
              // if overlap is significant enough
              is_merge_success = true;
            } else if(it_re->pos >= 3 && strlen(seq_obj.sequence_[it_fw->rid]) - seed_len_ > 3)  {
            // we need to search the the suffix array to find supporting bridging reads
              string search_seq = string(seq_obj.sequence_[it_re->rid]).substr(it_re->pos - 3, 3)
                + string(seq_obj.sequence_[it_fw->rid]).substr(it_fw->pos, seed_len_ + 3);
              //cout << "overlap: " << rp.overlap << endl;
              //cout << "Search_seq:  " << search_seq << endl;
              pair<IdxType, IdxType> range = seq_obj.SearchSFA(search_seq);
              if(range.second > range.first)  {
                is_merge_success = true;
              }
            }
          }
          if(is_merge_success)  {
          //  cout << "Pair matched:" << seq_obj.GetHeader(it_fw->rid) << "  " << seq_obj.GetHeader(it_re->rid) << endl;
            read_pair.push_back(rp);
            fw_taken[it_fw->rid] = true;
            re_taken[it_re->rid] = true;
            break;
          }
        }
        ++ it_re;
      }
    }
    ++ it_fw;
  }
  // recruit the rest of the unmatched ones
  auto it = fw_sp.begin();
  while(it != fw_sp.end()) {
    //cout << "in fw loop lv1" << endl;
    if(fw_taken.find(it->rid) == fw_taken.end())  {
      ReadPairType rp;
      rp.rid_fw = it->rid;
      rp.r_pos_fw = it->pos;
      rp.init_fw = true;
      rp.init_re = false;
      rp.overlap = 0;
      read_pair.push_back(rp);
      fw_taken[it->rid] = true;
      //cout << "Pair sole fw:" << seq_obj.GetHeader(it->rid) << endl;
    }
    ++ it;
  }
  it = re_sp.begin();
  while(it != re_sp.end()) {
    //cout << "in re loop lv1" << endl;
    if(re_taken.find(it->rid) == re_taken.end())  {
      ReadPairType rp;
      rp.rid_re = it->rid;
      rp.r_pos_re = it->pos;
      rp.init_fw = false;
      rp.init_re = true;
      rp.overlap = 0;
      read_pair.push_back(rp);
      re_taken[it->rid] = true;
      //cout << "Pair sole re:" << seq_obj.GetHeader(it->rid) << endl;
    }
    ++ it;
  }
  return;
}

void DatabaseIndex::DumpSeedExt(
    std::unordered_map<std::string, std::list<ReadPairType> >& ext_seed,
    std::string& out_file
) {
  // create and check validity of the output file
  ofstream out_fh(out_file.c_str(), ios_base::out | ios_base::binary);
  if(!out_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::DumpSeedExt cannot write index file " << out_file << endl;
    exit(1);
  }
  // write config information: seed length, size_of_RIDType, size_of_POSType, total entries
  int total_entries = 0;
  out_fh.write((char*) &seed_len_, sizeof(int));
  streampos total_size_pos = out_fh.tellp();  // record the position for size
  out_fh.write((char*) &total_entries, sizeof(int));
  // write stored information from the input map
  for(auto it = ext_seed.begin(); it != ext_seed.end(); ++ it) {
    // write the seed-mer in reduced alphabet
    out_fh.write((char*) it->first.c_str(), sizeof(char) * seed_len_);
    // write number of forward and reverse extensions, respectiv,y
    RIDType num_entries = (RIDType) ext_seed[it->first].size();
    out_fh.write((char*) &num_entries, sizeof(RIDType));
    // write the actual Info
    for(auto it_pf = ext_seed[it->first].begin(); 
        it_pf != ext_seed[it->first].end(); ++ it_pf
    ) {
      out_fh.write((char*) &(it_pf->rid_fw), sizeof(RIDType));
      out_fh.write((char*) &(it_pf->r_pos_fw), sizeof(POSType));
      out_fh.write((char*) &(it_pf->rid_re), sizeof(RIDType));
      out_fh.write((char*) &(it_pf->r_pos_re), sizeof(POSType));
      out_fh.write((char*) &(it_pf->init_fw), sizeof(bool));
      out_fh.write((char*) &(it_pf->init_re), sizeof(bool));
      out_fh.write((char*) &(it_pf->overlap), sizeof(int));
      ++ total_entries;
    }
  }
  // update the number of entries in the file
  out_fh.seekp(total_size_pos);
  //cout << "num total entries: " << total_entries << endl;
  out_fh.write((char*) &total_entries, sizeof(int));
  out_fh.close();
  return;
}

void DatabaseIndex::LoadSeedExt(
    std::string& in_file,
    std::unordered_map<std::string, std::list<ReadPairType> >& ext_seed
) {
  // try open and check validity of the input file
  ifstream in_fh(in_file.c_str(), ios_base::in);
  if(!in_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::LoadSeedExt cannot read index file " << in_file << endl;
    exit(1);
  }
  // load basic configurations: seed length, size_of_RIDType, size_of_POSType
  int loaded_seed_len, loaded_num_entries;
  in_fh.read((char*) &loaded_seed_len, sizeof(int));
  seed_len_ = loaded_seed_len;
  is_seed_len_set_ = true;
  in_fh.read((char*) &loaded_num_entries, sizeof(int));
  if(loaded_num_entries <= 0)  {
    in_fh.close();
    return;
  }
  //cout << "loaded num entires:  " << loaded_num_entries << endl;
  // load information per seed-mer
  char* mer_cstr = new char[loaded_seed_len + 1];
  while(!in_fh.eof())  {
    // read seed_mer  
    in_fh.read((char*) mer_cstr, sizeof(char) * loaded_seed_len);
    mer_cstr[loaded_seed_len] = '\0';
    string mer = mer_cstr;
    // read num of entries
    int num_entries;
    in_fh.read((char*) &num_entries, sizeof(RIDType));
    for(int i = 0; i < num_entries; ++ i) {
      ReadPairType sp;
      in_fh.read((char*) &sp.rid_fw, sizeof(RIDType));
      POSType pf, pr;
      in_fh.read((char*) &pf, sizeof(POSType));
      sp.r_pos_fw = (int) pf;
      in_fh.read((char*) &sp.rid_re, sizeof(RIDType));
      in_fh.read((char*) &pr, sizeof(POSType));
      sp.r_pos_re = (int) pr;
      in_fh.read((char*) &sp.init_fw, sizeof(bool));
      in_fh.read((char*) &sp.init_re, sizeof(bool));
      in_fh.read((char*) &sp.overlap, sizeof(int));
      //cout << "fw positions loaded:  " << (unsigned int) loc.pos << endl;
      ext_seed[mer].push_back(sp);
      -- loaded_num_entries;
    }
    if(loaded_num_entries <= 0)  {
      break;
    }
  }
  delete [] mer_cstr;
  if(loaded_num_entries < 0)  {
    cout << "Warning: DatabaseIndex::LoadSeedExt corrupted file, please rebuild index " << endl;
  }
  in_fh.close();
  return;
}

// build the index for read extension
void DatabaseIndex::CreateReadExt(
    int min_ext_coverage,
    SequenceBuild& seq_obj, 
    SequenceBuild& rev_seq_obj,
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read
) {
  //cout << "CreateReadExt Forward: *****************************" << endl;
  CreateReadExtWorker(min_ext_coverage, seq_obj, fw_ext_read);
  //cout << "CreateReadExt Reverse: *****************************" << endl;
  CreateReadExtWorker(min_ext_coverage, rev_seq_obj, re_ext_read);
  return;
}

void DatabaseIndex::CreateReadExtWorker(
    int min_ext_coverage, SequenceBuild& seq_obj, 
    std::unordered_map<RIDType, std::list<OverlapType> >& ext_read
) {
  if(!is_seed_len_set_ || !is_overlap_len_set_)  {
    cout << "Error: DatabaseIndex::CreateReadExtWorker: object not initialized with seed length/overlap length, only load functions are available" << endl;
    exit(0);
  }
  IdxType i = 1, n = seq_obj.suffix_array_->getSize();
  while(i < n) {
    if(seq_obj.suffix_array_->getSuffixLength(i) == overlap_len_)  {
      // first find all reads that ends with the overlap-mer
      //cout << "####################" << endl;
      list<RIDType> src_reads;
      do  {
        src_reads.push_back(seq_obj.suffix_array_->getId(i));
        //cout << "size:  " << seq_obj.num_seqs_ << "  " << i << endl;
        //cout << "src: " << seq_obj.suffix_array_->getId(i) << " " << seq_obj.sequence_[seq_obj.suffix_array_->getId(i)] << endl;
        //cout << "src seq: " << seq_obj.suffix_array_->getSuffix(i) << endl;
        ++ i;
      } while (seq_obj.suffix_array_->getSuffixLength(i) == overlap_len_ // the length of the suffix
          && seq_obj.suffix_array_->getLcp(i) == overlap_len_ // the sequence must be the same
      );
      -- i;
      // then find all reads that begin with/contain the overlap-mer
      list<OverlapType> overlaps;
      do {
        //for(int ii = i; ii <= seq_obj.key_array_[i] + 1; ++ ii) {
        //  cout << "sa seq:  " << seq_obj.suffix_array_->getSuffix(ii) << endl;
        //}
        IdxType old_i = i;
        i = seq_obj.key_array_[i] + 1;
        // record the extension information
        if(i - old_i > min_ext_coverage && seq_obj.suffix_array_->getSuffixLength(i - 1) > overlap_len_)  {
          OverlapType olp;
          olp.rid = seq_obj.suffix_array_->getId(i - 1);
          olp.len = seq_obj.suffix_array_->getPos(i - 1) + overlap_len_;
          overlaps.push_back(olp);
          //cout << "tgt: " << seq_obj.suffix_array_->getId(i - 1) << endl;
          //cout << "tgt seq: " << seq_obj.suffix_array_->getSuffix(i - 1) << endl;
        }
      } while(i < seq_obj.suffix_array_->getSize()  
          && seq_obj.suffix_array_->getLcpAt(i) >= overlap_len_
      );
      -- i;
      // record information if extension is possible
      if(!overlaps.empty())  {
        for(auto it = src_reads.begin(); it != src_reads.end(); ++ it) {
          for(auto itt = overlaps.begin(); itt != overlaps.end(); ++ itt) {
            ext_read[*it].push_back(*itt);
          }
        }
      }
    }
    ++ i;
  }
  return;
}


// write the extension inforamtion
void DatabaseIndex::DumpReadExt(
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read,
    std::string& out_file
) {
  //cout << "DumpReadExt called:  " << endl;
  // create and check validity of the output file
  ofstream out_fh(out_file.c_str(), ios_base::out | ios_base::binary);
  if(!out_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::DumpReadExt cannot write index file " << out_file << endl;
    exit(1);
  }
  int total_entries = 0;
  out_fh.write((char*) &overlap_len_, sizeof(int));
  streampos total_size_pos = out_fh.tellp();  // record the position for size
  out_fh.write((char*) &total_entries, sizeof(int));
  //cout << "Info:  " << overlap_len_ << "  " << total_entries << endl;
  // write stored information from the input map
  for(auto it = fw_ext_read.begin(); it != fw_ext_read.end(); ++ it) {
    // write the source read ID
    out_fh.write((char*) &(it->first), sizeof(RIDType));
    //cout << "********:  " << (int) it->first << endl;
    // write number of forward and reverse extensions, respectiv,y
    RIDType num_entries_fw = (RIDType) fw_ext_read[it->first].size();
    out_fh.write((char*) &num_entries_fw, sizeof(RIDType));
    RIDType num_entries_re = 0;
    out_fh.write((char*) &num_entries_re, sizeof(RIDType));
    // write the actual Info
    for(auto it_pf = fw_ext_read[it->first].begin(); it_pf != fw_ext_read[it->first].end(); ++ it_pf) {
      out_fh.write((char*) &(it_pf->rid), sizeof(RIDType));
      out_fh.write((char*) &(it_pf->len), sizeof(POSType));
      //cout << (int) it_pf->rid << " " << (int) it_pf->len << endl;
      ++ total_entries;
    }
    //for(auto it_pr = re_ext_read[it->first].begin(); it_pr != re_ext_read[it->first].end(); ++ it_pr) {
    //  out_fh.write((char*) &(it_pr->rid), sizeof(RIDType));
    //  out_fh.write((char*) &(it_pr->len), sizeof(POSType));
      //cout << (int) it_pr->rid << " " << (int) it_pr->len << endl;
    //  ++ total_entries;
    //}
  }
  for(auto it = re_ext_read.begin(); it != re_ext_read.end(); ++ it) {
    // write the source read ID
    out_fh.write((char*) &(it->first), sizeof(RIDType));
    //cout << "********:  " << (int) it->first << endl;
    // write number of forward and reverse extensions, respectiv,y
    RIDType num_entries_fw = (RIDType) 0;
    out_fh.write((char*) &num_entries_fw, sizeof(RIDType));
    RIDType num_entries_re = (RIDType) re_ext_read[it->first].size();
    out_fh.write((char*) &num_entries_re, sizeof(RIDType));
    // write the actual Info
    //for(auto it_pf = fw_ext_read[it->first].begin(); it_pf != fw_ext_read[it->first].end(); ++ it_pf) {
    //  out_fh.write((char*) &(it_pf->rid), sizeof(RIDType));
    //  out_fh.write((char*) &(it_pf->len), sizeof(POSType));
      //cout << (int) it_pf->rid << " " << (int) it_pf->len << endl;
    //  ++ total_entries;
    //}
    for(auto it_pr = re_ext_read[it->first].begin(); it_pr != re_ext_read[it->first].end(); ++ it_pr) {
      out_fh.write((char*) &(it_pr->rid), sizeof(RIDType));
      out_fh.write((char*) &(it_pr->len), sizeof(POSType));
      //cout << (int) it_pr->rid << " " << (int) it_pr->len << endl;
      ++ total_entries;
    }
  }
  // update the number of entries in the file
  out_fh.seekp(total_size_pos);
  //cout << "num total entries: " << total_entries << endl;
  out_fh.write((char*) &total_entries, sizeof(int));
  out_fh.close();
  return;
}

void DatabaseIndex::LoadReadExt(
    std::string& in_file,
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read
) {
  // try open and check validity of the input file
  ifstream in_fh(in_file.c_str(), ios_base::in);
  if(!in_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::LoadSeedExt cannot read index file " << in_file << endl;
    exit(1);
  }
  // load basic configurations: seed length, size_of_RIDType, size_of_POSType
  int loaded_overlap_len, loaded_num_entries;
  in_fh.read((char*) &loaded_overlap_len, sizeof(int));
  overlap_len_ = loaded_overlap_len;
  is_overlap_len_set_ = true;
  in_fh.read((char*) &loaded_num_entries, sizeof(int));
  if(loaded_num_entries <= 0)  {
    in_fh.close();
    return;
  }
  //cout << "loaded num entires:  " << loaded_num_entries << endl;
  // load information per seed-mer
  RIDType src_rid;
  while(!in_fh.eof())  {
    // read seed_mer  
    in_fh.read((char*) &src_rid, sizeof(RIDType));
    // read num of entries
    //cout << "load src rid*****************************:  " << src_rid << endl;
    int num_entries_fw, num_entries_re;
    in_fh.read((char*) &num_entries_fw, sizeof(RIDType));
    in_fh.read((char*) &num_entries_re, sizeof(RIDType));
    for(int i = 0; i < num_entries_fw; ++ i) {
      OverlapType loc;
      in_fh.read((char*) &loc.rid, sizeof(RIDType));
      in_fh.read((char*) &loc.len, sizeof(POSType));
      //cout << "load fw rid:  " << loc.rid << endl;
      fw_ext_read[src_rid].push_back(loc);
      -- loaded_num_entries;
    }
    for(int i = 0; i < num_entries_re; ++ i) {
      OverlapType loc;
      in_fh.read((char*) &loc.rid, sizeof(RIDType));
      in_fh.read((char*) &loc.len, sizeof(POSType));
      //cout << "load re rid:  " << loc.rid << endl;
      re_ext_read[src_rid].push_back(loc);
      -- loaded_num_entries;
    }
    if(loaded_num_entries <= 0)  {
      break;
    }
  }
  if(loaded_num_entries < 0)  {
    cout << "Warning: DatabaseIndex::LoadReadExt corrupted file, please rebuild index " << endl;
  }
  in_fh.close();
  return;
}

void DatabaseIndex::CreateHighScoreMer(
    std::unordered_map<std::string, std::list<std::string> >& high_score_match
)  {
  int mer_size = 3;
  double def_high_frac = 0.3;
  vector<char> alphabet = {
      'P', 'G', 'E', 'K', 'R', 'Q', 'D', 'S', 'N', 'T', 
      'H', 'C', 'I', 'V', 'W', 'Y', 'F', 'A', 'L', 'M'
  };
  ReducedAlphabet reduc_alph = ReducedAlphabet(Alphabet(alph_id_)); 
  ScoringFunction<int> score_scheme(PROTEIN, BLOSUM62, -1, -11);
  // allocate strings
  vector<string> mer_enum;
  mer_enum.resize((int) pow(alphabet.size(), mer_size));
  for(int i = 0; i < (int) mer_enum.size(); ++ i) {
    string foo(mer_size, ' ');
    int num = i;
    //cout << "num: " << num << endl;
    for(int p = mer_size - 1; p >= 0; -- p) {
      int d = pow(alphabet.size(), p);
      int q = (int) (num / d);
      foo[mer_size - p - 1] = alphabet[q];
      //cout << " pow, quotion, char:  " << d << " " << q << "  " << alphabet[q] << endl;
      num = num - q * d; 
    }
    mer_enum[i] = foo;
    //cout << foo << endl;
  }  
  // compute score
  double score_cutoff = (double) def_high_frac * score_scheme.GetAveMatch() * mer_size;
  for(int i = 0; i < (int) mer_enum.size(); ++ i) {
    for(int j = i; j < (int) mer_enum.size(); ++ j) {
      bool foo;
      if((double) score_scheme.CalMatchScore(mer_enum[i], mer_enum[j]) >= score_cutoff &&
          reduc_alph.Convert((string) mer_enum[i], foo) == reduc_alph.Convert((string) mer_enum[j], foo) && 
          foo
      )  {
        high_score_match[mer_enum[i]].push_back(mer_enum[j]);
        high_score_match[mer_enum[j]].push_back(mer_enum[i]);
      }
    }
  }
  //for(auto it = high_score_match_.begin(); it != high_score_match_.end(); ++ it) {
  //  cout << "high-score size: " << it->second.size() << endl;
  //}
  return;
}

void DatabaseIndex::DumpHighScoreMer(
    std::unordered_map<std::string, std::list<std::string> >& high_score_match,
    std::string& out_file
) {
    //cout << "DumpReadExt called:  " << endl;
  // create and check validity of the output file
  ofstream out_fh(out_file.c_str(), ios_base::out | ios_base::binary);
  if(!out_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::DumpReadExt cannot write index file " << out_file << endl;
    exit(1);
  }
  int mer_len = 3;
  int total_entries = 0;
  out_fh.write((char*) &mer_len, sizeof(int));
  streampos total_size_pos = out_fh.tellp();  // record the position for size
  out_fh.write((char*) &total_entries, sizeof(int));
  //cout << "Info:  " << overlap_len_ << "  " << total_entries << endl;
  // write stored information from the input map
  for(auto it = high_score_match.begin(); it != high_score_match.end(); ++ it) {
    // write the source read ID
    out_fh.write((char*) it->first.c_str(), sizeof(char) * mer_len);
    //cout << "********:  " << it->first << endl;
    // write number of forward and reverse extensions, respectiv,y
    int num_entries = it->second.size();
    out_fh.write((char*) &num_entries, sizeof(int));
    // write the actual Info
    for(auto it_pf = it->second.begin(); it_pf != it->second.end(); ++ it_pf) {
      out_fh.write((char*) it_pf->c_str(), sizeof(char) * mer_len);
      //cout << *it_pf << endl;
      ++ total_entries;
    }
  }
  out_fh.seekp(total_size_pos);
  //cout << "num total entries: " << total_entries << endl;
  out_fh.write((char*) &total_entries, sizeof(int));
  out_fh.close();
  return;
}

void DatabaseIndex::LoadHighScoreMer(
    std::string& in_file,
    std::unordered_map<std::string, std::list<std::string> >& high_score_match
)  {
  //cout << "LoadHighScoreMer called" << endl;
  ifstream in_fh(in_file.c_str(), ios_base::in);
  if(!in_fh.good())  {
    cout << "Fatal Error: DatabaseIndex::LoadSeedExt cannot read index file " << in_file << endl;
    exit(1);
  }
  // load basic configurations: seed length, size_of_RIDType, size_of_POSType
  int loaded_mer_len, loaded_num_entries;
  in_fh.read((char*) &loaded_mer_len, sizeof(int));
  in_fh.read((char*) &loaded_num_entries, sizeof(int));
  if(loaded_num_entries <= 0)  {
    in_fh.close();
    return;
  }
  //cout << "loaded mer len:  " << loaded_mer_len << endl;
  //cout << "loaded num entires:  " << loaded_num_entries << endl;
  // load information per seed-mer
  char* src_mer = new char[loaded_mer_len + 1];
  while(!in_fh.eof())  {
    // read seed_mer  
    in_fh.read((char*) src_mer, sizeof(char) * loaded_mer_len);
    src_mer[loaded_mer_len] = '\0';
    string src_str = src_mer;
    // read num of entries
    //cout << "load src seq*****************************:  " << src_str << endl;
    int num_entries;
    in_fh.read((char*) &num_entries, sizeof(int));
    char* tgt_mer = new char[loaded_mer_len + 1];
    for(int i = 0; i < num_entries; ++ i) {
      in_fh.read((char*) tgt_mer, sizeof(char) * loaded_mer_len);
      tgt_mer[loaded_mer_len] = '\0';
      string tgt_str = tgt_mer;
      //cout << "load tgt seq:  " << tgt_str << endl;
      high_score_match[src_str].push_back(tgt_str);
      -- loaded_num_entries;
    }
    delete [] tgt_mer;
    if(loaded_num_entries <= 0)  {
      break;
    }
  }
  delete [] src_mer;
  if(loaded_num_entries < 0)  {
    cout << "Warning: DatabaseIndex::LoadReadExt corrupted file, please rebuild index " << endl;
  }
  in_fh.close();
  return;
}
