#include "reachable_reads.h"

using namespace std;

ReachableReads::ReachableReads()  {
  return;
}

ReachableReads::ReachableReads(
    std::string& reduc_map_file, std::string& seed_ext_file, 
    std::string& read_ext_file, std::string& hsm_file)  {
  DatabaseIndex db_load;
  db_load.LoadReducedMap(reduc_map_file, reduc_alph_map_);
  db_load.LoadSeedExt(seed_ext_file, seed_ext_);
  db_load.LoadReadExt(read_ext_file, fw_read_ext_, re_read_ext_);
  db_load.LoadHighScoreMer(hsm_file, high_score_match_);
  
  alph_id_ = db_load.GetAlphID();
  seed_len_ = db_load.GetSeedLen();
  overlap_len_ = db_load.GetOverlapLen();
  return;
}

ReachableReads::~ReachableReads(void) {
  return;
}

int ReachableReads::GetOverlapLen(void) {
  return overlap_len_;
}

void ReachableReads::ComputeHighScoreMatch(
    std::vector<char> alphabet, double def_high_frac,
    int mer_size, ScoringFunction<int>& score_scheme,
    ReducedAlphabet& reduc_alph
)  {
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
  double score_cutoff = def_high_frac * score_scheme.GetAveMatch() * seed_len_;
  for(int i = 0; i < (int) mer_enum.size(); ++ i) {
    for(int j = i; j < (int) mer_enum.size(); ++ j) {
      bool foo;
      if((double) score_scheme.CalMatchScore(mer_enum[i], mer_enum[j]) >= score_cutoff &&
          reduc_alph.Convert(mer_enum[i], foo) == reduc_alph.Convert(mer_enum[j], foo)
      )  {
        high_score_match_[mer_enum[i]].push_back(mer_enum[j]);
        high_score_match_[mer_enum[j]].push_back(mer_enum[i]);
      }
    }
  }
  //for(auto it = high_score_match_.begin(); it != high_score_match_.end(); ++ it) {
  //  cout << "high-score size: " << it->second.size() << endl;
  //}
  return;
}

void ReachableReads::GetHighScoreSeq(
    int seed_len, int mer_len,
    std::string& seq, std::unordered_map<std::string, bool>& hs_seq
)  {
  assert(seed_len == (int) seq.length());
  int num_iter = (int) (seed_len / mer_len);  
  int num_olp = seed_len % mer_len;
  unordered_map<string, bool> seq_temp;
  seq_temp[""] = true;
  int index = 0;
  while(index < num_iter) {
    string s = seq.substr(mer_len * index, mer_len);
    if(high_score_match_[s].size() <= 0)  {
      //cout << "Returned here" << endl;
      return;
    }
    // append sequence
    unordered_map<string, bool> seq_app;
    for(auto it_c = seq_temp.begin(); it_c != seq_temp.end(); ++ it_c) {
      for(auto it_ap = high_score_match_[s].begin(); 
          it_ap != high_score_match_[s].end(); ++ it_ap
      ) {
        seq_app[it_c->first + *it_ap] = true;
      }
    }
    seq_temp = seq_app;
    ++ index;
  }
  // handle overlap case
  if(num_olp != 0)  {
    string s = seq.substr(seq.length() - mer_len, mer_len);
    for(auto it_c = seq_temp.begin(); it_c != seq_temp.end(); ++ it_c) {
      for(auto it_ap = high_score_match_[s].begin(); 
          it_ap != high_score_match_[s].end(); ++ it_ap
      ) {
        hs_seq[it_c->first + it_ap->substr(num_olp, mer_len - num_olp)] = true;
      }
    }
  } else  {
    hs_seq = seq_temp;
  }
  return;  
}

void ReachableReads::SelectSeeds(
    double seed_score_scale, ScoringFunction<int>& score_scheme, std::string& query, 
    std::map<int, std::list<SeedType> >& candidate_seeds
) {

  // convert the query into reduced alphabet
  bool is_seq_canonical;
  ReducedAlphabet reduc_alp((Alphabet) alph_id_);
  string query_reduc = reduc_alp.Convert(query, is_seq_canonical);
  for(unsigned int i = 0; i < query_reduc.length() - seed_len_; ++ i) {
    string seed_reduc_key = query_reduc.substr(i, seed_len_);
    string seed_ori_key = query.substr(i, seed_len_);
    auto it = reduc_alph_map_.find(seed_reduc_key);
    //cout << "SelectSeeds: query sequence: " << i << " " << seed_ori_key << endl;
    if(it == reduc_alph_map_.end())  {
      continue;
    }
    //cout << "Size in the set: " << it->second.size() << endl;
    unordered_map<string, bool> hs_seq;
    GetHighScoreSeq(seed_len_, 3, seed_ori_key, hs_seq);
    //cout << "hs_seq size: " << hs_seq.size() << " " << it->second.size() << endl;
    if(hs_seq.size() <= it->second.size())  {
      //cout << "Follow first" << endl;
      for(auto it_hs = hs_seq.begin(); it_hs != hs_seq.end(); ++ it_hs) {
        string hs = it_hs->first;
        if(it->second.find(hs) != it->second.end()) {
          //cout << hs << endl;
          int match_score = score_scheme.CalMatchScore(seed_ori_key, hs);
          if(match_score > (int) (seed_score_scale * score_scheme.GetAveMatch() * seed_len_))  {
            SeedType seed;
            seed.seed_seq = hs;
            seed.q_pos = i;
            candidate_seeds[match_score].push_back(seed);
          }
        }
      }
    } else  {
      //cout << "Follow second" << endl;
      for(auto it_sd = it->second.begin(); it_sd != it->second.end(); ++ it_sd) {
        // if the alignment score between the two seed sequences is high enough
        string db_ori_key = it_sd->first;
        int match_score = score_scheme.CalMatchScore(seed_ori_key, db_ori_key);
        //cout << "SelectSeeds: target sequence:  " << db_ori_key << "  " << match_score << endl;
        if(match_score > (int) (seed_score_scale * score_scheme.GetAveMatch() * seed_len_))  {
          //cout << db_ori_key << " " << match_score << endl;
          SeedType seed;
          seed.seed_seq = db_ori_key;
          seed.q_pos = i;
          candidate_seeds[match_score].push_back(seed);
        }
      }
    }
  }
  // remove seed redundancies
  RefineSeedSameScore(ceil(0.8 * seed_len_), query.length(), candidate_seeds);
  RefineSeedDiffScore(ceil(0.8 * seed_len_), candidate_seeds);
  // DEBUG print
  //for(auto it = candidate_seeds.rbegin(); it != candidate_seeds.rend(); ++ it) {
  //  cout << it->first << endl;
  //  for(auto it_s = it->second.begin(); it_s != it->second.end(); ++ it_s) {
  //    cout << it_s->seed_seq << "," << endl;
  //  }
  //  cout << endl;
  //}
  //cout << "CheckSize SelectSeeds: " << re_read_ext_[(RIDType) 2080].size() << endl;
  return;
}

void ReachableReads::RecruitFWReads(
    int num_steps, 
    RIDType seed_rid, std::unordered_map<RIDType, bool>& fw_reads
) {
  //cout << " ++++++++++ RecruitFWReads: begin" << endl;
  if(fw_reads.find(seed_rid) != fw_reads.end())  {
    return;
  }
  queue<ReachabilityType> traversed_reads;
  ReachabilityType se;
  se.rid = seed_rid;
  se.reached_step = 0;
  traversed_reads.push(se);
  // extend the reads in the stack
  while(!traversed_reads.empty()) {
    ReachabilityType top_read = traversed_reads.front();
    traversed_reads.pop();
    //cout << " RecruitFWReads: traversed reads and step: " << top_read.rid << "  " << top_read.reached_step << endl;
    if(top_read.reached_step <= num_steps)  { // if the read can be reached
      fw_reads[top_read.rid] = true;
      //cout << " RecruitFWReads: num of extensions:  " << fw_read_ext_[top_read.rid].size() << endl;
      //cout << " RecruitFWReads: reached step:  " << top_read.reached_step << "  " << num_steps << endl;
      if(top_read.reached_step < num_steps && 
          fw_read_ext_[top_read.rid].size() < 20)  { 
        // if further extension is permitted, enqueue
        for(auto it = fw_read_ext_[top_read.rid].begin(); 
            it != fw_read_ext_[top_read.rid].end(); ++ it) {
          if(fw_reads.find(it->rid) != fw_reads.end())  {
            continue;
          }
          ReachabilityType push_read;
          push_read.rid = it->rid;
          push_read.reached_step = top_read.reached_step + 1;
          traversed_reads.push(push_read); 
        }
      }
    }
  }
  //cout << " RecruitFWReads: end ----------" << endl;
  return;
}
  
void ReachableReads::RecruitREReads(
    int num_steps, 
    RIDType seed_rid, std::unordered_map<RIDType, bool>& re_reads
) {
  // note that the seed refers to its forward sequence
  //cout << " ++++++++++ RecruitREReads: begin" << endl;
  if(re_reads.find(seed_rid) != re_reads.end())  {
    return;
  }
  queue<ReachabilityType> traversed_reads;
  ReachabilityType se;
  se.rid = seed_rid;
  se.reached_step = 0;
  traversed_reads.push(se);
  // extend the reads in the stack
  while(!traversed_reads.empty()) {
    ReachabilityType top_read = traversed_reads.front();
    traversed_reads.pop();
    //cout << " RecruitREReads: traversed reads and step: " << top_read.rid << "  " << top_read.reached_step << endl;
    if(top_read.reached_step <= num_steps)  { // if the read can be reached
      //cout << " RecruitREReads: current read accepted" << endl;
      re_reads[top_read.rid] = true;
      //cout << " RecruitREReads: num of extensions:  " << re_read_ext_[top_read.rid].size() << endl;
      //cout << " RecruitREReads: reached step:  " << top_read.reached_step << "  " << num_steps << endl;
      if(top_read.reached_step < num_steps && 
          re_read_ext_[top_read.rid].size() < 20)  { 
        // if further extension is permitted, enqueue
        for(auto it = re_read_ext_[top_read.rid].begin(); 
            it != re_read_ext_[top_read.rid].end(); ++ it) {
          //cout << " RecruitREReads: new read added: " << it->rid << endl;
          if(re_reads.find(it->rid) != re_reads.end())  {
            continue;
          }
          ReachabilityType push_read;
          push_read.rid = it->rid;
          push_read.reached_step = top_read.reached_step + 1;
          traversed_reads.push(push_read); 
        }
      }
    }
  }
  //cout << " RecruitREReads: end ----------" << endl;
  return;
}  

void ReachableReads::CollectNeighbourReads(
    int num_steps,
    std::map<int, std::list<ReadPairType> >& seed_read_pairs, 
    std::unordered_map<RIDType, bool>& candidate_reads
) {
  unordered_map<RIDType, bool> recruited;
  for(auto it = seed_read_pairs.rbegin(); it != seed_read_pairs.rend(); ++ it) {
    for(auto it_s = it->second.begin(); it_s != it->second.end(); ++ it_s) {      
      if(it_s->init_fw)  {
        unordered_map<RIDType, bool> fw_reads;
        RecruitFWReads(num_steps, it_s->rid_fw, candidate_reads);
      }
      if(it_s->init_re)  {
        unordered_map<RIDType, bool> re_reads;
        RecruitREReads(num_steps, it_s->rid_re, candidate_reads);
      }
    }
  }
  //cout << "CheckSize CollectNeighbourReads: " << re_read_ext_[(RIDType) 2080].size() << endl;
  return;
}

void ReachableReads::GetSeedReads(
    SequenceBuild& seq_obj,
    std::map<int, std::list<SeedType> >& candidate_seeds,
    std::map<int, std::list<ReadPairType> >& seed_reads
) {
  for(auto it = candidate_seeds.rbegin(); it != candidate_seeds.rend(); ++ it) {
    //cout << "***********seed score: " << it->first << endl;
    for(auto it_s = it->second.begin(); it_s != it->second.end(); ++ it_s) {
      //cout << it_s->seed_seq << endl;
      list<ReadType> fw_seed_reads;
      list<ReadType> re_seed_reads;
      set<RIDType> fw_rids;
      set<RIDType> re_rids;
      // handle forward seeds
      if(seed_ext_.find(it_s->seed_seq) != seed_ext_.end())  {
        for(auto it_p = seed_ext_[it_s->seed_seq].begin(); 
            it_p != seed_ext_[it_s->seed_seq].end(); ++ it_p
        ) {
          it_p->score = it->first;
          it_p->q_pos_fw = it_p->q_pos_re = it_s->q_pos;
          seed_reads[it->first].push_back(*it_p);
          // record the information
          //if(fw_rids.find(it_p->rid) == fw_rids.end())  {
          //  ReadType r;
          //  r.rid = it_p->rid; r.score = it->first;
          //  r.q_pos = it_s->q_pos; r.r_pos = (int) it_p->pos;
          //  fw_seed_reads.push_back(r);
          //  fw_rids.insert(r.rid);
          //}
        }
      }
      /***********************************/
      //cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
      //for(auto it = fw_seed_reads.begin(); it != fw_seed_reads.end(); ++ it) {
      //  cout << "FW:  " << seq_obj.sequence_[it->rid] << " " << it->q_pos << " " << it->r_pos << endl;
      //}
      //for(auto it = re_seed_reads.begin(); it != re_seed_reads.end(); ++ it) {
      //  cout << "RE:  " << seq_obj.sequence_[it->rid] << " " << it->q_pos << " " << it->r_pos << endl;
      //}
      //MatchReadPairs(seq_obj, fw_seed_reads, re_seed_reads, seed_reads);
    }
  }
  // DEBUG print
  //for(auto it = seed_reads.rbegin(); it != seed_reads.rend(); ++ it) {
  //  cout << "seed reads:  " << it->first << endl;
  //  for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
  //    cout << it_r->rid << "  " << seq_obj.GetHeader(it_r->rid) << endl;
  //  }
  //}
  return;
}


bool ReachableReads::IsSeqCompatible(
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


void ReachableReads::RefineLinks(
    std::unordered_map<RIDType, bool>& read_candidates,
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_connect,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_connect
)  {
  for(auto it = read_candidates.begin(); it != read_candidates.end(); ++ it) {
    // for each extension overlap
    for(auto it_rext = fw_read_ext_[it->first].begin(); 
        it_rext != fw_read_ext_[it->first].end(); ++ it_rext) {
      // if the read to be extened is also selected, record such info
      if(read_candidates.find(it_rext->rid) != read_candidates.end())  {
        //cout << "Updated local FW list" << endl;
        fw_connect[it->first].push_back(*it_rext);
      }
    }
    for(auto it_rext = re_read_ext_[it->first].begin(); 
        it_rext != re_read_ext_[it->first].end(); ++ it_rext) {
      if(read_candidates.find(it_rext->rid) != read_candidates.end())  {
        //cout << "Updated local RE list" << endl;
        re_connect[it->first].push_back(*it_rext);
      }
    }
  }
  return;
}

void ReachableReads::InitContig(
    RIDType init_read,
    std::map<RIDType, bool>& read_table, 
    std::vector<std::deque<ReadConnectType> >& pre_contigs
) {
  // insert the read to initialize a contig
  ReadConnectType r;
  r.rid = init_read;
  r.overlap_len_prev = r.overlap_len_aft = 0;
  deque<ReadConnectType> q;
  q.push_back(r);
  pre_contigs.push_back(q);
  // delete the read from the list
  read_table.erase(init_read);
  return;
}

bool ReachableReads::CheckFWTerminate(
    std::list<OverlapType> ext_fw, 
    std::unordered_map<RIDType, std::list<OverlapType> >& re_rext
)  {
  if(ext_fw.size() > 1)  {
    return true;
  } else if(ext_fw.size() == 1)  {
    // searching the reverse direction of the extension read
    auto it_rext = re_rext.find(ext_fw.begin()->rid);
    if(it_rext != re_rext.end() && it_rext->second.size() > 1)  {
      return true;
    }
  }
  return false;
}

bool ReachableReads::CheckRETerminate(
    std::list<OverlapType> ext_re, 
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext
)  {
  if(ext_re.size() > 1)  {
    return true;
  } else if(ext_re.size() == 1)  {
    // searching the reverse direction of the extension read
    auto it_rext = fw_rext.find(ext_re.begin()->rid);
    if(it_rext != fw_rext.end() && it_rext->second.size() > 1)  {
      return true;
    }
  }
  return false;
}

void ReachableReads::ConnectFW(
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_rext, 
    std::vector<std::deque<ReadConnectType> >& pre_contigs, int index,
    std::map<RIDType, bool>& read_table, 
    std::unordered_map<int, std::list<int> >& fw_link,
    std::unordered_map<int, std::list<int> >& re_link
) {
  while(true) {
    RIDType last_rid = pre_contigs[index].back().rid;
    // check if there are multiple extension options for the last read
    auto it_ext = fw_rext.find(last_rid);
    if(it_ext == fw_rext.end() || it_ext->second.size() == 0)  {
      // it means there is no extension information available, return
      break;
    } else  {
      bool is_terminate = CheckFWTerminate(it_ext->second, re_rext);
      // check if need to terminate
      if(is_terminate)  {
        // initiate each of the extension as a new contig
        for(auto it = it_ext->second.begin(); it != it_ext->second.end(); ++ it) {
          if(read_table.find(it->rid) != read_table.end())  {
            InitContig(it->rid, read_table, pre_contigs);
            // record the linking relation
            fw_link[index].push_back(pre_contigs.size() - 1);
            re_link[pre_contigs.size() - 1].push_back(index);
          }
        }
        break;
      } else  {
        if(read_table.find(it_ext->second.begin()->rid) != read_table.end())  {
          // enque the extended read
          pre_contigs[index].back().overlap_len_aft = it_ext->second.begin()->len;
          ReadConnectType rx;
          rx.rid = it_ext->second.begin()->rid;
          rx.overlap_len_prev = it_ext->second.begin()->len;
          pre_contigs[index].push_back(rx);
          // remove the read from the read table
          read_table.erase(it_ext->second.begin()->rid);
        } else  {
          break;
        }
      }
    }
  }
  return;
}

void ReachableReads::ConnectRE(
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_rext,
    std::vector<std::deque<ReadConnectType> >& pre_contigs, int index,
    std::map<RIDType, bool>& read_table,
    std::unordered_map<int, std::list<int> >& fw_link,
    std::unordered_map<int, std::list<int> >& re_link
) {
  //cout << fw_rext.size() << endl;
  //cout << re_rext.size() << endl;
  while(true) {
    RIDType first_rid = pre_contigs[index].front().rid;
    //cout << "last_rid:  " << last_rid << endl;
    // check if there are multiple extension options for the last read
    auto it_ext = re_rext.find(first_rid);
    if(it_ext == re_rext.end() || it_ext->second.size() == 0)  {
      // it means there is no extension information available, return
      //cout << "********** no extension" << endl;
      break;
    } else  {
      bool is_terminate = CheckRETerminate(it_ext->second, fw_rext);
      // check if need to terminate
      if(is_terminate)  {
        //cout << "********** multiple extension" << endl;
        // initiate each of the extension as a new contig
        for(auto it = it_ext->second.begin(); it != it_ext->second.end(); ++ it) {
          if(read_table.find(it->rid) != read_table.end())  {
            InitContig(it->rid, read_table, pre_contigs);
            // record the linking relation
            re_link[index].push_back(pre_contigs.size() - 1);
            fw_link[pre_contigs.size() - 1].push_back(index);
          }
        }
        break;
      } else  {
        //cout << "********** sinlge extension" << endl;
        if(read_table.find(it_ext->second.begin()->rid) != read_table.end())  {
          // enque the extended read
          pre_contigs[index].front().overlap_len_aft = it_ext->second.begin()->len;
          ReadConnectType rx;
          rx.rid = it_ext->second.begin()->rid;
          rx.overlap_len_aft = it_ext->second.begin()->len;
          pre_contigs[index].push_front(rx);
          // remove the read from the read table
          read_table.erase(it_ext->second.begin()->rid);
        } else  {
          break;
        }
      }
    }
  }
  return;
}

void ReachableReads::RefineSeedSameScore(
    int min_overlap, int query_len, 
    std::map<int, std::list<SeedType> >& candidate_seeds
) {
  // the function detects and removes redundancies in the seeds 
  // (requires min_overlap betwee the seeds)
  //int mid_index = (int) (query_len / 2) + 1;
  for(auto it = candidate_seeds.begin(); it != candidate_seeds.end(); ++ it) {
    auto it_s1 = it->second.begin();
    while(it_s1 != it->second.end()) {
      // get the next iterator
      auto it_s2 = std::next(it_s1);
      while(it_s2 != it->second.end()) {
        int diff = abs(it_s2->q_pos - it_s1->q_pos);
        int len_overlap = seed_len_ - diff;
        //cout << it_s1->q_pos << " " << it_s2->q_pos << "  " << diff << "  " << len_overlap << endl;
        if((len_overlap >= min_overlap) &&
          ((it_s1->q_pos > it_s2->q_pos && 
            it_s1->seed_seq.substr(0, len_overlap) == \
            it_s2->seed_seq.substr(diff, len_overlap)) ||
           (it_s1->q_pos <= it_s2->q_pos && 
            it_s1->seed_seq.substr(diff, len_overlap) == \
            it_s2->seed_seq.substr(0, len_overlap))
          )
        )  {
          // redundancy detected, 
          // get the one that has more extensions
          if(seed_ext_[it_s1->seed_seq].size() > seed_ext_[it_s2->seed_seq].size() )  {
            // if the first seed is closer to middle, delete the second one
            //cout << "Deleted: " << it_s2->seed_seq << endl; 
            auto to_delete = it_s2;
            it_s2 = std::next(it_s2);
            it->second.erase(to_delete);         
          } else  {
            // if the second seed is closer, delete the first one and break loop
            //cout << "Deleted: " << it_s1->seed_seq << endl; 
            auto to_delete = it_s1;
            it_s1 = std::next(it_s1);
            it->second.erase(to_delete);
            // move iterator back, since it will be increased again at the end
            it_s1 = std::prev(it_s1);
             
            break;
          }
        } else  {
          ++ it_s2;
        }
      }
      ++ it_s1;
    }
  }
  return;
}

void ReachableReads::RefineSeedDiffScore(
    int min_overlap, 
    std::map<int, std::list<SeedType> >& candidate_seeds
) {
  // note that the first field in the map indicates position in the query
  unordered_map<int, list<string> > presented_seeds;
  auto it = candidate_seeds.rbegin();
  // loop from high-score seeds to low-score seeds
  while(it != candidate_seeds.rend()) {
    auto it_seed = it->second.begin();
    //cout << "score tried: " << it->first << endl;
    while(it_seed != it->second.end()) {
      // check if there are reads already taken
      bool is_delete = false;
      auto to_delete = it_seed;
      int max_diff = seed_len_ - min_overlap;
      for(int i = it_seed->q_pos - max_diff; i <= it_seed->q_pos + max_diff; ++ i) {
        auto it_str_list = presented_seeds.find(i);
        if(it_str_list != presented_seeds.end())  {
          int diff = abs(i - it_seed->q_pos);
          int len_overlap = seed_len_ - diff;
          for(auto it_str = it_str_list->second.begin(); \
            it_str != it_str_list->second.end(); ++ it_str
          ) {
            //cout << it_seed->q_pos << " " << i << " " << it_seed->seed_seq << " " << *it_str << endl;
            // check if the two sequences overlap
            if((i <= it_seed->q_pos && \
                it_str->substr(diff, len_overlap) == \
                it_seed->seed_seq.substr(0, len_overlap)) || 
               (i > it_seed->q_pos && \
                it_str->substr(0, len_overlap) == \
                it_seed->seed_seq.substr(diff, len_overlap))
            )  {
              //cout << *it_str << endl;
              // redundancy detected, remove the seed
              is_delete = true; break;
            }
          }
        }
        if(is_delete)  {
          break;
        }
      }
      // erase the seed if redundancy detected
      if(is_delete)  {
        //cout << "Delete_diffscore:  " << to_delete->seed_seq << endl;
        it_seed = std::next(it_seed);
        it->second.erase(to_delete);  
      } else  {
        presented_seeds[it_seed->q_pos].push_back(it_seed->seed_seq);
        it_seed = std::next(it_seed);
      }
    }
    ++ it;
  }
  auto it_fw = candidate_seeds.begin();
  while(it_fw != candidate_seeds.end()) {
    if(it_fw->second.size() == 0)  {
      auto to_delete = it_fw;
      it_fw = std::next(it_fw);
      candidate_seeds.erase(to_delete);
    } else  {
      it_fw = std::next(it_fw);
    }
  }
  return;
}

