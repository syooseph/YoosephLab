#include "assemble_extend.h"

using namespace std;

AssembleExtend::AssembleExtend(void)  {
  return;
}

AssembleExtend::~AssembleExtend(void) {
  return;
}

std::string AssembleExtend::EncodeHeader(
    std::string& barcode, ContigType& contig
) {
  string out_header = barcode + "::";
  std::stringstream outs;
  outs << contig.q_begin << "-" << contig.q_end << "::";
  outs << contig.e_value << "::" << contig.bit_score << "::" << contig.score;
  out_header += outs.str();
  return out_header;
}

bool AssembleExtend::DecodeHeader(
    std::string& header, std::string& barcode, ContigType& contig
) {
  return true;
}

void AssembleExtend::AssembleAllContigs(
    std::string& query, SequenceBuild& seq_obj,
    ReachableReads& link_obj, 
    ScoringFunction<int>& score_obj, int band_size,
    std::map<int, std::list<ReadPairType> >& seed_reads,
    int assembly_depth,
    double e_value_cutoff,
    std::list<ContigType>& contigs
) {
  // shrink the extension links
  //unordered_map<RIDType, list<OverlapType> > fw_rext, re_rext;
  //fw_rext = link_obj.fw_read_ext_;
  //re_rext = link_obj.re_read_ext_;
  //link_obj.RefineLinks(candidate_reads, fw_rext, re_rext);
  // construct the redundancy hash table
  unordered_map<RIDType, bool> fw_alned_reads;
  unordered_map<RIDType, bool> re_alned_reads;
  // for each seed read
  for(auto it = seed_reads.rbegin(); it != seed_reads.rend(); ++ it) {
    //cout << "size:  " << it->second.size() << endl;
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
      //ReadPosType seed_read_pos;
      //seed_read_pos.rid = it_r->rid;
      //seed_read_pos.begin = it_r->q_pos - it_r->r_pos + 1;
      //seed_read_pos.end = seed_read_pos.begin + strlen(seq_obj.sequence_[it_r->rid]) - 1;
      
      //cout << "Seed read ID:  " << it_r->rid << endl;
      //cout << "Seed read name:  " << seq_obj.GetHeader(it_r->rid) << endl;
      //cout << "Seed read sequence:  " << sseq << endl;
      //cout << "position " << it_r->r_pos_fw << endl;
      //cout << "////////////////////////////" << endl << endl;
      //string sseq = seq_obj.sequence_[it_r->rid_fw];
      //cout << "Seed sequence: " << sseq.substr(it_r->r_pos_fw, 6) << endl;
      //if(it_r->init_fw && it_r->init_re)  {
      //  cout << "EXTENSION_BOTH" << endl;
      //  cout << "Seed read name:  " << seq_obj.GetHeader(it_r->rid_fw) << " " << seq_obj.GetHeader(it_r->rid_re) << endl;
      //  cout << "Positions: " << it_r->r_pos_fw << "  " << it_r->r_pos_re << endl;
      //} else if(!it_r->init_fw && it_r->init_re) {
      //  cout << "EXTENSION_RE" << endl;
      //  cout << "Seed read name:  " << seq_obj.GetHeader(it_r->rid_re) << endl;
      //  cout << "Positions: " << it_r->r_pos_re << endl;
      //} else if(it_r->init_fw && !it_r->init_re) {
      //  cout << "EXTENSION_FW" << endl;
      //  cout << "Seed read name:  " << seq_obj.GetHeader(it_r->rid_fw) << endl;
      //  cout << "Positions: " << it_r->r_pos_fw << endl;
      //}
      
      stack<AlignInfoType> fw_aln, re_aln;
      list<ContigType> fw_contigs, re_contigs;
      if(it_r->init_fw && !IsReadRedundant(fw_alned_reads, it_r->rid_fw))  {
        UpdateRedundancyRecord(fw_alned_reads, it_r->rid_fw);
        InitFW(query, seq_obj, *it_r, score_obj, band_size, link_obj, fw_aln);
        ExtendSeedReadFW(
            query, seq_obj, score_obj, band_size, assembly_depth,
            fw_alned_reads, link_obj.fw_read_ext_, fw_aln, fw_contigs
        );
      }        
      if(it_r->init_re && !IsReadRedundant(re_alned_reads, it_r->rid_re))  {
        UpdateRedundancyRecord(re_alned_reads, it_r->rid_re);
        InitRE(query, seq_obj, *it_r, score_obj, band_size, link_obj, re_aln);
        ExtendSeedReadRE(
            query, seq_obj, score_obj, band_size, assembly_depth,
            re_alned_reads, link_obj.re_read_ext_, re_aln, re_contigs
        );
      }
      //for(auto itp = fw_contigs.begin(); itp != fw_contigs.end(); ++ itp) {
      //  cout << "New contig fw:  " << itp->score << " " << itp->sequence << endl;
      //}
      //for(auto itp = re_contigs.begin(); itp != re_contigs.end(); ++ itp) {
      //  cout << "New contig re:  " << itp->score << " " << string(itp->sequence.rbegin(), itp->sequence.rend()) << endl;
      //}
      // merge results
      MergeContigs(
          query, seq_obj, link_obj, score_obj, 
          e_value_cutoff, it->first,
          fw_contigs, re_contigs, contigs
      );
      //for(auto itp = contigs.begin(); itp != contigs.end(); ++ itp) {
      //  cout << "New contig merged:  " << itp->score << " " << itp->sequence << endl;
      //}
      
    }
  }
  return;
}

void AssembleExtend::MergeContigs(
    std::string& query, SequenceBuild& seq_obj, ReachableReads& link_obj,
    ScoringFunction<int>& score_obj, 
    double e_value_cutoff, int seed_score,
    std::list<ContigType>& fw_contigs, std::list<ContigType>& re_contigs,
    std::list<ContigType>& merged_contigs
) {
  // compute raw score cutoff
  long int m = 1 + (long int) query.length();
  long int n = 1 + (long int) seq_obj.GetDBSizeInMegaBase();
  n = n * 1000000;
  int score_cutoff = score_obj.ComputeRawScore(
      m, n, e_value_cutoff
  );
  //cout << "******************************************" << endl;
  //cout << "e_value cutoff:  " << e_value_cutoff << endl;
  //cout << "score cutoff:  " << score_cutoff << endl;
  // find the partial(forward or reverse) contig that has the best alignment score
  ContigType best_fw_contig, best_re_contig;
  int best_fw_score = 0, best_re_score = 0;
  auto best_fw_iter = fw_contigs.begin();
  auto best_re_iter = re_contigs.begin();
  for(auto it = fw_contigs.begin(); it != fw_contigs.end(); ++ it) {
    if(it->score > best_fw_score)  {
      best_fw_score = it->score;
      best_fw_contig = *it;
      best_fw_iter = it;
    }
  }
  for(auto it = re_contigs.begin(); it != re_contigs.end(); ++ it) {
    if(it->score > best_re_score)  {
      best_re_score = it->score;
      best_re_contig = *it;
      best_re_iter = it;
    }
  }
  //cout << "MergeContigs:  scores: " << best_fw_score << " " << best_re_score << endl;
  //cout << "MergeContigs:  seed_score: " << seed_score << endl;
  // merge the contigs
  if(best_fw_score > seed_score && best_re_score > seed_score)  {
    // take the best fw contig and attach to all re contigs
    int fw_cl = best_fw_contig.sequence.length() - link_obj.seed_len_;
    string fw_app_seq = best_fw_contig.sequence.substr(link_obj.seed_len_, fw_cl);
    for(auto it = re_contigs.begin(); it != re_contigs.end(); ++ it) {
      if(it == best_re_iter)  {
        continue;
      }
      ContigType c;
      c.score = it->score + best_fw_contig.score - seed_score;
      //cout << "MergeContigs: contig score:  " << c.score << " " << it->score << " " << best_fw_contig.score << " " << seed_score << endl;
      if(c.score >= score_cutoff)  {
        c.bit_score = score_obj.ComputeBitScore(c.score);
        c.sequence = string(it->sequence.rbegin(), it->sequence.rend()) + fw_app_seq;
        c.q_begin = it->q_begin;
        c.q_end = best_fw_iter->q_end;
        merged_contigs.push_back(c);  
        //if(c.q_begin > c.q_end)  {   
        //  cout << "MergeContigs (both):  " << c.q_begin << "  " << c.q_end << endl;
        //}
      }
    }
    // take the best re contig and attach to all fw_contigs
    int re_cl = best_re_contig.sequence.length() - link_obj.seed_len_;
    // reverse the contig sequence
    string re_app_seq = string(best_re_contig.sequence.rbegin(), best_re_contig.sequence.rend());
    re_app_seq = re_app_seq.substr(0, re_cl);
    for(auto it = fw_contigs.begin(); it != fw_contigs.end(); ++ it) {
      if(it == best_fw_iter)  {
        continue;
      }
      ContigType c;
      c.score = best_re_contig.score + it->score - seed_score;
      //cout << "MergeContigs: contig score:  " << c.score << " " << it->score << " " << best_re_contig.score << " " << seed_score << endl;
      if(c.score >= score_cutoff)  {
        c.bit_score = score_obj.ComputeBitScore(c.score);
        c.sequence = re_app_seq + it->sequence;
        c.q_begin = best_re_iter->q_begin;
        c.q_end = it->q_end;
        merged_contigs.push_back(c);    
        //if(c.q_begin > c.q_end)  {   
        //  cout << "MergeContigs (both):  " << c.q_begin << "  " << c.q_end << endl;
        //} 
        //cout << "MergeContigs (both):  " << c.sequence << endl;
      }
    }
    ContigType c;
    c.score = best_re_contig.score + best_fw_contig.score - seed_score;
    if(c.score >= score_cutoff)  {
      c.bit_score = score_obj.ComputeBitScore(c.score);
      c.sequence = re_app_seq + best_fw_contig.sequence;
      c.q_begin = best_re_iter->q_begin;
      c.q_end = best_fw_iter->q_end;
      merged_contigs.push_back(c);     
      //cout << "MergeContigs (both):  " << c.sequence << endl;
    }
  } else if(best_fw_score > seed_score && best_re_score <= seed_score) {
    // directly take all fw contigs
    for(auto it = fw_contigs.begin(); it != fw_contigs.end(); ++ it) {
      if(it->score >= score_cutoff)  {
        merged_contigs.push_back(*it);     
        //cout << "MergeContigs (fw):  " << it->sequence << endl;
      }
    }
  } else if(best_fw_score <= seed_score && best_re_score > seed_score) {
    // directly take all re contigs
    for(auto it = re_contigs.begin(); it != re_contigs.end(); ++ it) {
      if(it->score >= score_cutoff)  {
        it->sequence = string(it->sequence.rbegin(), it->sequence.rend());
        merged_contigs.push_back(*it);     
        //cout << "MergeContigs (re):  " << it->sequence << endl;
      }
    }
  }
  //cout << "end of mergecontigs" << endl;
  return;
}

double AssembleExtend::ComputeJaccardIndex(int i1, int j1, int i2, int j2) {
  double jc_index = 0.0;
  if(i1 > j1 || i2 > j2)  {
    cout << "Warning: AssembleExtend::ComputeJacardIndex: undefined \
      region input:" << i1 << "," << j1 << "-" << i2 << "," << j2 << endl;
  }
  if(j1 < i2 || j2 < i1)  {
    // no overlap
    jc_index = 0.0;
  } else if(i2 >= i1 && j2 <= j1) {
    // fully contained
    jc_index = (double) (j2 - i2 + 1) / (j1 - i1 + 1);
  } else if(i1 >= i2 && j1 <= j2) {
    // fully contained
    jc_index = (double) (j1 - i1 + 1) / (j2 - i2 + 1);
  } else if(i1 <= i2 && i2 <= j1 && j1 <= j2) {
    // partial overlap
    jc_index = (double) (j1 - i2 + 1) / (j2 - i1 + 1);
  } else if(i2 <= i1 && i1 <= j2 && j2 <= j1) {
    // partial overlap
    jc_index = (double) (j2 - i1 + 1) / (j1 - i2 + 1);
  }
  return jc_index;
}

bool AssembleExtend::IsReadRedundant(
    std::unordered_map<RIDType, bool>& alned_reads, 
    RIDType read_id
) {
  auto it = alned_reads.find(read_id);
  if(it == alned_reads.end())  {
    // if the read is not contained then there is no redundancy
    return false;
  } else  {
    return true;
    // need to check which region the read has aligned to
    //for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
    //  if(
    //    ComputeJaccardIndex(it_r->begin, it_r->end, read_pos.begin, read_pos.end) >= jc_index
    //  )  {
    //    return true;
    //  }
    //}
  }
  return true;
}

void AssembleExtend::UpdateRedundancyRecord(
    std::unordered_map<RIDType, bool>& alned_reads, 
    RIDType read_id
) {
  alned_reads[read_id] = true;
  /*****************
  RangeType r_range;
  r_range.begin = read_pos.begin; r_range.end = read_pos.end;
  auto it = alned_reads.find(read_pos.rid);
  if(it == alned_reads.end())  {
    // if the read is not contained then there is no redundancy
    alned_reads[read_pos.rid].push_back(r_range);
  } else  {
    bool updated = false;
    // need to check which region the read has aligned to
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
      if(ComputeJaccardIndex(it_r->begin, it_r->end, read_pos.begin, read_pos.end) \
          >= jc_index
      )  {
        // expand the range
        it_r->begin = it_r->begin < read_pos.begin ? it_r->begin : read_pos.begin;
        it_r->end = it_r->end > read_pos.end ? it_r->end : read_pos.end;
        updated = true;
      }
    }
    if(!updated)  {
      it->second.push_back(r_range);
    }
  }
  ****************/
  return;
}

std::string AssembleExtend::FetchQuerySeqFW(
    std::string& query, int q_begin, int len
) {
  if(q_begin < 0)  {
    cout << "Warning: AssembleExtend::FetchQuerySeqFW: query begin \
      out of range: " << q_begin << endl;
    return "";
  } else if(len <= 0 || q_begin >= (int) query.length()) {
    return "";
  }
  // check range
  int q_end = q_begin + len - 1;
  q_end = q_end < (int) query.length() ? q_end : (int) query.length() - 1;
  // return sequence
  return query.substr(q_begin, q_end - q_begin + 1);
}

std::string AssembleExtend::FetchQuerySeqRE(
    std::string& query, int q_end, int len
) {
  if(q_end >= (int) query.length())  {
    cout << "Warning: AssembleExtend::FetchQuerySeqRE: query end \
      out of range: " << q_end << endl;
    return "";
  } else if(len <= 0 || q_end < 0) {
    return "";
  }
  // check range
  int q_begin = q_end - len + 1;
  q_begin = q_begin >= 0 ? q_begin : 0;
  // return sequence
  return query.substr(q_begin, q_end - q_begin + 1);
}


void AssembleExtend::InitFW(
    std::string& query, SequenceBuild& seq_obj, ReadPairType& seed_read, 
    ScoringFunction<int>& score_obj, int band_size,
    ReachableReads& link_obj,
    std::stack<AlignInfoType>& fw_aln
) {
  
  // prepare sequences
  string sr_seq = seq_obj.sequence_[seed_read.rid_fw];
  //cout << "InitFW sr_seq: " << sr_seq << "  " << seed_read.r_pos << endl;
  string fw_partial_read = \
      sr_seq.substr(seed_read.r_pos_fw, sr_seq.length() - seed_read.r_pos_fw + 1);
  string fw_partial_query = FetchQuerySeqFW(\
    query, seed_read.q_pos_fw, fw_partial_read.length() + band_size - 1);
  //cout << "InitFW query_seq: " << fw_partial_query << endl;
  //cout << "InitFw target_seq:  " << fw_partial_read << endl;
  // do alignment and initialize forward extension
  SeqAlign<int> aln_init_fw(
      fw_partial_query, fw_partial_read, &score_obj, GLOBAL,
      band_size, band_size, false
  );
  aln_init_fw.Align();
  AlignInfoType fw_aln_info;
  fw_aln_info.contig = fw_partial_read;
  aln_init_fw.FillEdgeScore(fw_aln_info.edge_score);
  //cout << "InitFW edge score dimension: " << fw_aln_info.edge_score.begin_row << "  " << fw_aln_info.edge_score.begin_column << endl;
  fw_aln_info.is_fw = true;
  fw_aln_info.score = fw_aln_info.opt_score = aln_init_fw.GetBestScore();
  //cout << "initialize fw alignment score:  " << fw_aln_info.opt_score << endl;
  fw_aln_info.bit_score = fw_aln_info.opt_bit_score = \
      score_obj.ComputeBitScore(fw_aln_info.score);
  fw_aln_info.prev_bit_score = 0;
  fw_aln_info.opt_contig_len = fw_aln_info.contig.length();
  fw_aln_info.q_begin = seed_read.q_pos_fw;
  fw_aln_info.q_end = fw_aln_info.q_begin + fw_partial_query.length() - 1;
  fw_aln_info.q_end = fw_aln_info.q_end > (int) query.length() - 1 ? \
      (int) query.length() - 1 : fw_aln_info.q_end;
  fw_aln_info.opt_q_begin = fw_aln_info.q_begin;
  fw_aln_info.opt_q_end = fw_aln_info.q_end;
  //cout << "InitFW q_range:  " << fw_aln_info.q_begin << " " << fw_aln_info.q_end << endl;
  //cout << "InitFW q_seq:  " << fw_partial_query << endl;
  fw_aln_info.last_rid = seed_read.rid_fw;
  fw_aln_info.num_extensions = 0;
  fw_aln.push(fw_aln_info);
  
  return;
}

void AssembleExtend::InitRE(
    std::string& query, SequenceBuild& seq_obj, ReadPairType& seed_read, 
    ScoringFunction<int>& score_obj, int band_size,
    ReachableReads& link_obj,
    std::stack<AlignInfoType>& re_aln
) {
  
  // prepare sequences
  string sr_seq = seq_obj.sequence_[seed_read.rid_re];
  string re_partial_read = \
      sr_seq.substr(0, seed_read.r_pos_re + link_obj.seed_len_);
  string re_partial_query = FetchQuerySeqRE(\
    query, seed_read.q_pos_re + link_obj.seed_len_ - 1,\
    re_partial_read.length() + band_size - 1);
  re_partial_read = string(re_partial_read.rbegin(), re_partial_read.rend());
  re_partial_query = string(re_partial_query.rbegin(), re_partial_query.rend());
  //cout << "query_seq: " << re_partial_query << endl;
  // do alignment and initialize reverse extension
  //cout << "InitRE query_seq: " << re_partial_query << endl;
  //cout << "InitRE target_seq:  " << re_partial_read << endl;
  SeqAlign<int> aln_init_re(
      re_partial_query, re_partial_read, &score_obj, GLOBAL,
      band_size, band_size, false
  );
  aln_init_re.Align();
  AlignInfoType re_aln_info;
  re_aln_info.contig = re_partial_read;
  aln_init_re.FillEdgeScore(re_aln_info.edge_score);
  re_aln_info.is_fw = false;
  re_aln_info.score = re_aln_info.opt_score = aln_init_re.GetBestScore();
  //cout << "initialize re alignment score:  " << re_aln_info.opt_score << endl;
  re_aln_info.bit_score = re_aln_info.opt_bit_score = \
      score_obj.ComputeBitScore(re_aln_info.score);
  re_aln_info.prev_bit_score = 0;
  re_aln_info.opt_contig_len = re_aln_info.contig.length();
  re_aln_info.q_end = seed_read.q_pos_re + link_obj.seed_len_ - 1;
  re_aln_info.q_begin = re_aln_info.q_end - re_partial_query.length() + 1;
  re_aln_info.q_begin = re_aln_info.q_begin >= 0 ? re_aln_info.q_begin : 0;
  re_aln_info.opt_q_begin = re_aln_info.q_begin;
  re_aln_info.opt_q_end = re_aln_info.q_end;
  //cout << "InitRE:  " << re_aln_info.q_begin << " " << re_aln_info.q_end << endl;
  //cout << "InitRE q_range:  " << re_aln_info.q_begin << " " << re_aln_info.q_end << endl;
  //cout << "InitRE q_seq:  " << re_partial_query << endl;
  re_aln_info.last_rid = seed_read.rid_re;
  re_aln_info.num_extensions = 0;
  re_aln.push(re_aln_info);
  
  return;
}

bool cmp_info_score(const AlignInfoType& a, const AlignInfoType& b) {
  if(a.score > b.score)  {
    return true;
  }
  return false;
}

void AssembleExtend::ExtendSeedReadFW(
    std::string& query, SequenceBuild& seq_obj,
    ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
    std::unordered_map<RIDType, bool>& alned_reads,
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext, 
    std::stack<AlignInfoType>& fw_aln, std::list<ContigType>& fw_contigs
) {
  //cout << " (((((((((( ExtendSeedReadFW: begin" << endl;
  list<AlignInfoType> terminus;
  // get the top alignment to process
  while(!fw_aln.empty()) {
    //cout << "Seed extension:  " << endl;
    AlignInfoType cr = fw_aln.top();
    fw_aln.pop();
    list<AlignInfoType> extended_alns;
    // do alignment on all valid branches
    //cout << "In seed extension alignment size: " << cr.edge_score.begin_row << "  " << cr.edge_score.begin_column  << endl;
    ExtendFWPhase(
        query, seq_obj, score_obj, band_size, assembly_depth, 
        alned_reads, fw_rext, cr, extended_alns,
        terminus
    );
    if(extended_alns.empty())  {
      // terminate, record current path as a compelte contig
      terminus.push_back(cr);
    } else  {
      // if more extension possible, sort and put the high-score one on top
      //cout << "Extension found" << endl;
      extended_alns.sort(cmp_info_score);
      for(auto it = extended_alns.rbegin(); it != extended_alns.rend(); ++ it) {
        fw_aln.push(*it);
      }
    }
  }
  for(auto it = terminus.begin(); it != terminus.end(); ++ it) {
    // convert the terminus information into contig
    ContigType c;
    c.sequence = it->contig.substr(0, it->opt_contig_len);
    c.score = it->opt_score;
    c.bit_score = it->opt_bit_score;
    c.q_begin = it->opt_q_begin;
    c.q_end = it->opt_q_end;
    fw_contigs.push_back(c);
  }
  //cout << " ExtendSeedReadFW: end ))))))))))" << endl;
  return;
}

void AssembleExtend::ExtendFWPhase(
    std::string& query, SequenceBuild& seq_obj,
    ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
    std::unordered_map<RIDType, bool>& alned_reads,
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_rext, 
    AlignInfoType& src_aln, std::list<AlignInfoType>& emitted_alns,
    std::list<AlignInfoType>& terminus
)  {
  //cout << " [[[[[[[[[[ ExtendFWPhase: begin" << endl;
  if(!src_aln.is_fw)  {
    cout << "Warning: AssembleExtend::ExtendPhase: Forward extension \
      applied on reverse alignment" << endl;
    return;
  }
  auto it = fw_rext.find(src_aln.last_rid);
  if(it != fw_rext.end())  {
    string seed_contig_seq = src_aln.contig;
    string seed_query_seq = query.substr(
        src_aln.q_begin, src_aln.q_end - src_aln.q_begin + 1
    );
    //cout << "ExtendFWPhase: seed sequences:  " << endl << seed_query_seq << endl << seed_contig_seq << endl;
    //cout << "ExtendFWPhase: before for loop" << endl;
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
      //cout << "ExtendFWPhase: extending read:  " << it_r->rid << "  " << seq_obj.GetHeader(it_r->rid) << endl;
      string sr = seq_obj.sequence_[it_r->rid];
      // check read redundancy
      if(IsReadRedundant(alned_reads, it_r->rid))  {
        //cout << "ExtendFWPhase: Redundancy detected" << endl;
        continue;
      }
      // construct sequences for alignment      
      string extd_contig_seq = sr.substr((int) it_r->len, sr.length() - (int) it_r->len);
      string extd_query_seq = FetchQuerySeqFW(
          query, src_aln.q_begin + seed_query_seq.length(), 
          extd_contig_seq.length()
      );
      //cout << "ExtendFWPhase: extd sequences:  " << endl << extd_query_seq << endl << extd_contig_seq << endl;
      //cout << "ExtendFWPhase: after getting extend query sequence" << endl;
      if((int) seed_query_seq.length() + (int) extd_query_seq.length() <= \
          (int) seed_contig_seq.length() + (int) extd_contig_seq.length() - band_size)  {
        // extension reaches the end of the query
        //cout << "ExtendFWPhase: Extension reaches end" << endl;
        continue;
      }
      // do extend alignment
      AlignInfoType emit_aln;
      bool is_good_ext = DoExtendAlnFW(
          query, it_r->rid,
          seed_query_seq, seed_contig_seq, extd_query_seq, extd_contig_seq,
          score_obj, band_size, assembly_depth, src_aln, emit_aln    
      );
      //cout << " ExtendFWPhase: alignment result: " << emit_aln.contig << " " << emit_aln.score << " " << emit_aln.opt_score << endl;
      //cout << " ExtendFWPhase: bit_score: " << emit_aln.bit_score << " " << emit_aln.opt_bit_score << endl;
      if(is_good_ext)  {
        //cout << "ExtendFWPhase: good alignment inserted" << endl;
        emitted_alns.push_back(emit_aln);
        UpdateRedundancyRecord(alned_reads, it_r->rid);
      }
    }
  }
  //cout << " ExtendFWPhase: end ]]]]]]]]]]" << endl;
  return;
}

bool AssembleExtend::DoExtendAlnFW(
    std::string& query, RIDType ext_rid,
    std::string& seed_query, std::string& seed_contig,
    std::string& extd_query, std::string& extd_contig, 
    ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
    AlignInfoType& src_aln, AlignInfoType& emit_aln
)  {
  // do extend alignment
  SeqAlignExtend<int> aln_ext(
      seed_query, seed_contig, extd_query, extd_contig,
      src_aln.edge_score, &score_obj, GLOBAL, 
      band_size, band_size, false
  );
  aln_ext.Align();
  emit_aln.score = aln_ext.GetBestScore();
  emit_aln.bit_score = score_obj.ComputeBitScore(emit_aln.score);
  emit_aln.prev_bit_score = src_aln.bit_score;
  aln_ext.FillEdgeScore(emit_aln.edge_score);
  //cout << "DoExtendAlnFW: edge begins next: " << emit_aln.edge_score.begin_row << " " << emit_aln.edge_score.begin_column << endl;
  emit_aln.contig = seed_contig + extd_contig;
  emit_aln.q_begin = src_aln.q_begin;
  emit_aln.q_end = src_aln.q_end + extd_query.length();
  emit_aln.q_end = emit_aln.q_end > (int) query.length() - 1 ? \
      (int) query.length() - 1 : emit_aln.q_end;
  //cout << "DoExtendAlnFW q_range:  " << emit_aln.q_begin << " " << emit_aln.q_end << endl;
  //cout << "DoExtendAlnFW q_seq:  " << seed_query + extd_query << endl;
  emit_aln.is_fw = true;
  emit_aln.last_rid = ext_rid;
  emit_aln.num_extensions = src_aln.num_extensions + 1;
  
  // check if the current alignment is better than before, if yes, record it
  //cout << "DoExtendAlnFW: src opt_bit_score:  " << src_aln.opt_bit_score << endl;
  //cout << "DoExtendAlnFW: emit opt_bit_score:  " << emit_aln.opt_bit_score << endl;
  if(emit_aln.score > src_aln.opt_score)  {
    emit_aln.opt_bit_score = emit_aln.bit_score;
    emit_aln.opt_score = emit_aln.score;
    //cout << "DoExtendAlnFW: opt_score:  " << emit_aln.opt_score << endl;
    emit_aln.opt_contig_len = seed_contig.length() + extd_contig.length();
    emit_aln.opt_q_begin = emit_aln.q_begin;
    emit_aln.opt_q_end = emit_aln.q_end;
  } else  {
    emit_aln.opt_bit_score = src_aln.opt_bit_score;
    emit_aln.opt_score = src_aln.opt_score;
    emit_aln.opt_contig_len = src_aln.opt_contig_len;
    emit_aln.opt_q_begin = src_aln.opt_q_begin;
    emit_aln.opt_q_end = src_aln.opt_q_end;
  }
  // check score drop-off
  if(emit_aln.bit_score < emit_aln.opt_bit_score - 25 || 
      emit_aln.num_extensions > assembly_depth)  {
    return false;
  }
  return true;
}

void AssembleExtend::ExtendSeedReadRE(
    std::string& query, SequenceBuild& seq_obj,
    ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
    std::unordered_map<RIDType, bool>& alned_reads,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_rext, 
    std::stack<AlignInfoType>& re_aln, std::list<ContigType>& re_contigs
) {
  //cout << " (((((((((( ExtendSeedReadRE: begin" << endl;
  // get the top alignment to process
  list<AlignInfoType> terminus;
  while(!re_aln.empty()) {
    AlignInfoType cr = re_aln.top();
    re_aln.pop();
    list<AlignInfoType> extended_alns;
    // do alignment on all valid branches
    ExtendREPhase(
        query, seq_obj, score_obj, band_size, assembly_depth,
        alned_reads, re_rext, cr, extended_alns,
        terminus
    );
    if(extended_alns.empty())  {
      // terminate, record current path as a compelte contig
      terminus.push_back(cr);
    } else  {
      // if more extension possible, sort and put the high-score one on top
      //cout << "Extension found" << endl;
      extended_alns.sort(cmp_info_score);
      for(auto it = extended_alns.rbegin(); it != extended_alns.rend(); ++ it) {
        re_aln.push(*it);
      }
    }
  }
  for(auto it = terminus.begin(); it != terminus.end(); ++ it) {
    ContigType c;
    c.sequence = it->contig.substr(0, it->opt_contig_len);
    c.score = it->opt_score;
    c.bit_score = it->opt_bit_score;
    c.q_begin = it->opt_q_begin;
    c.q_end = it->opt_q_end;
    re_contigs.push_back(c);
  }
  //cout << " ExtendSeedReadRE: end ))))))))))" << endl;
  return;
}

void AssembleExtend::ExtendREPhase(
    std::string& query, SequenceBuild& seq_obj,
    ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
    std::unordered_map<RIDType, bool>& alned_reads,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_rext, 
    AlignInfoType& src_aln, std::list<AlignInfoType>& emitted_alns,
    std::list<AlignInfoType>& terminus
) {
  //cout << " [[[[[[[[[[ ExtendREPhase: begin" << endl;
  if(src_aln.is_fw)  {
    cout << "Warning: AssembleExtend::ExtendPhase: Reverse extension \
      applied on forward alignment" << endl;
    return;
  }
  //cout << " ExtendREPhase: src alignment extension read id: " << src_aln.last_rid << endl;
  //cout << "CheckSize: " << re_rext[(RIDType) 2080].size() << endl;
  auto it = re_rext.find(src_aln.last_rid);
  if(it != re_rext.end())  {
    string seed_contig_seq = src_aln.contig;
    //cout << "ExtendREPhase: before getting query" << endl;
    //cout << "ExtendREPhase: " << src_aln.q_begin << " " << src_aln.q_end - src_aln.q_begin + 1 << endl;
    string seed_query_seq = query.substr(
        src_aln.q_begin, src_aln.q_end - src_aln.q_begin + 1
    );
    seed_query_seq = string(seed_query_seq.rbegin(), seed_query_seq.rend());
    //cout << "ExtendREPhase: seed sequences:  " << endl << seed_query_seq << endl << seed_contig_seq << endl;
    //cout << "ExtendREPhase: before for loop" << endl;
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
      //cout << "ExtendREPhase: extending read:  " << it_r->rid << "  " << seq_obj.GetHeader(it_r->rid) << endl;
      string sr = seq_obj.sequence_[it_r->rid];
      // check read redundancy
      if(IsReadRedundant(alned_reads, it_r->rid))  {
        continue;
      }
      // construct sequences for alignment
      string r_frag = sr.substr(0, sr.length() - (int) it_r->len);      
      string extd_contig_seq = string(r_frag.rbegin(), r_frag.rend());
      //cout << "ExtendREPhase: before getting extend query sequence" << endl;
      string extd_query_seq = FetchQuerySeqRE(
          query, src_aln.q_end - seed_query_seq.length(), 
          extd_contig_seq.length()
      );
      //cout << "ExtendREPhase: extd sequences:  " << endl << extd_query_seq << endl << extd_contig_seq << endl;
      //cout << "ExtendREPhase: after getting extend query sequence" << endl;
      extd_query_seq = string(extd_query_seq.rbegin(), extd_query_seq.rend());
      if((int) seed_query_seq.length() + (int) extd_query_seq.length() <= \
          (int) seed_contig_seq.length() + (int) extd_contig_seq.length() - band_size)  {
        // extension reaches the end of the query
        //cout << " ExtendREPhase: Exit for reaching end of query sequence" << endl;
        //cout << " ExtendREPhase: lengths: " << seed_query_seq.length() << " " << extd_query_seq.length() << " " << seed_contig_seq.length() << " " << extd_contig_seq.length() << endl;
        //cout << " ExtendREPhase:  " << extd_query_seq << endl;
        //cout << " ExtendREPhase:  " << extd_contig_seq << endl;
        continue;
      }
      // do extend alignment
      AlignInfoType emit_aln;
      bool is_good_ext = DoExtendAlnRE(
          query, it_r->rid,
          seed_query_seq, seed_contig_seq, extd_query_seq, extd_contig_seq,
          score_obj, band_size, assembly_depth, src_aln, emit_aln    
      );
      //cout << " ExtendREPhase: alignment result: " << emit_aln.contig << " " << emit_aln.score << " " << emit_aln.opt_score << endl;
      //cout << " ExtendREPhase: bit_score: " << emit_aln.bit_score << " " << emit_aln.opt_bit_score << endl;
      if(is_good_ext)  {
        //cout << " ExtendREPhase: good alignment inserted" << endl;
        emitted_alns.push_back(emit_aln);
        UpdateRedundancyRecord(alned_reads, it_r->rid);
      }
    }
  }
  //cout << " ExtendREPhase: end ]]]]]]]]]]" << endl;
  return;
}

bool AssembleExtend::DoExtendAlnRE(
    std::string& query, RIDType ext_rid,
    std::string& seed_query, std::string& seed_contig,
    std::string& extd_query, std::string& extd_contig, 
    ScoringFunction<int>& score_obj, int band_size, int assembly_depth,
    AlignInfoType& src_aln, AlignInfoType& emit_aln
)  {
  //cout << "DoExtendAlnRE: begin" << endl;
  // do extend alignment
  SeqAlignExtend<int> aln_ext(
      seed_query, seed_contig, extd_query, extd_contig,
      src_aln.edge_score, &score_obj, GLOBAL, 
      band_size, band_size, false
  );
  aln_ext.Align();
  // compute alignment score and bit score
  emit_aln.score = aln_ext.GetBestScore();
  emit_aln.bit_score = score_obj.ComputeBitScore(emit_aln.score);
  emit_aln.prev_bit_score = src_aln.bit_score;
  aln_ext.FillEdgeScore(emit_aln.edge_score);
  emit_aln.contig = seed_contig + extd_contig;
  emit_aln.q_end = src_aln.q_end;
  emit_aln.q_begin = src_aln.q_begin - extd_query.length();
  emit_aln.q_begin = emit_aln.q_begin >= 0 ? emit_aln.q_begin : 0;
  //cout << "DoExtendAlnRE: " << src_aln.q_begin << " " << src_aln.q_end << "  " << emit_aln.q_begin << "  " << emit_aln.q_end << endl;
  //cout << "DoExtendAlnRE q_range:  " << emit_aln.q_begin << " " << emit_aln.q_end << endl;
  //cout << "DoExtendAlnRE q_seq:  " << seed_query + extd_query << endl;
  emit_aln.is_fw = false;
  emit_aln.last_rid = ext_rid;
  emit_aln.num_extensions = src_aln.num_extensions + 1;
  // check if the current alignment is better than before, if yes, record it
  if(emit_aln.score > src_aln.opt_score)  {
    emit_aln.opt_bit_score = emit_aln.bit_score;
    emit_aln.opt_score = emit_aln.score;
    emit_aln.opt_contig_len = seed_contig.length() + extd_contig.length();
    emit_aln.opt_q_begin = emit_aln.q_begin;
    emit_aln.opt_q_end = emit_aln.q_end;
  } else  {
    emit_aln.opt_bit_score = src_aln.opt_bit_score;
    emit_aln.opt_score = src_aln.opt_score;
    emit_aln.opt_contig_len = src_aln.opt_contig_len;
    emit_aln.opt_q_begin = src_aln.opt_q_begin;
    emit_aln.opt_q_end = src_aln.opt_q_end;
  }
  // check score drop-off
  if(emit_aln.bit_score < emit_aln.opt_bit_score - 25 || 
      emit_aln.num_extensions > assembly_depth)  {
    return false;
  }
  
  //cout << "DoExtendAlnRE: end" << endl;
  return true;
}

