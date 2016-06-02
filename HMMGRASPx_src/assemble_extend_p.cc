#include "assemble_extend_p.h"

using namespace std;

void AssembleExtendP::AssembleAllContigs(
    HMMProfile &profile_obj, SequenceBuild& seq_obj,
    ReachableReads& link_obj, int band_size,
    std::map<int, std::list<SeedMatchType> >& seed_reads,
    bool dup_seeds, int assembly_depth, double p_value_cutoff,
    std::list<ContigType>& contigs
) {
  // shrink the extension links
  //unordered_map<RIDType, list<OverlapType> > fw_rext, re_rext;
  //fw_rext = link_obj.fw_read_ext_;
  //re_rext = link_obj.re_read_ext_;
  //link_obj.RefineLinks(candidate_reads, fw_rext, re_rext);
  // construct the redundancy hash table

  uint8_t *fw_alned_reads = new uint8_t [(int) seq_obj.num_seqs_ / 8 + 1];
  uint8_t *re_alned_reads = new uint8_t [(int) seq_obj.num_seqs_ / 8 + 1];
  for(int i = 0; i <= (int) seq_obj.num_seqs_ / 8; ++ i)  {
    fw_alned_reads[i] = re_alned_reads[i] = 0;
  }
  AlignInfoType *fw_aln = new AlignInfoType [(assembly_depth + 1) * (assembly_depth + 1)];
  AlignInfoType *re_aln = new AlignInfoType [(assembly_depth + 1) * (assembly_depth + 1)];
  int fw_aln_index = 0;
  int re_aln_index = 0;
  
  p_value_cutoff = p_value_cutoff >= 0.49 ? 0.49 : p_value_cutoff;
  double p_relaxed = 1.0 - 2 * p_value_cutoff;
  double score = profile_obj.CalGumbleScore(p_relaxed);
  int score_cutoff = score * SCORE_SCALE;
  // for each seed read
  unordered_map<string, bool> filtered_reads;
  for(auto it = seed_reads.begin(); it != seed_reads.end(); ++ it) {
    //cout << "SEED INFO:  " << it->first << "  " << it->second.size() << endl;
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
      //cout << it_r->seed << " " << it_r->pos << endl;
      
      if(filtered_reads.find(it_r->seed) != filtered_reads.end() || 
          link_obj.seed_ext_.find(it_r->seed) == link_obj.seed_ext_.end()
      )  {
        continue;
      }
      
      for(auto it_rext = link_obj.seed_ext_[it_r->seed].begin(); 
          it_rext != link_obj.seed_ext_[it_r->seed].end(); ++ it_rext
      ) {
        // see if the seed read is fully contained in some contigs that have already been assembled
        //if(itf != filtered_reads.end() && itf->second.find(it_rext->rid) != itf->second.end()) {
          //cout << "Seed skipped..." << endl;
        //  continue;
        //}
        
        ReadPairType read_pair = *it_rext;
        read_pair.q_pos = it_r->pos;

        //cout << "!!!!!!!!!!!!!!!  " << it_r->seed << "  " << it_r->pos << endl;
        //if(read_pair.init_fw) cout << read_pair.rid_fw << " " << seq_obj.GetHeader(read_pair.rid_fw) << endl;
        //if(read_pair.init_re) cout << read_pair.rid_re << " " << seq_obj.GetHeader(read_pair.rid_re) << endl;         

        int fw_score_cutoff = (int) (score_cutoff * (read_pair.q_pos + link_obj.seed_len_) 
            / profile_obj.GetProfileLen());
        int re_score_cutoff = score_cutoff - fw_score_cutoff;
      
        
        fw_aln_index = re_aln_index = 0;
        if(!IsReadRedundant(fw_alned_reads, read_pair.rid) && !IsReadRedundant(re_alned_reads, read_pair.rid))  {
          UpdateRedundancyRecord(fw_alned_reads, read_pair.rid);
          UpdateRedundancyRecord(re_alned_reads, read_pair.rid);
          list<ContigType> fw_contigs;
          list<ContigType> re_contigs;
          list<ContigType> merged_contigs;
          // forward extension
          SeqAlignExtendP fw_align_obj(profile_obj, band_size);
          InitFW(seq_obj, read_pair, link_obj, fw_align_obj, fw_aln, fw_aln_index);
          ExtendSeedReadFW(
              seq_obj, assembly_depth,
              fw_alned_reads, link_obj.fw_read_ext_vt_, 
              fw_align_obj, dup_seeds, fw_score_cutoff,
              fw_aln, fw_aln_index, fw_contigs
          );
          // reverse extension
          SeqAlignExtendP re_align_obj(profile_obj, band_size);
          InitRE(seq_obj, read_pair, link_obj, re_align_obj, re_aln, re_aln_index);
          ExtendSeedReadRE(
              seq_obj, assembly_depth,
              re_alned_reads, link_obj.re_read_ext_vt_, 
              re_align_obj, dup_seeds, re_score_cutoff,
              re_aln, re_aln_index, re_contigs
          );
          // filter potential seed reads
          
          FilterSeedReads(seq_obj, link_obj, fw_contigs, filtered_reads);
          FilterSeedReads(seq_obj, link_obj, re_contigs, filtered_reads);
          
          // merge contigs from both directions
          score = profile_obj.CalGumbleScore(1.0 - p_value_cutoff);
          score_cutoff = score * SCORE_SCALE;
          MergeContigs(
              profile_obj, read_pair, seq_obj, link_obj, score_cutoff, it->first,
              fw_contigs, re_contigs, contigs
          );

          
        }      
      }
    }
  }
  delete [] re_alned_reads;
  delete [] fw_alned_reads;
  delete [] fw_aln;
  delete [] re_aln;
  return;
}

void AssembleExtendP::FilterSeedReads(
    SequenceBuild& seq_obj,
    ReachableReads& link_obj,
    std::list<ContigType>& contigs,
    std::unordered_map<std::string, bool> &filtered_seeds
) {
  int seed_len = link_obj.seed_len_;
  for(auto it = contigs.begin(); it != contigs.end(); ++ it) {
    for(int i = 0; i < (int) it->sequence.length() - seed_len; ++ i) {
      string seed = it->sequence.substr(i, seed_len);
      filtered_seeds[seed] = true;
    }
  }
  return;
}

void AssembleExtendP::InitFW(
    SequenceBuild& seq_obj, ReadPairType& seed_read,
    ReachableReads& link_obj,
    SeqAlignExtendP &align_obj, AlignInfoType *fw_aln, int &fw_aln_index
) {
  // prepare read sequences
  string sr_seq = seq_obj.sequence_[seed_read.rid];
  string fw_partial_read = \
      sr_seq.substr(seed_read.r_pos, sr_seq.length() - seed_read.r_pos + 1);
  // perform inital alignment
  AlnResultStruct res_fw;
  align_obj.AlignInitFW(seed_read.q_pos, fw_partial_read, res_fw);
  // fill-up the information, prepare for extensions
  fw_aln[fw_aln_index].is_fw = true;  
  fw_aln[fw_aln_index].ext_len = fw_aln[fw_aln_index].acc_len = fw_partial_read.length();
  fw_aln[fw_aln_index].score = fw_aln[fw_aln_index].opt_score = res_fw.m_score - res_fw.n_score;
  fw_aln[fw_aln_index].opt_contig_len = fw_aln[fw_aln_index].ext_len;
  fw_aln[fw_aln_index].last_rid = fw_aln[fw_aln_index].opt_rid = seed_read.rid;
  fw_aln[fw_aln_index].mbegin = seed_read.q_pos;
  fw_aln[fw_aln_index].m_filled = res_fw.m_filled; 
  fw_aln[fw_aln_index].s_filled = res_fw.s_filled;
  fw_aln[fw_aln_index].ext_id = 0;
  fw_aln[fw_aln_index].num_extensions = 0;
  fw_aln[fw_aln_index].recorded = false;
  ++ fw_aln_index;
  return;
}

void AssembleExtendP::InitRE(
    SequenceBuild& seq_obj, ReadPairType& seed_read,
    ReachableReads& link_obj,
    SeqAlignExtendP &align_obj, AlignInfoType *re_aln, int &re_aln_index
) {
  // prepare read sequences
  string sr_seq = seq_obj.sequence_[seed_read.rid];
  string re_partial_read = \
      sr_seq.substr(0, seed_read.r_pos + link_obj.seed_len_);
  re_partial_read = string(re_partial_read.rbegin(), re_partial_read.rend());
  // perform inital alignment
  //cout << "InitRE:  " << seed_read.q_pos_re << "  " << re_partial_read << endl;
  AlnResultStruct res_re;
  align_obj.AlignInitRE(seed_read.q_pos + link_obj.seed_len_ - 1, re_partial_read, res_re);
  //cout << "InitRE:  alignment done" << endl;
  // fill-up the information, prepare for extensions
  re_aln[re_aln_index].is_fw = false;  
  re_aln[re_aln_index].ext_len = re_aln[re_aln_index].acc_len = re_partial_read.length();
  re_aln[re_aln_index].score = re_aln[re_aln_index].opt_score = res_re.m_score - res_re.n_score;
  re_aln[re_aln_index].opt_contig_len = re_aln[re_aln_index].ext_len;
  re_aln[re_aln_index].last_rid = re_aln[re_aln_index].opt_rid = seed_read.rid; 
  re_aln[re_aln_index].mbegin = seed_read.q_pos + link_obj.seed_len_ - 1;
  re_aln[re_aln_index].m_filled = res_re.m_filled; 
  re_aln[re_aln_index].s_filled = res_re.s_filled;
  re_aln[re_aln_index].ext_id = 0;
  re_aln[re_aln_index].num_extensions = 0;
  re_aln[re_aln_index].recorded = false;
  ++ re_aln_index;
  return;
}

bool cmp_info_score(const AlignInfoType& a, const AlignInfoType& b) {
  if(a.score > b.score)  {
    return true;
  }
  return false;
}

void AssembleExtendP::ExtendSeedReadFW(
    SequenceBuild& seq_obj, int assembly_depth,
    uint8_t *alned_reads,
    std::unordered_map<RIDType, std::vector<OverlapType> >& fw_rext,
    SeqAlignExtendP &fw_align_obj,
    bool dup_seeds, int score_cutoff,
    AlignInfoType *fw_aln, int &fw_aln_index,
    std::list<ContigType>& fw_contigs
) {
  //cout << " (((((((((( ExtendSeedReadFW: begin" << endl;
  // begin of loop for each possible extensions
  if(fw_aln_index <= 0)  {
    return;
  }
  RIDType source_rid = fw_aln[fw_aln_index - 1].last_rid;
  while(fw_aln_index > 0) {
    //cout << "alignment record processed" << endl;
    // get the top alignment to process
    
    /*
    cout << "<******" << endl;
    cout << "NUM EXTENSIONS:  " << fw_aln.back().num_extensions << endl;
    cout << "CONTIG:  " << fw_aln.back().contig << endl;
    cout << "FILLED REGION: " << fw_aln.back().m_filled << "  " << fw_aln.back().s_filled << endl;
    cout << "LAST READ: " << fw_aln.back().last_rid << endl;
    cout << "POSSIBLE EXT SIZE: " << fw_rext[fw_aln.back().last_rid].size() << endl;
    if(fw_aln.back().ext_ptr != fw_rext[fw_aln.back().last_rid].end())  {
      cout << "TO EXTD: " << fw_aln.back().ext_ptr->rid << endl; 
    }
    cout << "******>" << endl;
    */
    auto ext_ptr = fw_rext.find(fw_aln[fw_aln_index - 1].last_rid);
    if(fw_aln[fw_aln_index - 1].num_extensions > assembly_depth || 
       ext_ptr == fw_rext.end() ||
       fw_aln[fw_aln_index - 1].ext_id == (int) (ext_ptr->second).size() ||
       IsReadRedundant(alned_reads, (ext_ptr->second)[fw_aln[fw_aln_index - 1].ext_id].rid) ||
       fw_aln[fw_aln_index - 1].s_filled - fw_aln[fw_aln_index - 1].m_filled > fw_align_obj.band_size_)  {
      
      /*
      cout << "termination criterion met" << endl;
      if(fw_aln.back().num_extensions > assembly_depth)  {
        cout << "extension depth" << endl;
      }
      if(fw_aln.back().ext_ptr == fw_rext[fw_aln.back().last_rid].end())  {
        cout << "no extension path" << endl;
      }
      if(IsReadRedundant(alned_reads, fw_aln.back().last_rid))  {
        cout << "redundant read" << endl;
      }
      */
      
      // termination criteria, no need to extend the current path
      if(!fw_aln[fw_aln_index - 1].recorded)  {
        int recorded_len = 0;
        for(int k = fw_aln_index - 1; k >= 0; -- k)  {
          // if the alignment is recorded, then all alignments before it should be recorded
          // no need to go further, stop and record the maximal length that has been recorded
          if(fw_aln[k].recorded)  {
            recorded_len = (int) fw_aln[k].acc_len;
            break;
          } else if(!dup_seeds || (int) fw_aln[k].acc_len <= fw_aln[fw_aln_index - 1].opt_contig_len)  {
          // otherwise if the alignment is less than contig length to be record,
          // then we know the alignment should be recorded, mark it
            fw_aln[k].recorded = true;
          // record that the current extension read has been used
            UpdateRedundancyRecord(alned_reads, fw_aln[k].last_rid);
          }
        }
        // if the current path contains some part that has not been recorded
        if(recorded_len < fw_aln[fw_aln_index - 1].opt_contig_len 
            && fw_aln[fw_aln_index - 1].opt_score < score_cutoff
        )  {
          ContigType c;
          //c.sequence = fw_aln.back().contig.substr(0, fw_aln.back().opt_contig_len);
          c.sequence = "";
          for(int k = 0; k < fw_aln_index; ++ k)  {
            // get the forward sequences
            string s = seq_obj.sequence_[fw_aln[k].last_rid];
            c.sequence += s.substr(s.length() - fw_aln[k].ext_len, fw_aln[k].ext_len);
          }
          //c.score = fw_aln.back().opt_score;
          c.score = fw_aln[fw_aln_index - 1].score;
          //c.fw_rid = fw_aln.back().opt_rid;
          c.fw_rid = fw_aln[fw_aln_index - 1].last_rid;
          c.re_rid = source_rid;
          fw_contigs.push_back(c);
          //cout << "++++++ FW +++++" << endl;
          //cout << c.sequence << endl;
          //cout << seq_obj.sequence_[c.fw_rid] << endl;
        }
      }
      // remove the alignment record from the stack
      -- fw_aln_index;
      // go next iteration to check the top record
      continue;
    }
    // otherwise we need to perform alignment to extend the record
    OverlapType olp = (ext_ptr->second)[fw_aln[fw_aln_index - 1].ext_id];
    string sr = seq_obj.sequence_[olp.rid];
    string ext_seq = sr.substr(
        (int) olp.len, sr.length() - (int) olp.len
    );
    // modify the current alignment record
    RIDType last_rid = olp.rid;
    ++ fw_aln[fw_aln_index - 1].ext_id;
    // run alignment
    AlnResultStruct res_fw;
    fw_align_obj.AlignExtdFW(
        fw_aln[fw_aln_index - 1].mbegin,    
        fw_aln[fw_aln_index - 1].m_filled, fw_aln[fw_aln_index - 1].s_filled,
        ext_seq, res_fw
    );
    // based on the results record information for next extension
    fw_aln[fw_aln_index].is_fw = true;
    fw_aln[fw_aln_index].ext_len = ext_seq.length();
    fw_aln[fw_aln_index].acc_len = fw_aln[fw_aln_index - 1].acc_len + fw_aln[fw_aln_index].ext_len;
    fw_aln[fw_aln_index].score = res_fw.m_score - res_fw.n_score;
    fw_aln[fw_aln_index].m_filled = res_fw.m_filled; 
    fw_aln[fw_aln_index].s_filled = res_fw.s_filled;
    //cout << "FW results: " << next_ext.score << " " << next_ext.contig << endl;
    // check if the extension is desirable
    if(fw_aln[fw_aln_index].score > fw_aln[fw_aln_index - 1].opt_score + 5 * SCORE_SCALE)  {
      // the extension is not desirable, do not record it
      //cout << "!!! Undesirable extension... quit... " << next_ext.score << "  " << fw_aln.back().score << endl;
      continue;
    }
    if(fw_aln[fw_aln_index].score < fw_aln[fw_aln_index - 1].opt_score)  {  // note that less score is better
      fw_aln[fw_aln_index].opt_score = fw_aln[fw_aln_index].score;
      fw_aln[fw_aln_index].opt_contig_len = fw_aln[fw_aln_index].acc_len;
      fw_aln[fw_aln_index].opt_rid = last_rid;
    } else  {
      fw_aln[fw_aln_index].opt_score = fw_aln[fw_aln_index - 1].opt_score;
      fw_aln[fw_aln_index].opt_contig_len = fw_aln[fw_aln_index - 1].opt_contig_len;
      fw_aln[fw_aln_index].opt_rid = fw_aln[fw_aln_index - 1].opt_rid;
    }
    fw_aln[fw_aln_index].mbegin = fw_aln[fw_aln_index - 1].mbegin;
    
    fw_aln[fw_aln_index].last_rid = last_rid;
    fw_aln[fw_aln_index].ext_id = 0;
    fw_aln[fw_aln_index].num_extensions = fw_aln[fw_aln_index - 1].num_extensions + 1;   
    fw_aln[fw_aln_index].recorded = false;   
    // push new alignment record into the stack
    ++ fw_aln_index;
    
    //cout << "FW results opt: " << next_ext.opt_score << endl;
  }
  //cout << " ExtendSeedReadFW: end ))))))))))" << endl;
  return;
}

void AssembleExtendP::ExtendSeedReadRE(
    SequenceBuild& seq_obj, int assembly_depth,
    uint8_t *alned_reads,
    std::unordered_map<RIDType, std::vector<OverlapType> >& re_rext,
    SeqAlignExtendP &re_align_obj,
    bool dup_seeds, int score_cutoff,
    AlignInfoType *re_aln, int &re_aln_index,
    std::list<ContigType>& re_contigs
) {
  //cout << " (((((((((( ExtendSeedReadRE: begin" << endl;
  if(re_aln_index <= 0)  {
    return;
  }
  RIDType source_rid = re_aln[re_aln_index - 1].last_rid;
  while(re_aln_index > 0) {
    /*
    cout << "<******" << endl;
    cout << "NUM EXTENSIONS:  " << re_aln.back().num_extensions << endl;
    cout << "CONTIG:  " << re_aln.back().contig << endl;
    cout << "FILLED REGION: " << re_aln.back().m_filled << "  " << re_aln.back().s_filled << endl;
    cout << "LAST READ: " << re_aln.back().last_rid << endl;
    cout << "POSSIBLE EXT SIZE: " << re_rext[re_aln.back().last_rid].size() << endl;
    if(re_aln.back().ext_ptr != re_rext[re_aln.back().last_rid].end())  {
      cout << "TO EXTD: " << re_aln.back().ext_ptr->rid << endl; 
    }
    cout << "******>" << endl;
    */
    auto ext_ptr = re_rext.find(re_aln[re_aln_index - 1].last_rid);
    // get the top alignment to process
    if(re_aln[re_aln_index - 1].num_extensions > assembly_depth || 
       ext_ptr == re_rext.end() ||
       re_aln[re_aln_index - 1].ext_id == (int) (ext_ptr->second).size() ||
       IsReadRedundant(alned_reads, (ext_ptr->second)[re_aln[re_aln_index - 1].ext_id].rid) ||
       re_aln[re_aln_index - 1].s_filled - re_aln[re_aln_index - 1].m_filled > re_align_obj.band_size_)  {
      
      /*
      cout << "termination criterion met" << endl;
      if(re_aln.back().num_extensions > assembly_depth)  {
        cout << "extension depth" << endl;
      }
      if(re_aln.back().ext_ptr == re_rext[re_aln.back().last_rid].end())  {
        cout << "no extension path" << endl;
      }
      if(IsReadRedundant(alned_reads, re_aln.back().last_rid))  {
        cout << "redundant read" << endl;
      }
      */
      
      // termination criteria, no need to extend the current path
      if(!re_aln[re_aln_index - 1].recorded)  {
        //cout << "&&& TRYING BACK TRACK" << endl;
        int recorded_len = 0;
        for(int k = re_aln_index - 1; k >= 0; -- k)  {
          // if the alignment is recorded, then all alignments before it should be recorded
          // no need to go further, stop and record the maximal length that has been recorded
          //cout << it->contig << " " << (int) it->recorded << endl;
          if(re_aln[k].recorded)  {
            recorded_len = (int) re_aln[k].acc_len;
            break;
          } else if(!dup_seeds || (int) re_aln[k].acc_len <= re_aln[re_aln_index - 1].opt_contig_len)  {
          // otherwise if the alignment is less than contig length to be record,
          // then we know the alignment should be recorded, mark it
            re_aln[k].recorded = true;
          // record that the current extension read has been used
            UpdateRedundancyRecord(alned_reads, re_aln[k].last_rid);
          }
        }
        //cout << "&&&: " << recorded_len << "  " << re_aln.back().opt_contig_len << endl;
        // if the current path contains some part that has not been recorded
        if(recorded_len < re_aln[re_aln_index - 1].opt_contig_len 
            && re_aln[re_aln_index - 1].opt_score < score_cutoff
        )  {
          //cout << "&&& PATH RECORD BEGIN" << endl;
          ContigType c;
          //c.sequence = re_aln.back().contig.substr(0, re_aln.back().opt_contig_len);
          c.sequence = ""; 
          for(int k = 0; k < re_aln_index; ++ k)  {
            string s = seq_obj.sequence_[re_aln[k].last_rid];
            s = s.substr(0, re_aln[k].ext_len);
            c.sequence += string(s.rbegin(), s.rend());
          }
          //c.score = re_aln.back().opt_score;
          c.score = re_aln[re_aln_index - 1].score;
          //c.re_rid = re_aln.back().opt_rid;
          c.re_rid = re_aln[re_aln_index - 1].last_rid;
          c.fw_rid = source_rid;
          re_contigs.push_back(c);
          //cout << "++++++ RE +++++" << endl;
          //cout << c.sequence << endl;
          //cout << seq_obj.sequence_[c.fw_rid] << endl;
        }
      }
      // remove the alignment record from the stack
      -- re_aln_index;
      // go next iteration to check the top record
      continue;
    }
    // otherwise we need to perform alignment to extend the record
    OverlapType olp = (ext_ptr->second)[re_aln[re_aln_index - 1].ext_id];
    string sr = seq_obj.sequence_[olp.rid];
    string ext_seq = sr.substr(
        0, sr.length() - (int) olp.len
    );
    ext_seq = string(ext_seq.rbegin(), ext_seq.rend());
    // modify the current alignment record
    RIDType last_rid = olp.rid;
    ++ re_aln[re_aln_index - 1].ext_id;
    // run alignment
    AlnResultStruct res_re;
    re_align_obj.AlignExtdRE(
        re_aln[re_aln_index - 1].mbegin,    
        re_aln[re_aln_index - 1].m_filled, re_aln[re_aln_index - 1].s_filled,
        ext_seq, res_re
    );
    // based on the results record information for next extension
    re_aln[re_aln_index].is_fw = false;  
    re_aln[re_aln_index].ext_len = ext_seq.length();
    re_aln[re_aln_index].acc_len = re_aln[re_aln_index - 1].acc_len + re_aln[re_aln_index].ext_len;
    re_aln[re_aln_index].score = res_re.m_score - res_re.n_score;
    re_aln[re_aln_index].m_filled = res_re.m_filled; 
    re_aln[re_aln_index].s_filled = res_re.s_filled;
    //cout << "RE results: " << next_ext.score << " " << next_ext.contig << endl;
    if(re_aln[re_aln_index].score > re_aln[re_aln_index - 1].opt_score + 5 * SCORE_SCALE)  {
      continue;
    }
    if(re_aln[re_aln_index].score < re_aln[re_aln_index - 1].opt_score)  {  // note that less score is better
      re_aln[re_aln_index].opt_score = re_aln[re_aln_index].score;
      re_aln[re_aln_index].opt_contig_len = re_aln[re_aln_index].acc_len;
      re_aln[re_aln_index].opt_rid = last_rid;
    } else  {
      re_aln[re_aln_index].opt_score = re_aln[re_aln_index - 1].opt_score;
      re_aln[re_aln_index].opt_contig_len = re_aln[re_aln_index - 1].opt_contig_len;
      re_aln[re_aln_index].opt_rid = re_aln[re_aln_index - 1].opt_rid;
    }
    re_aln[re_aln_index].mbegin = re_aln[re_aln_index - 1].mbegin;
    re_aln[re_aln_index].last_rid = last_rid;
    re_aln[re_aln_index].ext_id = 0;
    re_aln[re_aln_index].num_extensions = re_aln[re_aln_index - 1].num_extensions + 1;  
    re_aln[re_aln_index].recorded = false;     
    ++ re_aln_index;
    
    //cout << "RE results opt: " << next_ext.opt_score << endl;
  }
  //cout << " ExtendSeedReadRE: end ))))))))))" << endl;
  return;
}

void AssembleExtendP::MergeContigs(
    HMMProfile &q_profile, ReadPairType &seed_pair,
    SequenceBuild& seq_obj, ReachableReads& link_obj,
    int score_cutoff, int seed_score,
    std::list<ContigType>& fw_contigs, std::list<ContigType>& re_contigs,
    std::list<ContigType>& merged_contigs
) {
  //cout << "******************************************" << endl;
  //cout << "e_value cutoff:  " << e_value_cutoff << endl;
  //cout << "score cutoff:  " << score_cutoff << endl;
  // find the partial(forward or reverse) contig that has the best alignment score
  ContigType best_fw_contig, best_re_contig;
  int best_fw_score = INF, best_re_score = INF;
  auto best_fw_iter = fw_contigs.begin();
  auto best_re_iter = re_contigs.begin();
  for(auto it = fw_contigs.begin(); it != fw_contigs.end(); ++ it) {
    if(it->score < best_fw_score)  {
      best_fw_score = it->score;
      best_fw_contig = *it;
      best_fw_iter = it;
    }
  }
  for(auto it = re_contigs.begin(); it != re_contigs.end(); ++ it) {
    if(it->score < best_re_score)  {
      best_re_score = it->score;
      best_re_contig = *it;
      best_re_iter = it;
    }
  }
  //cout << "MergeContigs:  scores: " << best_fw_score << " " << best_re_score << endl;
  //cout << "MergeContigs:  seed_score: " << seed_score << endl;
  // merge the contigs
  if(best_fw_score <= seed_score && best_re_score <= seed_score)  {
    // take the best fw contig and attach to all re contigs
    //cout << "MERGE BOTH" << endl;
    int fw_cl = best_fw_contig.sequence.length() - link_obj.seed_len_;
    string fw_app_seq = best_fw_contig.sequence.substr(link_obj.seed_len_, fw_cl);
    for(auto it = re_contigs.begin(); it != re_contigs.end(); ++ it) {
      if(it == best_re_iter)  continue;
      ContigType c;
      c.score = it->score + best_fw_contig.score - seed_score;
      c.fw_rid = best_fw_contig.fw_rid;
      c.re_rid = it->re_rid;
      //cout << "MergeContigs: contig score:  " << c.score << " " << it->score << " " << best_fw_contig.score << " " << seed_score << endl;
      if(c.score <= score_cutoff)  {
        c.sequence = string(it->sequence.rbegin(), it->sequence.rend()) + fw_app_seq;
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
      if(it == best_fw_iter)  continue;
      ContigType c;
      c.score = best_re_contig.score + it->score - seed_score;
      c.fw_rid = it->fw_rid;
      c.re_rid = best_re_contig.re_rid;
      //cout << "MergeContigs: contig score:  " << c.score << " " << it->score << " " << best_re_contig.score << " " << seed_score << endl;
      if(c.score <= score_cutoff)  {
        c.sequence = re_app_seq + it->sequence;
        merged_contigs.push_back(c);    
        //if(c.q_begin > c.q_end)  {   
        //  cout << "MergeContigs (both):  " << c.q_begin << "  " << c.q_end << endl;
        //} 
        //cout << "MergeContigs (both):  " << c.sequence << endl;
      }
    }
    ContigType c;
    c.score = best_re_contig.score + best_fw_contig.score - seed_score;
    if(c.score <= score_cutoff)  {
      c.sequence = re_app_seq + best_fw_contig.sequence;
      c.fw_rid = best_fw_contig.fw_rid;
      c.re_rid = best_re_contig.re_rid;
      merged_contigs.push_back(c);     
      //cout << "MergeContigs (both):  " << c.sequence << endl;
    }
  } else if(best_fw_score <= seed_score && best_re_score > seed_score) {
    // directly take all fw contigs
    //cout << "FORWARD" << endl;
    string s = seq_obj.sequence_[seed_pair.rid];
    //cout << "SSSEEEQQQ: " << s << endl;
    string seed_tipp = s.substr(0, seed_pair.r_pos);
    for(auto it = fw_contigs.begin(); it != fw_contigs.end(); ++ it) {
      if(it->score <= score_cutoff)  {
        it->sequence = seed_tipp + it->sequence;
        merged_contigs.push_back(*it);     
        //cout << "MergeContigs (fw):  " << it->sequence << endl;
      }
    }
  } else if(best_fw_score > seed_score && best_re_score <= seed_score) {
    // directly take all re contigs
    //cout << "REVERSE" << endl;
    string s = seq_obj.sequence_[seed_pair.rid];
    int n = seed_pair.r_pos + link_obj.GetSeedLen();
    string seed_tipp = s.substr(n, s.length() - n);
    for(auto it = re_contigs.begin(); it != re_contigs.end(); ++ it) {
      if(it->score <= score_cutoff)  {
        it->sequence = string(it->sequence.rbegin(), it->sequence.rend()) + seed_tipp;
        merged_contigs.push_back(*it);     
        //cout << "MergeContigs (re):  " << it->sequence << endl;
      }
    }
  }
  //cout << "end of mergecontigs" << endl;
  return;
}


void AssembleExtendP::ProgressiveExtension(
    SequenceBuild& seq_obj,
    ReachableReads& link_obj,
    int n,
    std::list<ContigType>& contigs
) {
  for(auto itc = contigs.begin(); itc != contigs.end(); ++ itc) {
    unordered_map<RIDType, bool> inc_reads;
    // extend to the N terminus
    RIDType pivot_re = itc->re_rid;
    inc_reads[pivot_re] = true;
    //cout << "@@@@@@@@@@@@@@@@@@@@@@ begin of RE extension: " << endl;
    int fw_n = n, re_n = n;   // number of extension steps
    while(fw_n > 0 && link_obj.re_read_ext_[pivot_re].size() == 1) {
      OverlapType o_info = link_obj.re_read_ext_[pivot_re].back();
      if(inc_reads.find(o_info.rid) != inc_reads.end())  {
        break;
      } else  {
        string ext_seq = seq_obj.sequence_[o_info.rid];
        //cout << "!!!pivot  " << pivot_re << endl;
        //cout << "!!!ext " << o_info.rid << endl;
        //cout << "!!!pivotseq  " << seq_obj.sequence_[pivot_re] << endl;
        //cout << "!!!extseq  " << ext_seq << endl;
        //cout << "!!!olen  " << (unsigned int) o_info.len << endl;
        //cout << "!!!contig_seq  " << itc->sequence << endl;
        itc->sequence = ext_seq.substr(0, ext_seq.length() - o_info.len) + itc->sequence;
        //cout << "###  " << itc->sequence << endl;
        inc_reads[o_info.rid] = true;
        pivot_re = o_info.rid;
      }
      -- fw_n;
    }
    //cout << "end of RE extension" << endl;
    // extend to the C terminus
    //cout << "@@@@@@@@@@@@@@@@@@@@@@ begin of FW extension: " << endl;
    RIDType pivot_fw = itc->fw_rid;
    inc_reads[pivot_fw] = true;
    while(re_n > 0 && link_obj.fw_read_ext_[pivot_fw].size() == 1) {
      OverlapType o_info = link_obj.fw_read_ext_[pivot_fw].back();
      if(inc_reads.find(o_info.rid) != inc_reads.end())  {
        break;
      } else  {
        string ext_seq = seq_obj.sequence_[o_info.rid];
        //cout << "!!!pivot  " << pivot_fw << endl;
        //cout << "!!!ext " << o_info.rid << endl;
        //cout << "!!!pivotseq  " << seq_obj.sequence_[pivot_fw] << endl;
        //cout << "!!!extseq  " << ext_seq << endl;
        //cout << "!!!olen  " << (unsigned int) o_info.len << endl;
        //cout << "!!!contig_seq  " << itc->sequence << endl;
        itc->sequence = itc->sequence + ext_seq.substr(o_info.len, ext_seq.length() - o_info.len);
        //cout << "###  " << itc->sequence << endl;
        inc_reads[o_info.rid] = true;
        pivot_fw = o_info.rid;
      }
      -- re_n;
    }
    //cout << "end of FW extension" << endl;
  }
  return;
}
