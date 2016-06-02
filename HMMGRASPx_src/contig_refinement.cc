#include "contig_refinement.h"

using namespace std;

ContigRefinement::ContigRefinement()  {
  return;
}

ContigRefinement::~ContigRefinement() {
  return;
}

bool _cmp_contig(ContigType& c1, ContigType& c2)  {
  if(c1.score > c2.score)  {
    return true;
  }
  return false;
}

bool _cmp_merpos_fw(MerPosType& m1, MerPosType& m2) {
  if(m1.fw_len > m2.fw_len)  {
    return true;
  }
  return false;
}

bool _cmp_merpos_re(MerPosType& m1, MerPosType& m2) {
  if(m1.re_len > m2.re_len)  {
    return true;
  }
  return false;
}

void ContigRefinement::RefineContigs(
    std::string& query, SequenceBuild& seq_obj, ScoringFunction<int>& score_obj, 
    int band_size, double e_value, int mer_len, std::list<ContigType>& contigs, 
    std::list<ContigType>& refined_contigs
) {
  if(contigs.size() <= 0)  {
    return;
  }
  //cout << "RefineContigs called" << endl;
  // sort the contigs based on alignment score
  contigs.sort(_cmp_contig);
  contig_holder_.resize(contigs.size());
  int i = 0;
  for(auto it = contigs.begin(); it != contigs.end(); ++ it) {
    //cout << "contig:  " << it->sequence << endl;
    it->valid = true;
    contig_holder_[i] = *it;
    ++ i;
  }
  // build hash table for the mers
  unordered_map<string, list<MerPosType> > mer_contig;
  IndexContigs(mer_len, mer_contig);
  // merge the contigs based on overlap
  MergeAllContigs(
      query, seq_obj, score_obj, band_size, e_value, 
      mer_len, mer_contig, refined_contigs
  );
  refined_contigs.sort(_cmp_contig);
  //cout << "RefineContigs end" << endl;
  return;
}

void ContigRefinement::RefineContigsNoAln(
    SequenceBuild &seq_obj, int mer_len, 
    std::list<ContigType> &contigs, std::list<ContigType> &refined_contigs
) {
  if(contigs.size() <= 0)  {
    return;
  }
  contigs.sort(_cmp_contig);
  contig_holder_.resize(contigs.size());
  int i = 0;
  for(auto it = contigs.rbegin(); it != contigs.rend(); ++ it) {
    //cout << "contig:  " << it->sequence << endl;
    it->valid = true;
    contig_holder_[i] = *it;
    ++ i;
  }
  unordered_map<string, list<MerPosType> > mer_contig;
  IndexContigs(mer_len, mer_contig);
  // greedily select longest contig
  for(i = 0; i < (int) contig_holder_.size(); ++ i) {
    if(contig_holder_[i].valid)  {
      ContigType c = contig_holder_[i];
      list<int> incorporated_contigs;
      IncorporateContig(false, mer_len, mer_contig, i, c, incorporated_contigs);
      refined_contigs.push_back(c); 
      for(auto it_inc = incorporated_contigs.begin(); 
        it_inc != incorporated_contigs.end(); ++ it_inc) {
        contig_holder_[*it_inc].valid = false;
      }
      contig_holder_[i].valid = false;
    }
  }
  refined_contigs.sort(_cmp_contig);
  return;
}

void ContigRefinement::IndexContigs(
    int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& mer_contig 
) {
  //cout << "IndexContigs called" << endl;
  // count all mers in contigs
  for(int i = 0; i < (int) contig_holder_.size(); ++ i) {
    //cout << i << endl;
    for(int j = 0; j < (int) contig_holder_[i].sequence.length() - mer_len; ++ j) {
      //cout << contig_holder_[i].sequence << endl;
      string s = contig_holder_[i].sequence;
      string mer = s.substr(j, mer_len);
      MerPosType mp;
      mp.cid = i; mp.pos = j; 
      mp.fw_len = s.length() - j; mp.re_len = j + mer_len; 
      mer_contig[mer].push_back(mp);
    }
  }
  //cout << "IndexContigs end" << endl;
  return;
}

void ContigRefinement::MergeAllContigs(
    std::string query, SequenceBuild& seq_obj, ScoringFunction<int>& score_obj, 
    int band_size, double e_value,
    int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& mer_contig,
    std::list<ContigType>& refined_contigs
) {
  // greedy merge: start with contig with highest score, merge with longest overlap
  int score_cutoff = score_obj.ComputeRawScore(
      query.length(), seq_obj.GetDBSizeInMegaBase() * 1000000, e_value
  );
  for(int i = 0; i < (int) contig_holder_.size(); ++ i) {
    if(contig_holder_[i].valid)  {
      ContigType c = contig_holder_[i];
      list<int> incorporated_contigs;
      IncorporateContig(true, mer_len, mer_contig, i, c, incorporated_contigs);
      // redo alignment for the contig
      string q_par = query.substr(c.q_begin, c.q_end - c.q_begin + 1);
      //cout << "******************" << endl;
      //cout << "q length:  " << query.length() << endl;
      //cout << "index: " << c.q_begin << " " << c.q_end << endl; 
      //cout << "MergeAllContigs: ref:  " << c.sequence << endl;
      //cout << "MergeAllContigs: tgt:  " << q_par << endl;
      SeqAlign<int> realn(
          q_par, c.sequence, &score_obj, 
          LOCAL, band_size * 2, band_size * 2, false
      );
      realn.Align();
      realn.TraceBack();
      realn.GetAlignment(
          c.al_print.nuc_match, 
          c.al_print.seq1, c.al_print.seq2, c.al_print.symbol
      );
      c.score = realn.GetBestScore();
      if(c.score >= score_cutoff)  {
        
        c.bit_score = score_obj.ComputeBitScore(c.score);
        c.e_value = score_obj.ComputeEValue(
            query.length(), seq_obj.GetDBSizeInMegaBase() * 1000000,
            c.score
        );
        //cout << c.score << "  " << c.e_value << " " << query.length() << "  " << seq_obj.GetDBSizeInMegaBase() * 1000000 << endl;
        
        refined_contigs.push_back(c); 
        for(auto it_inc = incorporated_contigs.begin(); 
            it_inc != incorporated_contigs.end(); ++ it_inc) {
          contig_holder_[*it_inc].valid = false;
        }
        contig_holder_[i].valid = false;
      }
    }
  }
  return;
}

bool _cmp_contig_overlap(ContigOverlapType& c1, ContigOverlapType& c2)  {
  if(c1.count >= c2.count)  {
    return true;
  }
  return false;
}

void ContigRefinement::IncorporateContig(
    bool handle_index, int mer_len, 
    std::unordered_map<std::string, std::list<MerPosType> >& kmer_map,
    int contig_index, ContigType& current_contig, std::list<int>& incorporated_contigs
)  {
  //cout << " IncorporateContig: begin" << endl;
  // count all contigs that share k-mers with the current contig
  unordered_map<int, ContigOverlapType> similar_contigs;
  for(int i = 0; i < (int) current_contig.sequence.length() - mer_len; ++ i) {
    string mer = current_contig.sequence.substr(i, mer_len);
    for(auto it = kmer_map[mer].begin(); it != kmer_map[mer].end(); ++ it) {
      // find if the contig exists
      if(it->cid == contig_index || !contig_holder_[it->cid].valid)  {
        continue;
      }
      auto it_contig = similar_contigs.find(it->cid);
      if(it_contig != similar_contigs.end())  {
        // update the information
        it_contig->second.count ++;
        it_contig->second.ref_end = i;
        it_contig->second.target_end = it->pos;
      } else  {
        // insert the record
        ContigOverlapType cl_record;
        cl_record.cid = it->cid;
        cl_record.count = 1;
        cl_record.ref_begin = cl_record.ref_end = i;
        cl_record.target_begin = cl_record.target_end = it->pos;
        similar_contigs.insert({it->cid, cl_record});
      }
    }
  }
  // filter and sort the overlap contigs
  list<ContigOverlapType> filtered_contigs;
  for(auto it = similar_contigs.begin(); it != similar_contigs.end(); ++ it) {
    if(it->second.ref_end - it->second.ref_begin + 1 == it->second.count &&
       it->second.target_end - it->second.target_begin + 1 == it->second.count)  {
      filtered_contigs.push_back(it->second);
    }
  }
  filtered_contigs.sort(_cmp_contig_overlap);
  // try to merge the contig
  int ref_index = 0;
  for(auto it = filtered_contigs.begin(); it != filtered_contigs.end(); ++ it) {
    //cout << " IncorporateContig:  target contig:" << endl;
    //cout << contig_holder_[it->cid].sequence << endl;
    if(TryMergeOverlapedContigs(handle_index, mer_len, *it, current_contig, ref_index))  {
      // mark the target contig as invalid
      incorporated_contigs.push_back(it->cid);
    }
    //cout << " IncorporateContig:  merged contig:" << endl;
    //cout << current_contig.sequence << endl;
    //cout << " IncorporateContig:  recruited contig ids: " << endl;
    //for(auto it_cid = incorporated_contigs.begin(); it_cid != incorporated_contigs.end(); ++ it_cid) {
    //  cout << " @@: " << *it_cid << endl; 
    //}
  }
  return;
}

bool ContigRefinement::TryMergeOverlapedContigs(
    bool handle_index, int mer_len, ContigOverlapType& overlap_record, 
    ContigType& ref_contig, int& ref_index
) {
  //cout << "**********************************" << endl;
  //cout << "   TryMergeOverlapContigs: begin" << endl;
  //cout << "   TryMergeOverlapContigs: ref_index:  " << ref_index << endl;
  ContigType target_contig = contig_holder_[overlap_record.cid];
  //cout << " TryMergeOverlapContigs: ref sequence: " << ref_contig.sequence << endl;
  //cout << " TryMergeOverlapContigs: tgt sequence: " << target_contig.sequence << endl;
  //cout << " TryMergeOverlapContigs: ref ends: " << ref_contig.q_begin << "  " << ref_contig.q_end << endl;
  //cout << " TryMergeOverlapContigs: tgt ends: " << target_contig.q_begin << "  " << target_contig.q_end << endl;
  //cout << " TryMergeOverlapContigs: overlap positions:  " << overlap_record.ref_begin << "  " << overlap_record.ref_end << "  " << overlap_record.target_begin << " " << overlap_record.target_end << endl;
  //cout << "   TryMergeOverlapContigs: overlap seq1" << endl;
  //cout << ref_contig.sequence.substr(ref_index + overlap_record.ref_begin, overlap_record.ref_end - overlap_record.ref_begin + mer_len + 1) << endl;
  //cout << "   TryMergeOverlapContigs: overlap seq2" << endl;
  //cout << target_contig.sequence.substr(overlap_record.target_begin, overlap_record.target_end - overlap_record.target_begin + mer_len + 1) << endl;
  //cout << "   TryMergeOverlapContigs: indexes:  " << overlap_record.ref_begin << "  " << overlap_record.ref_end << "  " << overlap_record.target_begin << " " << overlap_record.target_end << endl; 
  string ref_head = ref_contig.sequence.substr(0, ref_index + overlap_record.ref_begin);
  string target_head = target_contig.sequence.substr(0, overlap_record.target_begin);
  string ref_tail = ref_contig.sequence.substr(ref_index + overlap_record.ref_end + mer_len,
      ref_contig.sequence.length() - ref_index - overlap_record.ref_end - mer_len
  );
  string target_tail = target_contig.sequence.substr(overlap_record.target_end + mer_len,
      target_contig.sequence.length() - overlap_record.target_end - mer_len
  );
  //cout << "   TryMergeOverlapContigs: ref_head: " << ref_head << endl;
  //cout << "   TryMergeOverlapContigs: target_head: " << target_head << endl;
  //cout << "   TryMergeOverlapContigs: ref_tail: " << ref_tail << endl;
  //cout << "   TryMergeOverlapContigs: target_tail: " << target_tail << endl; 
  string consensus_head, consensus_tail;
  // if consensus can be constructed
  if(!FindConsensusHead(ref_head, target_head, consensus_head) || 
     !FindConsensusTail(ref_tail, target_tail, consensus_tail))  {
    return false;
  }
  //cout << "   TryMergeOverlapContigs: consensus_head:  " << consensus_head << endl;
  //cout << "   TryMergeOverlapContigs: consensus_tail:  " << consensus_tail << endl;
  // if yes, update contig
  ContigType merged_contig;
  merged_contig.sequence = consensus_head + target_contig.sequence.substr(
         overlap_record.target_begin, overlap_record.target_end - overlap_record.target_begin + mer_len
      ) + consensus_tail;
  if(handle_index)  {
    //cout << " TryMergeOverlapContigs: sequence: " << merged_contig.sequence << endl;    
    merged_contig.q_begin = ref_contig.q_begin;
    merged_contig.q_end = ref_contig.q_end;
    if(target_head.length() > ref_head.length())  {
      // update the query begin and reference index
      ref_index += target_head.length() - ref_head.length();
      if(target_contig.q_begin < merged_contig.q_begin)  {
        merged_contig.q_begin = target_contig.q_begin;
      }
    } 
    if(target_tail.length() > ref_tail.length())  {
      if(target_contig.q_end > merged_contig.q_end)  {
        merged_contig.q_end = target_contig.q_end;
      }
    } 
  }
  // update the reference contig
  //cout << " TryMergeOverlapContigs: q ends: " << merged_contig.q_begin << " " << merged_contig.q_end << endl;
  ref_contig = merged_contig;
  return true;
}

bool ContigRefinement::FindConsensusHead(
    std::string& seq1, std::string& seq2, std::string& consensus
)  {
  //cout << "   FindConsensusHead:  begin" << endl;
  if(seq2 == "" || (seq1.length() >= seq2.length() && 
     seq1.substr(seq1.length() - seq2.length(), seq2.length()) == seq2))  {
    consensus = seq1;
    //cout << "   FindConsensusHead:  end" << endl;
    return true;
  } else if(seq1 == "" || (seq2.length() >= seq1.length() && 
     seq2.substr(seq2.length() - seq1.length(), seq1.length()) == seq1)) {
    consensus = seq2;
    //cout << "   FindConsensusHead:  end" << endl;
    return true;
  }
  //cout << "   FindConsensusHead:  end" << endl;
  return false;
}

bool ContigRefinement::FindConsensusTail(
    std::string& seq1, std::string& seq2, std::string& consensus
)  {
  //cout << "   FindConsensusTail:  begin" << endl;
  if(seq2 == "" || (seq1.length() >= seq2.length() && 
     seq1.substr(0, seq2.length()) == seq2))  {
    consensus = seq1;
    //cout << "   FindConsensusTail:  end" << endl;
    return true;
  } else if(seq1 == "" || (seq2.length() >= seq1.length() && 
     seq2.substr(0, seq1.length()) == seq1)) {
    consensus = seq2;
    //cout << "   FindConsensusTail:  end" << endl;
    return true;
  }
  //cout << "   FindConsensusTail:  end" << endl;
  return false;
  
}

void ContigRefinement::FormatForPrint(
    std::string& stem,
    std::list<ContigType>& refined_contigs, 
    std::string& print_content
) {
  AssembleExtend header_foo;
  int index = 0;
  stringstream print_aln_out;
  for(auto it = refined_contigs.begin(); it != refined_contigs.end(); ++ it) {
    /**************************/
    //cout << "positions: " << it->q_begin << " " << it->q_end << endl;
    //for(auto it_m = it->al_print.nuc_match.begin(); it_m != it->al_print.nuc_match.end(); ++ it_m) {
    //  cout << it_m->first << "  " << it_m->second << endl;
    //}
    /**************************/
    stringstream outs;
    outs << stem << "||contig_" << index;
    string barcode = outs.str();  
    string header = header_foo.EncodeHeader(barcode, *it);
    ++ index;
    // finding starting and ending locations of the alignment;
    int q_start = 999999; 
    int c_start = 999999;
    for(auto it_m = it->al_print.nuc_match.begin(); 
        it_m != it->al_print.nuc_match.end(); ++ it_m
    ) {
      if(it_m->first < q_start) q_start = it_m->first;
      if(it_m->second < c_start) c_start = it_m->second;
    }
    q_start += it->q_begin;
    // count gaps, positives, identieis etc
    vector<int> seq1_index(it->al_print.symbol.length(), 0);
    vector<int> seq2_index(it->al_print.symbol.length(), 0);
    int s1x = 0, s2x = 0, num_gap = 0, num_positive = 0, num_identical = 0;
    int n = it->al_print.symbol.length();
    for(int i = 0; i < n; ++ i) {
      if(it->al_print.seq1[i] != '-' && it->al_print.seq2[i] != '-')  {
        s1x ++; s2x ++;
        if(it->al_print.symbol[i] == it->al_print.seq1[i] && 
            it->al_print.symbol[i] == it->al_print.seq2[i])  {
          ++ num_identical; ++ num_positive;
        } else if(it->al_print.symbol[i] != ' ') {
          ++ num_positive;
        }
      } else if(it->al_print.seq1[i] != '-') {
        s1x ++; num_gap ++;
      } else if(it->al_print.seq2[i] != '-') {
        s2x ++; num_gap ++;
      }
      seq1_index[i] = s1x; seq2_index[i] = s2x; 
    }
    // format output
    double perc_gap = 100 * num_gap / n;
    double perc_positive = 100 * num_positive / n;
    double perc_identical = 100 * num_identical / n;
    print_aln_out << "> " << header << endl;
    print_aln_out << "Length=" << n << endl << endl;
    print_aln_out << " Score = " << it->bit_score << " bits " << "(" << it->score << "),";
    print_aln_out << " Expect = " << it->e_value << ",";
    print_aln_out << " Method: N/A." << endl;
    print_aln_out << " Identities = " << num_identical << "/" << n << " (" << perc_identical << "\%),";
    print_aln_out << " Positives = " << num_positive << "/" << n << " (" << perc_positive << "\%),";
    print_aln_out << " Gaps = " << num_gap << "/" << n << " (" << perc_gap << "\%),";
    print_aln_out << endl << endl;
    int max_id_num = q_start + seq1_index[n - 1] + 1;
    max_id_num = c_start + seq2_index[n - 1] + 1 > max_id_num ? c_start + seq2_index[n - 1] + 1 : max_id_num;
    string full_len_str = to_string(max_id_num);
    //cout << full_len_str << " " << q_start << " " << seq1_index[n - 1] << endl;
    for(int i = 0; i < n; i += 60) {
      int al = n - i > 60 ? 60 : n - i;
      string lqb_s = FixedWidthString(full_len_str.length(), q_start + seq1_index[i]); 
      string lqe_s = FixedWidthString(full_len_str.length(), q_start + seq1_index[i + al - 1]);
      string lcb_s = FixedWidthString(full_len_str.length(), c_start + seq2_index[i]);
      string lce_s = FixedWidthString(full_len_str.length(), c_start + seq2_index[i + al - 1]);
      string sym_holder(full_len_str.length() + 9, ' ');
      print_aln_out << "Query  " << lqb_s << "  " << it->al_print.seq1.substr(i, al) << "  " << lqe_s << endl;
      print_aln_out << sym_holder << it->al_print.symbol.substr(i, al) << endl;
      print_aln_out << "Sbjct  " << lcb_s << "  " << it->al_print.seq2.substr(i, al) << "  " << lce_s << endl;
      print_aln_out << endl;
    }
    print_aln_out << endl;
  }
  print_content = print_aln_out.str();
  return;
}

std::string ContigRefinement::FixedWidthString(int len, int num)  {
  assert(len > 0);
  string num_s = to_string(num); 
  if(len < (int) num_s.length())  {
    cout << "Error: ContigRefinement::FixedWidthString: length is too short to include all digits" << endl;
    cout << "max_num_digits:  " << len << ",  output: " << num_s << endl;
    //exit(0);
  } else  {
    num_s += string(len - num_s.length(), ' ');
  }
  return num_s;
}
