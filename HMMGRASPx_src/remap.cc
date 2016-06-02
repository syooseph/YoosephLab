#include "remap.h"

using namespace std;

ReMap::ReMap() {
  contigs_loaded_ = false;
  num_seqs_ = 0;
  return;
}

ReMap::~ReMap() {
  if(contigs_loaded_)  {
    delete [] header_;
    delete [] sequence_;
  }
  return;
}

std::string ReMap::GetContigInfo(int contig_id)  {
  if(contig_id < num_seqs_)  {
    return string(header_[contig_id]);
  }
  return "Unknown_contig";
}

void ReMap::LoadContigs(
    std::string& contig_file, std::list<ContigType>& loaded_contigs
) {
  vector<string> files_in;
  files_in.push_back(contig_file);
  num_seqs_ = (unsigned int) seq::totalSequenceCount(files_in);
  header_ = new char* [num_seqs_];
  sequence_ = new char* [num_seqs_];
  seq::loadSequences(files_in, header_, sequence_, TAGSEQ);
  // construct the ContigTypes
  for(int i = 0; i < num_seqs_; ++ i) {
    ContigType c;
    c.sequence = sequence_[i];
    loaded_contigs.push_back(c);
  }
  contigs_loaded_ = true;
  return;
}

int ReMap::CountUngappedError(std::string& s1, std::string& s2) {
  int ne = 0;
  if(s1.length() != s1.length())  {
    cout << "Error: ReMap::CountUnGappedError:  ungapped matching applied on sequences with different lengths" << endl;
    return 9999;
  }
  for(unsigned int i = 0; i < s1.length(); ++ i) {
    if(s1[i] != s2[i])  {
      ++ ne;
    }
  }
  return ne;
}

int ReMap::Get5PSeq(
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

int ReMap::Get3PSeq(
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



void ReMap::MapToContig(
    SequenceBuild& seq_obj, std::list<ContigType>& contigs, 
    int seed_len, double portion, int num_errors, bool gap_allowed, 
    std::list<MapReadType>& mapped_reads
) {
  //cout << "begin of MapToContig" << endl;
  unordered_map<RIDType, MapInfoType> touched_reads;  // read ID and contig ID
  int contig_id = 0;
  for(auto it = contigs.begin(); it != contigs.end(); ++ it) {
    if((int) it->sequence.length() < seed_len)  {
      continue;
    }
    for(unsigned int i = 0; i < it->sequence.length() - seed_len; ++ i) {
      string mer = it->sequence.substr(i, seed_len); 
      pair<IdxType, IdxType> range = seq_obj.suffix_array_->searchWithLCPs(
          (SfaChar*) mer.c_str(), seed_len
      );
      
      for(IdxType j = range.first; j <= range.second; ++ j) {
        if(!gap_allowed)  {
          
          RIDType read_id = seq_obj.suffix_array_->getId(j);
          //cout << "******************" << endl;
          //cout << "search of: " << mer << " in " << read_id << endl;
          
          string read_seq = seq_obj.sequence_[read_id];
          
          /***** obsolete version
          string c_5p, c_3p, r_5p, r_3p;
          
          cout << it->sequence << endl << read_seq << endl;
          cout << j << "  " << (int) seq_obj.suffix_array_->getPos(j) << endl;
          int len_5p = Get5PSeq(
              0, it->sequence, i, read_seq, 
              seq_obj.suffix_array_->getPos(j), c_5p, r_5p
          );
          int len_3p = Get3PSeq(
              0, it->sequence, i + seed_len, read_seq, 
              seq_obj.suffix_array_->getPos(j) + seed_len, c_3p, r_3p
          );
          cout << "extensions:  " << len_5p << "  " << len_3p << endl;
          cout << "5p seqs: " << c_5p << "  " << r_5p << endl;
          cout << "3p seqs: " << c_3p << "  " << r_3p << endl;
          */
          int d_errors, m_begin, m_end;
          MapUngapped(
              it->sequence, read_seq,
              i, seq_obj.suffix_array_->getPos(j),
              seed_len, num_errors,     
              d_errors, m_begin, m_end
          );
          //cout << "New mapp results:  " << d_errors << "  " << m_begin << " " << m_end << endl;
          // compute aligned portion and count num. of errors
          double aligned_portion = (double) (m_end - m_begin + 1) / read_seq.length();
          if(aligned_portion < portion)  {
            // mapped protion does not meet criterion, skip for next read
            continue;
          }
          MapInfoType m_info;
          m_info.contig_id = contig_id;
          m_info.e_value = it->e_value;
          m_info.num_errors = d_errors;
          m_info.aligned_portion = aligned_portion;
          if(touched_reads.find(read_id) == touched_reads.end() || 
             CompMapQuality(true, touched_reads[read_id], m_info))  {
            touched_reads[read_id] = m_info;
          }
        } else  {
          cout << "Error! Remap::MapToContig: gapped alignment mode is currently not supported" << endl;
          exit(0);
        }
      }
    }    
    ++ contig_id;
  }
  for(auto it = touched_reads.begin(); it != touched_reads.end(); ++ it) {
    MapReadType m;
    m.rid = it->first;
    m.contig_id = it->second.contig_id;
    m.num_errors = it->second.num_errors;
    m.aligned_portion = it->second.aligned_portion;
    mapped_reads.push_back(m);
  }
  return;
}

void ReMap::MapUngapped(
    std::string &r_seq, std::string &t_seq,   // the input sequences
    int r_anchor, int t_anchor,               // the beginning of the anchor
    int seed_len,                             // the length of the seed match
    int num_sub_errors,                       // the number of allowed substitution errors
    int &detected_errors,                     // Outputs: number of detected errors
    int &mapped_begin, int &mapped_end        // Outputs: the indexes in r_seq
)  {
  //cout << "anchors: " << r_anchor << "  " << t_anchor << "  " << seed_len << endl;
  // the left and right extension length with up to i errors
  int *llen = new int [num_sub_errors + 1];
  int *rlen = new int [num_sub_errors + 1];
  llen[0] = rlen[0] = 0;
  // fill the error matrix
  int le = 0;
  int ir = r_anchor - 1, it = t_anchor - 1;
  while(ir >= 0 && it >= 0)  {
    //cout << "check char:  " << r_seq[ir] << " " << t_seq[it] << endl;
    if(r_seq[ir] == t_seq[it])  {
      //cout << "char match!" << endl;
      ++ llen[le];  // increase the length with le errors
    } else  {
      //cout << "char mis-match!" << endl;
      ++ le;        // increase the number of errors
      //cout << "le value:  " << le << endl;
      // record current length, should add the current mismatch
      if(le <= num_sub_errors) {
        //cout << "error extended:  " << endl;
        llen[le] = llen[le - 1] + 1;
      } else  {
        break;  
      }
    }
    -- ir; -- it;
  }
  //for(int k = 0; k <= num_sub_errors; ++ k) {
  //  cout << "!  " << llen[k] << endl;
  //}
  int re = 0;
  ir = r_anchor + seed_len; it = t_anchor + seed_len;
  while(ir < (int) r_seq.length() && it < (int) t_seq.length()) {
    //cout << "check char:  " << r_seq[ir] << " " << t_seq[it] << endl;
    if(r_seq[ir] == t_seq[it])  {
      //cout << "char match!" << endl;
      ++ rlen[re];
    } else  {
      //cout << "char mis-match!" << endl;
      ++ re;
      if(re <= num_sub_errors)  {
        rlen[re] = rlen[re - 1] + 1;
      } else  {
        break; 
      }
    }
    ++ ir; ++ it;
  }
  //for(int k = 0; k <= num_sub_errors; ++ k) {
  //  cout << "!  " << rlen[k] << endl;
  //}
  // find the longest mapping
  int i, j;
  int max_len = -1;   // the maximal extension length
  for(i = (le <= num_sub_errors ? le : num_sub_errors) ; i >= 0; -- i) {
    j = num_sub_errors - i;
    j = j < re ? j : re;
    
    int len = llen[i] + rlen[j];
    
    if(len > max_len)  {
      //cout << "index: " << i << " " << j << endl;
      //cout << "ext len: " << len << endl;
      max_len = len;
      mapped_begin = r_anchor - llen[i];
      mapped_end = r_anchor + seed_len + rlen[j] - 1;
      detected_errors = i + j;
    }
  }
  // collect memory
  delete [] llen;
  delete [] rlen;
  return;
}
