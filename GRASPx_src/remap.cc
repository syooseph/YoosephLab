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

double ReMap::CalErrorScore(int num_errors, double aligned_portion)  {
  // 1 is added as pseudo-count
  return (double) (num_errors + 1) / aligned_portion;
}

bool ReMap::CompMapQuality(bool is_qual, MapInfoType &m1, MapInfoType &m2) {
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
      //cout << "search of: " << mer << endl;
      for(IdxType j = range.first; j <= range.second; ++ j) {
        if(!gap_allowed)  {
          //cout << "******************" << endl;
          RIDType read_id = seq_obj.suffix_array_->getId(j);
          if(touched_reads.find(read_id) != touched_reads.end())  {
            continue;
          }
          string read_seq = seq_obj.sequence_[read_id];
          string c_5p, c_3p, r_5p, r_3p;
          
          //cout << it->sequence << endl << read_seq << endl;
          //cout << j << "  " << (int) seq_obj.suffix_array_->getPos(j) << endl;
          int len_5p = Get5PSeq(
              0, it->sequence, i, read_seq, 
              seq_obj.suffix_array_->getPos(j), c_5p, r_5p
          );
          int len_3p = Get3PSeq(
              0, it->sequence, i + seed_len, read_seq, 
              seq_obj.suffix_array_->getPos(j) + seed_len, c_3p, r_3p
          );
          // compute aligned portion and count num. of errors
          double aligned_portion = (double) (len_5p + len_3p + seed_len) / read_seq.length();
          double num_errors = CountUngappedError(c_5p, r_5p) + CountUngappedError(c_3p, r_3p);
          MapInfoType m_info;
          m_info.contig_id = contig_id;
          m_info.e_value = it->e_value;
          m_info.error_score = CalErrorScore(num_errors, aligned_portion);
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
    mapped_reads.push_back(m);
  }
  return;
}
