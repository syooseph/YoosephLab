#include "read_alignment.h"

using namespace std;

ReadAlignment::ReadAlignment(void)  {
  return;
}

ReadAlignment::~ReadAlignment(void) {
  return;
}

std::string ReadAlignment::GetQueryFragment(
    std::string& query, int begin, 
    int end, int band_size
) {
  // expend to band size
  begin = begin - (int) (band_size / 2);
  end = end + (int) (band_size / 2);
  // put ensure valid ranges
  begin = begin < 0 ? 0 : begin;
  end = end >= (int) query.length() ? (int) query.length() - 1 : end;
  if(begin >= (int) query.length() || end < 0 || begin > end)  {
    return "";
  }
  // return the sequence
  return query.substr(begin, end - begin + 1);
}

void ReadAlignment::ComputeAlignmentScore(
    std::string& query, SequenceBuild& seq_obj,
    ScoringFunction<int>& score_obj, 
    int band_size, std::list<ReadType>& candidate_reads
) {
  // compute alignment for each read regarding their local regions
  for(auto it = candidate_reads.begin(); it != candidate_reads.end(); ++ it) {
    string seq_read = seq_obj.sequence_[it->rid];
    string seq_query = GetQueryFragment(query, it->q_begin, it->q_end, band_size);
    //cout << ">read_" << it->rid << endl << seq_read << endl;
    if(!seq_read.empty() && !seq_query.empty())  {
      SeqAlign<int> align(seq_read, seq_query, &score_obj, LOCAL);
      align.Align();
      int aln_score = align.GetBestScore();
      it->aln_e_value = score_obj.ComputeEValue(
          seq_read.length(), query.length(), aln_score
      );
    } else  {
      it->aln_e_value = _MAX_VALUE;
    }
  }
  return;
}

bool cmp_read_evalue(const ReadType& a, const ReadType& b)  {
  if(a.aln_e_value < b.aln_e_value)  {
    return true;
  }
  return false;
}

void ReadAlignment::SortReadsOnEvalue(std::list<ReadType>& candidate_reads)  {
  candidate_reads.sort(cmp_read_evalue);
  return;
}
