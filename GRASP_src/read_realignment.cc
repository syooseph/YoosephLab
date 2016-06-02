#include "read_realignment.h"

using namespace std;

ReadRealignment::ReadRealignment()  {
  return;
}

ReadRealignment::~ReadRealignment() {
  return;
}

ReadRealignment::ReadRealignment(
    const std::string& in_query_seq,
    const int& in_query_begin, const int& in_query_end,
    const std::string& in_assembled_seq,
    ScoringFunction<AlignScoreType>& in_score_scheme,
    const unsigned int& in_down_band, const unsigned int& in_right_band,
    std::list<ReadAssignmentType>& in_path_reads
) {
  query_seq_ = in_query_seq;
  query_begin_ = in_query_begin;
  query_end_ = in_query_end;
  assembled_seq_ = in_assembled_seq;
  score_scheme_ = &in_score_scheme;
  right_band_ = in_right_band;
  down_band_ = in_down_band;
  path_reads_ = in_path_reads;
  return;
}

void ReadRealignment::Align() {
  unordered_map<int, int> nuc_match;
  string aligned_query, aligned_assembled, aligned_symbol;
  AlignScoreType opt_score;
  double opt_bit_score;
  ReAlignSequences(opt_score, opt_bit_score, nuc_match, aligned_query, aligned_assembled, aligned_symbol);
  //FormatOutput(opt_score, opt_bit_score, nuc_match, aligned_query, aligned_assembled, aligned_symbol);
  return;
}

bool _cmp_read_assignment(const ReadAssignmentType& first, const ReadAssignmentType& second) {
  if(first.offset_to_start < second.offset_to_start)  {
    return true;
  } else  {
    return false;
  }
}

/*
void ReadRealignment::FormatOutput(
    AlignScoreType& opt_score, double& opt_bit_score,
    std::unordered_map<int, int>& nuc_match, 
    std::string& aligned_query, std::string& aligned_assembled, std::string& aligned_symbols
)  {
  // compute the aligned index
  int match_begin = nuc_match[0];
  int match_end = nuc_match[assembled_seq_.length() - 1];
  int aligned_query_begin = query_begin_;
  if(match_begin > 0)  {
    aligned_query_begin += match_begin;
  }
  int aligned_query_end = query_end_;
  if(match_end > 0 && match_end < query_end_ - query_begin_)  {
    int diff = query_end_ - query_begin_ - match_end;
    aligned_query_end -= diff;
  }
  // build a map of the aligned nucleotides in the actual alignment
  map<int, int> nuc_in_alignment; // the first int is the index in the assembled sequence, and the second int is the index in the alignment
  int i, seq_index = -1;
  for(i = 0; i < (int) aligned_assembled.length(); ++ i) {
    if(aligned_assembled[i] != '-')  {
     ++ seq_index;
     nuc_in_alignment[seq_index] = i;
    }
  }
  // TODO: construct the alignment of the individual reads
  path_reads_.sort(_cmp_read_assignment);
  list<string> reads_alignment;
  for(auto it = path_reads_.begin(); it != path_reads_.end(); ++ it) {
    int aln_begin = nuc_in_alignment[it->offset_to_start];
    int aln_end = nuc_in_alignment[it->offset_to_start + it->full_seq_len - 1];
    string fill_alignment = aligned_assembled.substr(aln_begin, aln_end - aln_begin + 1);
    string flanking_prefix = string(aln_begin, '.');
    string flanking_suffix = string(aligned_assembled.length() - aln_end - 1, '.');
    reads_alignment.push_back(flanking_prefix + fill_alignment + flanking_suffix);
  }
  
  cout << "Alignment score: " << opt_bit_score << " bits (" << opt_score << ")" << endl;
  cout << "Aligned regions: " << aligned_query_begin << " " << aligned_query_end << endl;
  for(i = 0; i < (int) aligned_query.length(); i += 100)	{
		cout << "\t" << aligned_query.substr(i, 100) << endl;
		cout << "\t" << aligned_symbols.substr(i, 100) << endl;
		cout << "\t" << aligned_assembled.substr(i, 100) << endl;
		for(auto it = reads_alignment.begin(); it != reads_alignment.end(); ++ it) {
      cout << "\t" << it->substr(i, 100) << endl;
    }
		std::cout << std::endl;
	}
  return;
}
*/

void ReadRealignment::ReAlignSequences(
    AlignScoreType& opt_score, double& opt_bit_score,
    std::unordered_map<int, int>& nuc_match, 
    std::string& aligned_query, std::string& aligned_assembled, std::string& aligned_symbols
)  {
  string partial_query_seq = query_seq_.substr(query_begin_, query_end_ - query_begin_ + 1);
  int seq_len_diff = assembled_seq_.length() - partial_query_seq.length();
  if(seq_len_diff < 0)  {
    seq_len_diff = 0;
  }
  SeqAlign<AlignScoreType> realign(
      assembled_seq_, partial_query_seq, score_scheme_, 
      SEMIGLOBAL, down_band_ + seq_len_diff, down_band_ + right_band_ + seq_len_diff, false
  );

  realign.Align();
  realign.TraceBack();
  realign.GetAlignment(nuc_match, aligned_assembled, aligned_query, aligned_symbols);
  opt_score = realign.GetBestGlobalScore();
  opt_bit_score = score_scheme_->ComputeBitScore(opt_score);
  return;
}
