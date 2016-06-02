#include "index_sample.h"
#include "gsa.h"
#include "sfa.h"
#include "interval_array.h"
#include "sequence.h"
#include "seq_align.h"
#include "scoring_function.h"
#include <string>

using namespace std;

typedef int AlignScoreType;
static int seed_len = 6;
static int n_back = 20;

int main()  {
  vector<string> files_in(1);
  files_in[0] = "../Works/GuidedAssemble/Data/marine.core.faa";
  int sample_num_reads_ = (unsigned int) seq::totalSequenceCount(files_in);
  char **sample_seqs_ = new char *[sample_num_reads_];
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, sample_seqs_, SEQONLY);
  delete tag_foo;
  //cout << "here" << endl;
  // load the suffix array
  string lcp_file_name = "./WorkSpace/marine.core.faa.lcp";
  string mcp_file_name = "./WorkSpace/marine.core.faa.mcp";
  string gsa_file_name = "./WorkSpace/marine.core.faa.gsa";
  GSA *suffix_array_;
  suffix_array_ = new GSA();
  suffix_array_->load(
    lcp_file_name.c_str(), mcp_file_name.c_str(), gsa_file_name.c_str()
  );
  suffix_array_->setSequences(sample_seqs_);
  suffix_array_->setReadCount(sample_num_reads_);
  unordered_map<int, bool> unique_reads;
  for(int i = 1; i < suffix_array_->getSize() - 1; ++ i)  {
    //cout << (int) suffix_array_->getLcp(static_cast<size_t>(i)) << endl;
    if(suffix_array_->getSuffixLength(i) >= n_back && 
        (unsigned int) suffix_array_->getLcp(i) < n_back && 
        (unsigned int) suffix_array_->getLcp(i + 1) < n_back)  {
      //cout << i << "	" << (unsigned int) suffix_array_->getLcp(i) << "	" << (unsigned int) suffix_array_->getLcp(i + 1) << endl;
      unique_reads[suffix_array_->getId(i)] = true;
    }
  }
  for(auto it = unique_reads.begin(); it != unique_reads.end(); ++ it)  {
    cout << it->first << endl;
  }
  return 0;
  IndexSample encoder(n_back, ALL20);
  unordered_map<KmerType, unordered_map<KmerType, bool> > adj_list;
  unordered_map<KmerType, bool> searched_kmers;
  //cout << "here" << endl;
  for(int i = 0; i < sample_num_reads_; ++ i) {
    string pivot_seq = sample_seqs_[i];
    if(i % 100000 == 0) cerr << i << " sequence processed" << endl; 
    for(int j = 0; j < pivot_seq.length() - n_back + 1; ++ j)  {
      //if(j != 0 && j != pivot_seq.length() - 6)  {
      //  continue;
      //}
      string sample_seq = pivot_seq.substr(j, n_back);
      KmerType encoded_kmer = encoder.encode_kmer(sample_seq);
      if(searched_kmers.find(encoded_kmer) != searched_kmers.end())  {
        continue;
      }
      searched_kmers[encoded_kmer] = true;
      BoundType seed_range = suffix_array_->searchWithLCPs((SfaChar*) sample_seq.c_str(), sample_seq.length());
      for(int j = seed_range.first; j <= seed_range.second; ++ j) {
        int read_ID = (int) suffix_array_->getId(j);
        string full_seq = suffix_array_->getSequence_explicit(read_ID);
        KmerType encoded_fw_nmer = encoder.encode_kmer(full_seq.substr(0, n_back));
        KmerType encoded_bk_nmer = encoder.encode_kmer(full_seq.substr(full_seq.length() - n_back, n_back));
        adj_list[encoded_kmer][encoded_fw_nmer] = true;
        adj_list[encoded_fw_nmer][encoded_kmer] = true;
        adj_list[encoded_kmer][encoded_bk_nmer] = true;
        adj_list[encoded_bk_nmer][encoded_kmer] = true;
      }
    }
    //if(num_final_reads >= final_reads.size())  {
    //  break;
    //}
    //candidate_reads = reads_next_stage;
    //num_final_reads = final_reads.size();
    //cerr << num_final_reads << endl; 
  }
  //for(auto it_fr = final_reads.begin(); it_fr != final_reads.end(); ++ it_fr) {
  //  cout << it_fr->first << endl;
  //}
  //for(auto it_fr = final_reads.begin(); it_fr != final_reads.end(); ++ it_fr)  {
  //  cout << ">read_" << it_fr->first << endl;
  //  cout << sample_seqs_[it_fr->first] << endl; 
  //}
  //cout << "graph	{" << endl;
  //for(auto it = adj_list.begin(); it != adj_list.end(); ++ it)  {
  //  for(auto it_j = it->second.begin(); it_j != it->second.end(); ++ it_j)  {
  //    cout << "	" << it->first << " -> " << it_j->first << ";" << endl; 
  //  }
  //}
  //cout << "}" << endl;
  for(auto it = adj_list.begin(); it != adj_list.end(); ++ it)  {
    cout << it->second.size() << endl;
  }

  return 0;
}
