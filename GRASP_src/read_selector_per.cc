#include "index_sample.h"
#include "gsa.h"
#include "sfa.h"
#include "interval_array.h"
#include "sequence.h"
#include "seq_align.h"
#include "scoring_function.h"

using namespace std;

typedef int AlignScoreType;
static int seed_len = 6;
static int n_back = 10;

int main(int argc, char** argv)  {
  // load the sequence
  vector<string> files_in;
  files_in.resize(1);
  files_in[0] = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/saliva.faa";
  int sample_num_reads_ = (unsigned int) seq::totalSequenceCount(files_in);
  char **sample_seqs_ = new char *[sample_num_reads_];
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, sample_seqs_, SEQONLY);
  delete tag_foo;
  // load the indexing file
  string index_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.idx";
  string position_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.pos";
  string lcp_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.lcp";
  string mcp_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.mcp";
  string gsa_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.gsa";
  GSA *suffix_array_;
  suffix_array_ = new GSA();
  suffix_array_->load(
    lcp_file_name.c_str(), mcp_file_name.c_str(), gsa_file_name.c_str()
  );
  suffix_array_->setSequences(sample_seqs_);
  suffix_array_->setReadCount(sample_num_reads_);
  //cerr << "All sequences loaded" << endl;
  for(int i = 0; i < sample_num_reads_; ++ i)	{
    unordered_map<KmerType, int> searched_kmers;
    unordered_map<int, int> final_reads;
    unordered_map<int, int> candidate_reads;
    // putting the current read as a candidate
    candidate_reads[i] = 1;
    int num_final_reads = 0;
    IndexSample encoder(n_back, ALL20);
    while(1) {
      int num_final_reads = final_reads.size();
      // check each candidate read
      unordered_map<int, int> reads_next_stage;
      int reads_processed = 0;
      int reads_in_next = 0;
      for(auto it_r = candidate_reads.begin(); it_r != candidate_reads.end(); ++ it_r) {
        // if the reads has been parsed, skip
        ++ reads_processed;
        //cerr << "read processed: " << reads_processed << "\t" << reads_next_stage.size()  << endl;
        if(final_reads.find(it_r->first) != final_reads.end())  {
          continue;
        }
        final_reads[it_r->first] = 1; 
        string pivot_seq = sample_seqs_[it_r->first];
        for(int j = 0; j < pivot_seq.length() - n_back + 1; ++ j)  {
          if(j != 0 && j != pivot_seq.length() - n_back)  {
            continue;
          }
          string sample_seq = pivot_seq.substr(j, n_back);
          KmerType encoded_kmer = encoder.encode_kmer(sample_seq);
          if(searched_kmers.find(encoded_kmer) != searched_kmers.end())  {
            continue;
          }
          searched_kmers[encoded_kmer] = 1;
          BoundType seed_range = suffix_array_->searchWithLCPs((SfaChar*) sample_seq.c_str(), sample_seq.length());
          for(int j = seed_range.first; j <= seed_range.second; ++ j) {
            int read_ID = (int) suffix_array_->getId(j);
            if(final_reads.find(read_ID) == final_reads.end())  {
              reads_next_stage[read_ID] = 1;
            }
          }
        }
      }
      if(((double) num_final_reads / (double) final_reads.size()) >= 0.95)  {
        break;
      }
      candidate_reads = reads_next_stage;
      num_final_reads = final_reads.size();
      //cerr << num_final_reads << endl;
    
    }
    cout << "read_recruitment:\t" << i << "\t" << final_reads.size() << endl;
    //for(auto it_fr = final_reads.begin(); it_fr != final_reads.end(); ++ it_fr)  {
    //  cout << ">read_" << it_fr->first << endl;
    //  cout << sample_seqs_[it_fr->first] << endl; 
    //}
  }
  return 0;
}
