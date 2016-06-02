#include "kmer_distance.h"
#include "sequence.h"
#include "gsa.h"
#include "time.h"
#include "seq_align.h"
#include <unordered_map>
#include <iostream>
#include <cstring>
#include <string>

typedef int AlignScoreType;
static int kmer_size = 10;
static double conv_rate = 0.8;

using namespace std;

int main(int argc, char **argv)  {
  // read the query sequence
  
  string qfile = argv[1];
  vector<string> qfiles_in;
  qfiles_in.push_back(qfile);
  seq::totalSequenceCount(qfiles_in);
  char **query_seqs_ = new char* [1];
  char **q_tag_foo = NULL;
  seq::loadSequences(qfiles_in, q_tag_foo, query_seqs_, SEQONLY);
  string query = query_seqs_[0];
  // load the target sequences
  vector<string> files_in;
  files_in.resize(1);
  files_in[0] = argv[2];
  int sample_num_reads = (unsigned int) seq::totalSequenceCount(files_in);
  char **sample_seqs = new char *[sample_num_reads];
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, sample_seqs, SEQONLY);
  delete tag_foo;
  // build kmer closure implicitly
  double start_time = mytime();
  KmerDistance index_closure(kmer_size, ALL20);
  index_closure.KmerAdjacency(files_in[0]); 
  double edge_end_time = mytime();
  //printElapsed(start_time, edge_end_time, "Constructing kmer-relation graph");
  // define the kmers as seeds
  start_time = mytime();
  string index_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.idx";
  string position_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.pos";
  unordered_map<KmerType, list<PositionType> > kmer_positions;
  IndexSample index_loader;
  index_loader.ReadKmerIndex(
      index_file_name, position_file_name,
      query, kmer_positions
  );
  unordered_map<string, bool> seed_kmers;
  IndexSample encoder(kmer_size, ALL20);
  ScoringFunction<AlignScoreType> *score_scheme_ = new ScoringFunction<AlignScoreType>(PROTEIN, BLOSUM62, -1,
 -11);
  for(int i = 0; i < query.length() - 6 + 1; ++ i)  {
      KmerType encoded_kmer = index_loader.encode_kmer(query.substr(i, 6));
      for(auto it_j = kmer_positions[encoded_kmer].begin(); 
          it_j != kmer_positions[encoded_kmer].end(); ++ it_j
      ) {
        string read_seq = sample_seqs[it_j->rid];
        string matched_kmer = read_seq.substr(it_j->pos, 6);
        string original_kmer = query.substr(i, 6);
        // check alignment score
        AlignScoreType est_score = 0;
        for(int i2 = 0; i2 < original_kmer.length(); ++ i2) {
          est_score += score_scheme_->CheckMatchScore(original_kmer[i2], matched_kmer[i2]);
        }
        if(est_score >= (AlignScoreType) (0.6 * score_scheme_->GetAveMatch() * 6))  {
          seed_kmers[read_seq.substr(0, kmer_size)] = true;
          seed_kmers[read_seq.substr(read_seq.length() - kmer_size, kmer_size)] = true;
        }
      }
  }
  unordered_map<string, bool> closure_kmers;
  index_closure.GetKmerClosure(0.8, seed_kmers, closure_kmers);
  double build_end_time = mytime();
  //printElapsed(start_time, build_end_time, "Resolve kmers in the closure");
  start_time = mytime();
  // define the sequences based on the kmers presented
  unordered_map<RIDType, bool> closure_reads;
  string lcp_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.lcp";
  string mcp_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.mcp";
  string gsa_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.faa.gsa";
  GSA *suffix_array;
  suffix_array = new GSA();
  suffix_array->load(
    lcp_file_name.c_str(), mcp_file_name.c_str(), gsa_file_name.c_str()
  );
  suffix_array->setSequences(sample_seqs);
  suffix_array->setReadCount(sample_num_reads);
  for(auto it = closure_kmers.begin(); it != closure_kmers.end(); ++ it) {
    BoundType seed_range = suffix_array->searchWithLCPs((SfaChar*) it->first.c_str(), it->first.length());
    for(int j = seed_range.first; j <= seed_range.second; ++ j) {
      int read_ID = (int) suffix_array->getId(j);
      closure_reads[read_ID] = 1;
    }
  }
  double fetch_end_time = mytime();
  //printElapsed(start_time, fetch_end_time, "Loading sequence based on kmers in closure");
  cerr << "Num of final recruited reads:  " << closure_reads.size() << "/" << sample_num_reads << endl;
  start_time = mytime();
  // load reads and build a new suffix array on it
  int selected_num_reads = closure_reads.size();
  char **selected_seqs = new char*[selected_num_reads];
  int index = 0;
  for(auto it = closure_reads.begin(); it != closure_reads.end(); ++ it) {
    selected_seqs[index] = new char[strlen(sample_seqs[it->first]) + 1];
    strcpy(selected_seqs[index], sample_seqs[it->first]);
    ++ index;
  }
  GSA* reduced_SA = new GSA(selected_seqs, selected_num_reads, true);
  char **reversed_selected_seqs = new char *[selected_num_reads];
  seq::reverseSequences(selected_num_reads, selected_seqs, reversed_selected_seqs);
  GSA* reverse_reduced_SA = new GSA(reversed_selected_seqs, selected_num_reads, true);
  // print the sequences
  for(int i = 0; i < selected_num_reads; ++ i) {
    cout << ">read_" << i << endl << selected_seqs[i] << endl; 
  }
  // build a suffix array on it
  
  double sa_end_time = mytime();
  printElapsed(start_time, sa_end_time, "Build suffix array on selected reads");
  return 0;
}
