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
static int n_back = 6;

int main(int argc, char** argv)  {
  string input = argv[1];
  //cout << input << endl;
  //getline(cin, input);
  int num_iter = stoi(input);
  string query = "MKSYQKTIFGEFQGQNIFRFTFENDLGYRLSVMNYGATILEYQTPDKKGQFANVILGFDQFEDYIGNSPKHGASIGPVAGRIAGASFELGGEIYHLEANNGQNCNHSGSTGWDSTVFQVEEVTDEGLVLFTERADGTGGFPGHLKVWISYTLSEKGELEISYQVQTDRDTLINPTNHSYFNLSADFAQSIDDHVFQVDSLGFYPIAEDGVPAKETEISDFVKHLQQAMLLKDLFAEKDEQVRLVSGLDHPFALKPGHETAGFLYHQESGRFLTFKTEAPCLVVYTANCVDEAIRFGGQTMKQHNGVALEMQALPDAIHSHQKDQVIVRAGQEFTSTTTYHAIAK";
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
  unordered_map<KmerType, list<PositionType> > kmer_positions;
  IndexSample index_loader;
  index_loader.ReadKmerIndex(
      index_file_name, position_file_name, 
      query, kmer_positions
  );
  // load the suffix array
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
  cerr << "Suffix array loaded" << endl;
  // check how many reads are included
  unordered_map<int, int> candidate_reads;
  ScoringFunction<AlignScoreType> *score_scheme_ = new ScoringFunction<AlignScoreType>(PROTEIN, BLOSUM62, -1, -11);
  int i;
  for(i = 0; i < query.length() - seed_len + 1; ++ i) {
    string current_kmer = query.substr(i, seed_len);
    KmerType encoded_kmer = index_loader.encode_kmer(current_kmer);
    list<PositionType> seed_candidates = kmer_positions[encoded_kmer];
    unordered_map<string, int> translated_kmers;
    for(auto it = seed_candidates.begin(); it != seed_candidates.end(); ++ it) {
      string seed_seq = sample_seqs_[(int) it->rid];
      string sample_kmer = seed_seq.substr(it->pos, seed_len);
      translated_kmers[sample_kmer] = (int) it->rid;
    }
    for(auto it = translated_kmers.begin(); it != translated_kmers.end(); ++ it) {
      // gets the exact k-mer sequence in the sample
      //string seed_seq = suffix_array_->getSuffix_explicit(it->rid, it->pos);
      
      AlignScoreType est_score = 0;
      for(int i2 = 0; i2 < current_kmer.length(); ++ i2) {
        est_score += score_scheme_->CheckMatchScore(current_kmer[i2], (it->first)[i2]);
      }
      // if the alignment score is high enough
      if(est_score >= (AlignScoreType) (0.6 * score_scheme_->GetAveMatch() * seed_len))  {
        candidate_reads[it->second] = 1;
      }
    }
  }
  cerr << "All sequences loaded" << endl;
  int num_rounds = (int) (query.length() / 30);
  num_rounds = num_iter;
  unordered_map<KmerType, int> searched_kmers;
  unordered_map<int, int> final_reads;
  int num_final_reads = 0;
  //IndexSample encoder(n_back, ALL20);
  while(num_rounds > 0) {
    //cerr << "new round begins with:\t" << candidate_reads.size() << endl;
    -- num_rounds;
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
        KmerType encoded_kmer = index_loader.encode_kmer(sample_seq);
        if(searched_kmers.find(encoded_kmer) != searched_kmers.end())  {
          continue;
        }
        searched_kmers[encoded_kmer] = 1;
        list<PositionType> candidate_reads = kmer_positions[encoded_kmer];
        for(auto it = candidate_reads.begin(); it != candidate_reads.end(); ++ it) {
          int read_ID = it->rid;
          if(final_reads.find(read_ID) == final_reads.end())  {
            reads_next_stage[read_ID] = 1;
          }
        }
      }
    }
    //if(num_final_reads >= final_reads.size())  {
    //  break;
    //}
    candidate_reads = reads_next_stage;
    num_final_reads = final_reads.size();
    cerr << num_final_reads << endl;
    
  }
  //for(auto it_fr = final_reads.begin(); it_fr != final_reads.end(); ++ it_fr) {
  //  cout << it_fr->first << endl;
  //}
  //for(auto it_fr = final_reads.begin(); it_fr != final_reads.end(); ++ it_fr)  {
  //  cout << ">read_" << it_fr->first << endl;
  //  cout << sample_seqs_[it_fr->first] << endl; 
  //}
  cerr << "Total kmer searched:	" << searched_kmers.size() << endl;
  return 0;
}
