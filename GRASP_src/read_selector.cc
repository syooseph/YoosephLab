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
  string qfile = argv[1];
  vector<string> qfiles_in;
  qfiles_in.push_back(qfile);
  seq::totalSequenceCount(qfiles_in);
  char **query_seqs_ = new char* [1];
  char **q_tag_foo = NULL;
  seq::loadSequences(qfiles_in, q_tag_foo, query_seqs_, SEQONLY);
  string query = query_seqs_[0];
  //cout << query << endl;
  string input = argv[2];
  //cout << input << endl;
  //getline(cin, input);
  int num_iter = stoi(input);
  // load the sequence
  vector<string> files_in;
  files_in.resize(1);
  files_in[0] = "./saliva.1M.faa";
  int sample_num_reads_ = (unsigned int) seq::totalSequenceCount(files_in);
  char **sample_seqs_ = new char *[sample_num_reads_];
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, sample_seqs_, SEQONLY);
  delete tag_foo;
  // load the indexing file
  string index_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.1M.faa.idx";
  string position_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.1M.faa.pos";
  unordered_map<KmerType, list<PositionType> > kmer_positions;
  IndexSample index_loader;
  index_loader.ReadKmerIndex(
      index_file_name, position_file_name, 
      query, kmer_positions
  );
  // load the suffix array
  string lcp_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.1M.faa.lcp";
  string mcp_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.1M.faa.mcp";
  string gsa_file_name = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/saliva.1M.faa.gsa";
  GSA *suffix_array_;
  suffix_array_ = new GSA();
  suffix_array_->load(
    lcp_file_name.c_str(), mcp_file_name.c_str(), gsa_file_name.c_str()
  );
  suffix_array_->setSequences(sample_seqs_);
  suffix_array_->setReadCount(sample_num_reads_);
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
  //cerr << "All sequences loaded" << endl;
  //int num_rounds = (int) (query.length() / 30);
  //num_rounds = num_iter;
  //unordered_map<KmerType, int> searched_kmers;
  //unordered_map<int, int> final_reads;
  int num_final_reads = 0;
  IndexSample encoder(n_back, ALL20);
  unordered_map<int, int> all_reads_touched;
  //for(auto it_r = candidate_reads.begin(); it_r != candidate_reads.end(); ++ it_r) {
  int num_total_edges = 0;
  for(int rid = 0; rid < sample_num_reads_; ++ rid)  {
    //cerr << "new round begins with:\t" << candidate_reads.size() << endl;
    //-- num_rounds;
    // check each candidate read
    int num_rounds = num_iter;
    unordered_map<KmerType, int> searched_kmers;
    unordered_map<int, int> final_reads;
    unordered_map<int, int> reads_next_stage;
    reads_next_stage[rid] = 1;
    int reads_processed = 0;
    int reads_in_next = 0;
    int prev_read_set_size = final_reads.size();
    while(num_rounds > 0) {
      // if the reads has been parsed, skip
      //-- num_rounds;
      unordered_map<int, int> reads_current_stage;
      for(auto it_c = reads_next_stage.begin(); it_c != reads_next_stage.end(); ++ it_c)	{
        all_reads_touched[it_c->first] = 1;
        //cerr << "read processed: " << reads_processed << "\t" << reads_next_stage.size()  << endl;
        if(final_reads.find(it_c->first) != final_reads.end())  {
          continue;
        }
        final_reads[it_c->first] = 1; 
        string pivot_seq = sample_seqs_[it_c->first];
        for(int j = 0; j < pivot_seq.length() - n_back + 1; ++ j)  {
          //if(j != 0 && j != pivot_seq.length() - n_back)  {
          //  continue;
          //}
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
              reads_current_stage[read_ID] = 1;
            }
          }
        }
      }
      reads_next_stage = reads_current_stage;
      if(final_reads.size() <= prev_read_set_size)  {
        break;
      }
      prev_read_set_size = final_reads.size(); 
      //cout << "  " << final_reads.size() << endl;
    }
    //if(num_final_reads >= final_reads.size())  {
    //  break;
    //}
    //candidate_reads = reads_next_stage;
    //double ratio = (double) num_final_reads / (double) final_reads.size();
    //cout << num_final_reads << "\t" << final_reads.size() << "\t" << ratio << "\t" << searched_kmers.size() << endl;
    //num_final_reads = final_reads.size();
    num_total_edges += final_reads.size();
    cout << "final reads:  " << final_reads.size() << endl;
  }
  cout << "num total edges:  " << num_total_edges << endl;
  cout << "size of all reads touched:  " << all_reads_touched.size() << endl;
  //for(auto it_fr = final_reads.begin(); it_fr != final_reads.end(); ++ it_fr) {
  //  cout << it_fr->first << endl;
  //}
  //for(auto it_fr = final_reads.begin(); it_fr != final_reads.end(); ++ it_fr)  {
  //  cout << ">read_" << it_fr->first << endl;
  //  cout << sample_seqs_[it_fr->first] << endl; 
  //}
  //cout << "kmer tested:	" << searched_kmers.size() << endl;
  return 0;
}
