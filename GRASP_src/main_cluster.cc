#include "read_cluster.h"
#include "seq_align.h"
#include "timer.h"

typedef int AlignScoreType;

using namespace std;


int main()  {
  string index_stem = "./WorkSpace/saliva.faa";
  string file_name = "./saliva.faa";
  string query_file = "./WorkSpace/SGO_0064.fa";
  int n_back = 10;
  int seed_len = 6;
  double start = mytime();
  string lcp_file = index_stem + ".lcp";
  string mcp_file = index_stem + ".mcp";
  string gsa_file = index_stem + ".gsa";
  string idx_file = index_stem + ".idx";
  string pos_file = index_stem + ".pos";
  // define partition of reads
  ReadCluster read_partition(n_back, file_name, lcp_file, mcp_file, gsa_file);
  //cout << "Done loading" << endl;
  read_partition.RecruitConnectedReads();
  //cout << "Done traversing suffix array and define connectivity" << endl;
  double end_con = mytime();
  //printElapsed(start, end_con, "Define connectivity");
  read_partition.InterpretConnections();
  //cout << "Done partition kmers" << endl;
  double end_partition = mytime();
  //printElapsed(end_con, end_partition, "Building kmer partition");
  // load query sequence
  vector<string> files_in;
  files_in.push_back(query_file);
  int query_num_reads = (unsigned int) seq::totalSequenceCount(files_in);
  char **query_seqs = new char *[query_num_reads];
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, query_seqs, SEQONLY);
  delete tag_foo;
  string query = query_seqs[0];
  // load kmer index
  unordered_map<KmerType, list<PositionType> > kmer_positions;
  IndexSample index_loader;
  index_loader.ReadKmerIndex(
      idx_file, pos_file, query, kmer_positions
  );
  double end_load_kmer = mytime();
  //printElapsed(end_partition, end_load_kmer, "Load kmer index");
  ScoringFunction<AlignScoreType> *score_scheme_ = new ScoringFunction<AlignScoreType>(PROTEIN, BLOSUM62, -1, -11);
  unordered_map<KmerType, bool> recruited_n_back_mer;

  unordered_map<int, bool> tested_reads;
  list<string> kmer_repo;
  IndexSample encoder(n_back, ALL20);
  IndexSample seed_encoder(seed_len, ALL20);
  int i;
  for(i = 0; i < query.length() - seed_len + 1; ++ i) {
    string current_kmer = query.substr(i, seed_len);
    KmerType encoded_kmer = index_loader.encode_kmer(current_kmer);
    list<PositionType> seed_candidates = kmer_positions[encoded_kmer];
    unordered_map<KmerType, bool> tested_seed_kmer;
    for(auto it = seed_candidates.begin(); it != seed_candidates.end(); ++ it) {
      if(tested_reads.find(it->rid) != tested_reads.end())  {
        continue;
      }    
      string seed_seq = read_partition.GetReadSequence((int) it->rid);
      string sample_kmer = seed_seq.substr(it->pos, seed_len);
      KmerType seed_encode = seed_encoder.encode_kmer(sample_kmer);
      if(tested_seed_kmer.find(seed_encode) != tested_seed_kmer.end())  {
        continue;
      }
      
      AlignScoreType est_score = 0;
      for(int i2 = 0; i2 < current_kmer.length(); ++ i2) {
        est_score += score_scheme_->CheckMatchScore(current_kmer[i2], sample_kmer[i2]);
      }
      // if the alignment score is high enough
      if(est_score >= (AlignScoreType) (0.6 * score_scheme_->GetAveMatch() * seed_len))  {
        tested_reads[it->rid] = true;
        string fw = seed_seq.substr(0, n_back);
        string re = seed_seq.substr(seed_seq.length() - n_back, n_back);
        KmerType fw_encode = encoder.encode_kmer(fw);
        KmerType re_encode = encoder.encode_kmer(re);
        if(recruited_n_back_mer.find(fw_encode) == recruited_n_back_mer.end())  {
          recruited_n_back_mer[fw_encode] = true;
          kmer_repo.push_back(fw);
        }
        if(recruited_n_back_mer.find(re_encode) == recruited_n_back_mer.end())  {
          recruited_n_back_mer[re_encode] = true;
          kmer_repo.push_back(re);
        }
      } else  {
        tested_seed_kmer[seed_encode] = true;
      }
    }
  }
  double end_build_seed = mytime();
  //printElapsed(end_load_kmer, end_build_seed, "Building seed kmer to search");
  read_partition.GetClosureSequences(kmer_repo);
  double end_print = mytime();
  //printElapsed(end_build_seed, end_print, "Retriving sequences and write to file");
  return 0;
  
}
