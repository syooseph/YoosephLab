#include "guided_assemble.h"
#include "extend_functor.h"
#include "certificate_functor.h"
#include "realignment_functor.h"

using namespace std;

GuidedAssemble::GuidedAssemble()  {
  info_loaded_ = false;
  return;
}

GuidedAssemble::GuidedAssemble(
    const std::string& in_working_directory,
    const std::string& in_result_directory,
    const std::string& in_query_sequence_file,
    const std::string& in_read_mfasta_file
)  {
  info_loaded_ = false;
  working_directory_ = in_working_directory;
  result_directory_ = in_result_directory;
  query_sequence_file_ = in_query_sequence_file;
  read_mfasta_file_ = in_read_mfasta_file;
  // setting default values
  num_threads_ = 1;
  //seed_len_ = 6;
  //alphabet_ = GBMR10;
  gap_open_ = -11;
  gap_extension_ = -1;
  scoring_matrix_ = BLOSUM62;
  right_band_ = 20; 
  down_band_ = 20;
  evalue_cutoff_ = 10;
  n_back_check_ = 10;
  dropoff_ = 25;
  write_certificate_ = false;
  return;
}

GuidedAssemble::GuidedAssemble(
      const std::string& in_working_directory,
      const std::string& in_result_directory,
      const std::string& in_query_sequence_file,
      const std::string& in_read_mfasta_file,
      const AssembleParameterType<AlignScoreType>& set_parameters
) {
  info_loaded_ = false;
  working_directory_ = in_working_directory;
  result_directory_ = in_result_directory;
  query_sequence_file_ = in_query_sequence_file;
  read_mfasta_file_ = in_read_mfasta_file;
  // setting default values
  num_threads_ = set_parameters.num_threads;
  //seed_len_ = set_parameters.seed_len;
  //alphabet_ = set_parameters.alphabet;
  gap_open_ = set_parameters.gap_open;
  gap_extension_ = set_parameters.gap_extension;
  scoring_matrix_ = set_parameters.scoring_matrix;
  right_band_ = ceil(set_parameters.band_size / 2); 
  down_band_ = ceil(set_parameters.band_size / 2);
  evalue_cutoff_ = set_parameters.evalue_cutoff;
  n_back_check_ = set_parameters.n_back_check;
  dropoff_ = set_parameters.dropoff;
  write_certificate_ = set_parameters.write_certificate;
  return;
}

GuidedAssemble::~GuidedAssemble() {
  if(info_loaded_)  {
    ClearForwardSeq();
    ClearReverseSeq();
    //ClearTag();
    ClearSuffixArrays();
  }
  return;
}

void GuidedAssemble::ClearSuffixArrays(void)  {
  suffix_array_->clear();
  reverse_suffix_array_->clear();
  delete suffix_array_;
  delete reverse_suffix_array_;
  return;
}


void GuidedAssemble::ClearForwardSeq(void)  {
  unsigned int i;
  for(i = 0; i < sample_num_reads_; ++ i) {
    delete [] sample_seqs_[i];
  }
  delete [] sample_seqs_;
  return;
}

void GuidedAssemble::ClearReverseSeq(void)  {
  unsigned int i;
  for(i = 0; i < sample_num_reads_; ++ i) {
    delete [] reversed_sample_seqs_[i];
  }
  delete [] reversed_sample_seqs_;
  return;
}

/*
void GuidedAssemble::ClearTag(void) {
  unsigned int i;
  for(i = 0; i < sample_num_reads_; ++ i) {
    delete [] sample_tags_[i];
  }
  delete [] sample_tags_;
  return;
}
*/

void GuidedAssemble::PrepareRunEnvironment(void)  {
  info_loaded_ = true;
  score_scheme_ = new ScoringFunction<AlignScoreType>(PROTEIN, scoring_matrix_, gap_extension_, gap_open_);
  PrepareWorkSpace();
  PrepareQuerySequence();
  PrepareSampleSequence();
  PrepareIndex();
  PrepareSuffixArray();
  PrepareReversedSuffixArray();
  return;
}

void GuidedAssemble::PrepareWorkSpace(void) {
  //if(working_directory_.empty())  {
  //  working_directory_ = "./";
  //}
  //if(*(-- working_directory_.end()) != '/')  {
  //  working_directory_ += "/";
  //}
  if(!boost::filesystem::exists(working_directory_))  {
    boost::filesystem::path dir(working_directory_);
    if(!boost::filesystem::create_directory(dir)) {
      cout << "GuidedAssemble::PrepareWorkSpace : Create working directory fails." << endl;
      exit(0);
    }  
  }
  if(!boost::filesystem::exists(result_directory_))  {
    boost::filesystem::path dir(result_directory_);
    if(!boost::filesystem::create_directory(dir)) {
      cout << "GuidedAssemble::PrepareWorkSpace : Create result directory fails." << endl;
      exit(0);
    }
  }
 
  return;
}

void GuidedAssemble::PrepareQuerySequence(void) {
  vector<string> files_in;
  files_in.resize(1);
  files_in[0] = query_sequence_file_;
  query_num_seqs_ = (unsigned int) seq::totalSequenceCount(files_in);
  char **tags = new char *[query_num_seqs_];
  char **seqs = new char *[query_num_seqs_];
  seq::loadSequences(files_in, tags, seqs, TAGSEQ);
  query_sequence_info_.resize(query_num_seqs_);
  query_sequence_.resize(query_num_seqs_);
  for(unsigned int i = 0; i < query_num_seqs_; ++ i) {
    query_sequence_info_[i] = string(tags[i]);
    query_sequence_[i] = string(seqs[i]);
    delete [] tags[i];
    delete [] seqs[i];
  }
  delete [] tags;
  delete [] seqs;
  return;
}

void GuidedAssemble::PrepareSampleSequence(void)  {
  vector<string> files_in;
  files_in.resize(1);
  files_in[0] = read_mfasta_file_;
  sample_num_reads_ = (unsigned int) seq::totalSequenceCount(files_in);
  //sample_tags_ = new char *[sample_num_reads_];
  sample_seqs_ = new char *[sample_num_reads_];
  //cerr << "num of reads:  " << sample_num_reads_ << endl;
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, sample_seqs_, SEQONLY);
  delete tag_foo;
  reversed_sample_seqs_ = new char *[sample_num_reads_];
  //cerr << "num of reads:  " << sample_num_reads_ << endl;
  seq::reverseSequences(sample_num_reads_, sample_seqs_, reversed_sample_seqs_);
  sample_size_ = 0.0;
  for(unsigned int i = 0; i < sample_num_reads_; ++ i) {
    sample_size_ += (double) strlen(sample_seqs_[i]) / 1000000.0;
  }
  return;
}

std::string GuidedAssemble::GetFileStem(const std::string& path)  {
  // going backward untill the '\/' character
  int i;
  for(i = path.length() - 1; i >= 0; -- i) {
    if(path[i] == '/')  {
      break;
    }
  }
  return path.substr(i + 1, path.length() - i - 1);
}

void GuidedAssemble::PrepareIndex(void) {
  double start_time = mytime();
  string read_stem = GetFileStem(read_mfasta_file_);
  
  string index_file_name = working_directory_ + '/' + read_stem + ".idx";
  string position_file_name = working_directory_ + '/' + read_stem + ".pos";
  //cerr << "index_file_name:	" << index_file_name << endl;
  /*
  IndexSample index_reads(seed_len_, alphabet_);
  if(!index_reads.IsCompatibleWithSetting(index_file_name, position_file_name))  {
    index_reads.WriteIndexFile(read_mfasta_file_, index_file_name, position_file_name);
    double write_end_time = mytime();
    printElapsed(start_time, write_end_time, "write indexing files");
    start_time = write_end_time;
  }
  */
  if(!boost::filesystem::exists(index_file_name) || !boost::filesystem::exists(position_file_name)) {
    cerr << "Indexing files not found: have you run \'Build\' first?" << endl;
    exit(0);    
  }
  // load positions for all kmers presented in the query sequence(s)
  IndexSample index_loader;
  index_loader.ReadKmerIndex(
      index_file_name, position_file_name, 
      query_sequence_, kmer_positions_
  );
  index_loader.BuiltParamLoader(index_file_name, position_file_name, alphabet_, seed_len_);
  cerr << "Alphabet and seed length:  " << alphabet_ << "  " << seed_len_  << endl;
  double load_end_time = mytime();
  printElapsed(start_time, load_end_time, "load indexing files");
  return;
}

void GuidedAssemble::PrepareSuffixArray(void)  {
  double start_time = mytime();
  string read_stem = GetFileStem(read_mfasta_file_);
  string sfa_file_name = working_directory_ + '/' + read_stem + ".sfa";
  string doc_file_name = working_directory_ + '/' + read_stem + ".doc";
  string lcp_file_name = working_directory_ + '/' + read_stem + ".lcp";
  string mcp_file_name = working_directory_ + '/' + read_stem + ".mcp";
  string gsa_file_name = working_directory_ + '/' + read_stem + ".gsa";
  if( //!boost::filesystem::exists(sfa_file_name) ||
     //!boost::filesystem::exists(doc_file_name) ||
     !boost::filesystem::exists(lcp_file_name) ||
     !boost::filesystem::exists(mcp_file_name) ||
     !boost::filesystem::exists(gsa_file_name))  {
    cerr << "Indexing files not found: have you run \'Build\' first?" << endl;
    exit(0);
    // create suffix arrays from scratch
    /*
    suffix_array_ = new GSA(sample_seqs_, sample_num_reads_, true);
    suffix_array_->dump(
        sfa_file_name.c_str(), doc_file_name.c_str(), 
        lcp_file_name.c_str(), mcp_file_name.c_str(), 
        gsa_file_name.c_str()
    );
    double build_end_time = mytime();
    printElapsed(start_time, build_end_time, "build suffix array");
    */
  } else  {
    // load the existing suffix array
    //cerr << "Forward suffix array built, loading indexing files..." << endl;
    suffix_array_ = new GSA();
    suffix_array_->load(
        //sfa_file_name.c_str(), doc_file_name.c_str(),
        lcp_file_name.c_str(), mcp_file_name.c_str(), gsa_file_name.c_str()
    );
    suffix_array_->setSequences(sample_seqs_);
    suffix_array_->setReadCount(sample_num_reads_);
    double load_end_time = mytime();
    printElapsed(start_time, load_end_time, "load suffix array");
  }
  return;
}

void GuidedAssemble::PrepareReversedSuffixArray(void)  {
  double start_time = mytime();
  string read_stem = GetFileStem(read_mfasta_file_);
  string sfa_file_name = working_directory_ + '/' + read_stem + ".re.sfa";
  string doc_file_name = working_directory_ + '/' + read_stem + ".re.doc";
  string lcp_file_name = working_directory_ + '/' + read_stem + ".re.lcp";
  string mcp_file_name = working_directory_ + '/' + read_stem + ".re.mcp";
  string gsa_file_name = working_directory_ + '/' + read_stem + ".re.gsa";
  if( //!boost::filesystem::exists(sfa_file_name) ||
     //!boost::filesystem::exists(doc_file_name) ||
     !boost::filesystem::exists(lcp_file_name) ||
     !boost::filesystem::exists(mcp_file_name) ||
     !boost::filesystem::exists(gsa_file_name))  {
    // create suffix arrays from scratch
    cerr << "Indexing files not found: have you run \'Build\' first?" << endl;
    exit(0);
    /*
    reverse_suffix_array_ = new GSA(reversed_sample_seqs_, sample_num_reads_, true);
    reverse_suffix_array_->dump(
        sfa_file_name.c_str(), doc_file_name.c_str(), 
        lcp_file_name.c_str(), mcp_file_name.c_str(), 
        gsa_file_name.c_str()
    );
    reverse_suffix_array_->setSequences(reversed_sample_seqs_);
    double build_end_time = mytime();
    printElapsed(start_time, build_end_time, "build reverse suffix array");
    */
  } else  {
    // load the existing suffix array
    //cerr << "Reverse suffix array built, loading indexing files..." << endl;
    reverse_suffix_array_ = new GSA();
    reverse_suffix_array_->load(
        //sfa_file_name.c_str(), doc_file_name.c_str(),
        lcp_file_name.c_str(), mcp_file_name.c_str(), gsa_file_name.c_str()
    );
    //reverse_suffix_array_->setSequences(reversed_sample_seqs_);
    reverse_suffix_array_->setSequences(reversed_sample_seqs_);
    //ClearReverseSeq();
    reverse_suffix_array_->setReadCount(sample_num_reads_);
    double load_end_time = mytime();
    printElapsed(start_time, load_end_time, "load reverse suffix array");
  }
  return;
}

void GuidedAssemble::DefineNonRepeatRegions(const string sequence, list<RegionType>& regions) {
  // this is a reserved void function, which will remove the low-complexity region from 
  // the sequence for seeding
  pair<int, int> single_region = {0, sequence.length() - 1};
  regions.push_back(single_region);
  return;
}

bool _cmp_kmer_info(KmerInfoType a, KmerInfoType b) {
  if(a.kmer_coverage >= b.kmer_coverage)  {
    return true;
  } else  {
    return false;
  }
}

void GuidedAssemble::DefineKmersToSearch(
    const string sequence, 
    const list<RegionType>& regions, 
    list<KmerInfoType>& kmers
) {
  //unordered_map<RIDType, bool> seed_reads;
  unsigned int i;
  for(auto it = regions.begin(); it != regions.end(); ++ it) {
    for(i = it->first; i <= it->second - seed_len_ + 1; ++ i) {
      string seed_str = sequence.substr(i, seed_len_);
      BoundType seed_range = suffix_array_->searchWithLCPs(
        (SfaChar*) seed_str.c_str(), seed_str.length()
      );
      KmerInfoType seed_info;
      seed_info.kmer_begin = i;
      seed_info.kmer_coverage = seed_range.second - seed_range.first + 1;
      kmers.push_back(seed_info);
      //for(int m = seed_range.first; m <= seed_range.second; ++ m) {
      //  seed_reads[suffix_array_->getId(m)] = true;
      //}
    }
  }
  /***
  cout << "size of seed reads:  " << seed_reads.size() << endl;
  for(auto it = seed_reads.begin(); it != seed_reads.end(); ++ it) {
    cout << (unsigned int) it->first << endl;
  }
  ***/
  kmers.sort(_cmp_kmer_info);
  return;
}


unsigned int GuidedAssemble::GetNumConcurrentThreads(void) {
  SLock r_lock(mutex_threads_);
  int n = num_running_threads_;
  r_lock.unlock();
  return n;
}

void GuidedAssemble::IncreaseNumThreads(void) {
  ULock w_lock(mutex_threads_);
  //cerr << "increasing:  " << num_running_threads_ << endl;
  ++ num_running_threads_;
  //cerr << "done increasing:  " << num_running_threads_ << endl;
  w_lock.unlock();
  return;
}

void GuidedAssemble::DecreaseNumThreads(void) {
  ULock w_lock(mutex_threads_);
  //cerr << "decreasing:  " << num_running_threads_ << endl;
  -- num_running_threads_;
  //cerr << "done decreasing:  " << num_running_threads_ << endl;
  w_lock.unlock();
  return;
}

void GuidedAssemble::PartitionQuery(
    const std::string& sequence, 
    map<int, int>& clump_map
) {
  // do the self-alignment
  SeqAlign<AlignScoreType> self_alignment(sequence, sequence, score_scheme_, LOCAL);
  self_alignment.Align();
  int cutoff_score = score_scheme_->ComputeRawScore(
      (int) (sample_size_ * 1000000 / sample_num_reads_),  
      (int) (sample_size_ * 1000000), evalue_cutoff_ * 1000
  );
  // get all the hsps
  list<vector<int> > hsp_range;
  self_alignment.ListAllHSPPositions(cutoff_score, hsp_range);
  // add the break point vector
  vector<int> break_points(sequence.length(), 0);
  for(auto it = hsp_range.begin(); it != hsp_range.end(); ++ it) {
    if((*it)[0] != (*it)[1])  { //  we don't care the alignment in the diagnoal
      int mid = (int) (((*it)[0] + (*it)[1]) / 2) + down_band_;
      if(mid > (int) sequence.length() - 1)  {
        mid = sequence.length() - 1;
      }
      break_points[mid] = 1;
    }
  }
  // construct the position-clump hash
  int i, clump_id = 0;
  for(i = 0; i < (int) sequence.length(); ++ i) {
    if(break_points[i] == 1)  {
      ++ clump_id;
    }    
    clump_map[i] = clump_id;
  }
  
  return;
}

bool _cmp_graph_critical_vertices(const pair<int, int>& item_a, const pair<int, int>& item_b) {
  if(item_a.second > item_b.second)  {
    return true;
  }
  return false;
}

bool _cmp_realignment(const CertificateType<AlignScoreType>& item_a, const CertificateType<AlignScoreType>& item_b) {
  if(item_a.alignment_score > item_b.alignment_score || 
      (item_a.alignment_score == item_b.alignment_score && item_a.mapped_locations.size() > item_b.mapped_locations.size()))  {
    return true;
  }
  return false;
}

AlignScoreType GuidedAssemble::EvalHamming(const string seq_a, const string seq_b)  {
  AlignScoreType sum_score = 0;
  int l = seq_a.length() > seq_b.length() ? seq_b.length() : seq_a.length();
  for(int i = 0; i < l; ++ i) {
    sum_score += score_scheme_->CheckMatchScore(seq_a[i], seq_b[i]);
  }
  return sum_score;
}

void GuidedAssemble::AssembleAndAlign(unordered_map<RIDType, ReadMapType>& mapped_reads) {
  // load k-mer positions from indexing file
  unsigned int i;
  IndexSample index_caller(seed_len_, alphabet_);
  
  for(i = 0; i < query_sequence_.size(); ++ i) {
    double start_time = mytime();
    map<int, int> clump_map;
    PartitionQuery(query_sequence_[i], clump_map);
    double partition_end_time = mytime();
    printElapsed(start_time, partition_end_time, "partition query");
    
    start_time = mytime();    
    int num_seeds = 0;
    list<RegionType> regions_to_seed;
    DefineNonRepeatRegions(query_sequence_[i], regions_to_seed);
    list<KmerInfoType> seed_kmers; // the first int is the starting location of the seed
        // in the query sequence, the second int is the coverage of the seed
    DefineKmersToSearch(query_sequence_[i], regions_to_seed, seed_kmers);
    unordered_map<unsigned int, bool> covered_region;
    map<int, unordered_map<int, set<string> > > seed_pool;
    for(auto it_kmerinfo = seed_kmers.begin(); it_kmerinfo != seed_kmers.end(); ++ it_kmerinfo) {
      // getting the k-mer sequence in the query
      string current_kmer = query_sequence_[i].substr(it_kmerinfo->kmer_begin, seed_len_);
      KmerType encoded_kmer = index_caller.encode_kmer(current_kmer);
      // getting positions of k-mer sequences in the sample
      list<PositionType> seed_candidates = kmer_positions_[encoded_kmer];
      unordered_map<string, int> translated_kmers;
      for(auto it = seed_candidates.begin(); it != seed_candidates.end(); ++ it) {
        string seed_seq = sample_seqs_[(int) it->rid];
        string sample_kmer = seed_seq.substr(it->pos, seed_len_);
        translated_kmers[sample_kmer] = 1;
      }
      for(auto it = translated_kmers.begin(); it != translated_kmers.end(); ++ it) {
        AlignScoreType est_score = EvalHamming(current_kmer, it->first);
        // if the alignment score is high enough
        if(est_score >= (AlignScoreType) (0.6 * score_scheme_->GetAveMatch() * seed_len_))  {
          //cout << it->first << endl;
          ++ num_seeds;
          seed_pool[static_cast<int>(est_score)][it_kmerinfo->kmer_begin].insert(it->first);
        }
      }
    }
    double build_end_time = mytime();
    printElapsed(start_time, build_end_time, "seeds selection");
    //cout << query_sequence_[0] << endl;
    //exit(0);
    
    start_time = mytime();
    vector<AssemblyGraph*> graph_holder(num_seeds);
    vector<list<VertexPairType>* > seed_vertices_holder(num_seeds);
    vector<bool> *space_indicator = new vector<bool>(num_seeds);
    boost::threadpool::pool extension_pool(num_threads_);
    unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> > *rough_mapped_region = 
        new unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >;
    LockType* mutex_mapped_region = new LockType;
    LockType* mutex_output = new LockType;
    LockType* mutex_indicator = new LockType;
    int phase_counter = 0;
    vector<int> begin_locs;
    int q_index = i;
    for(auto it = seed_pool.rbegin(); it != seed_pool.rend(); ++ it) {
      // for each score, from high to low
      for(auto it_s = it->second.begin(); it_s != it->second.end(); ++ it_s) {
        // for each position in the query sequence
        for(auto it_sp = it_s->second.begin(); it_sp != it_s->second.end(); ++ it_sp) {          
          graph_holder[phase_counter] = new AssemblyGraph;
          seed_vertices_holder[phase_counter] = new list<VertexPairType>;
          //cout << i << "  " << q_index << endl;
          boost::shared_ptr<ExtendFunctor> job(
            new ExtendFunctor(
                phase_counter, this, q_index, it_s->first, *it_sp, &clump_map,
                graph_holder[phase_counter], seed_vertices_holder[phase_counter], 
                rough_mapped_region, mutex_mapped_region, &mapped_reads, mutex_output,
                space_indicator, mutex_indicator 
            )
          );
          boost::threadpool::schedule(extension_pool, boost::bind(&ExtendFunctor::run, job));
          ++ phase_counter;
          begin_locs.push_back(it_s->first);
        }
      }
    }
    extension_pool.wait();
    double extend_end_time = mytime();
    printElapsed(start_time, extend_end_time, "seed extension");
    cout << "Num traversed reads: " << rough_mapped_region->size() << endl;
    delete rough_mapped_region;
    delete mutex_mapped_region;
    
    
    start_time = mytime();
    boost::threadpool::pool certificate_pool(num_threads_);
    list<AssembledPathType> path_holder;
    LockType* mutex_path_holder = new LockType;  
    for(unsigned int index = 0; index < space_indicator->size(); ++ index) {
      if((*space_indicator)[index])  {
        set<RIDType> foo_set;
        boost::shared_ptr<CertificateFunctor> job(
          new CertificateFunctor(
              this, i, graph_holder[index], seed_vertices_holder[index],
              &foo_set, &path_holder, mutex_path_holder, &mapped_reads, mutex_output        
          )
        );
        boost::threadpool::schedule(certificate_pool, boost::bind(&CertificateFunctor::run, job));
      }
    }
    certificate_pool.wait();
    double reconstruct_end_time = mytime();
    printElapsed(start_time, reconstruct_end_time, "path reconstruction");
    
    //exit(0);
    
    
    // collect memory
    for(int idx = 0; idx < num_seeds; ++ idx) {
      if((*space_indicator)[idx])  {
        delete graph_holder[idx];
        delete seed_vertices_holder[idx];
      }
    }
    delete space_indicator;
    delete mutex_output;
    delete mutex_path_holder;
    delete mutex_indicator;
    
    WriteRecruitedReads(q_index, mapped_reads);
    //cout << "end of loop: " << i << endl;

    if(!write_certificate_)  {
      continue;
    }
    //ClearSuffixArrays();
    
    //continue;
    
    start_time = mytime();
    list<CertificateType<AlignScoreType> > realignment_holder;
    LockType* mutex_realignment_holder = new LockType;
    boost::threadpool::pool realignment_pool(num_threads_);
    int counter = 0;
    for(auto it = path_holder.begin(); it != path_holder.end(); ++ it) {
      if(it->assembled_sequence.length() < 60)  {
        continue;
      }
      if(it->mapped_begin < 0)  {
        it->mapped_begin = 0;
      }
      if(it->mapped_end >= (int) query_sequence_[i].length())  {
        it->mapped_end = query_sequence_[i].length() - 1;
      }
      string query_partial = query_sequence_[i].substr(it->mapped_begin, it->mapped_end - it->mapped_begin + 1);
      boost::shared_ptr<RealignmentFunctor> job(
        new RealignmentFunctor(
            this, query_partial, it->assembled_sequence,
            &(*it), &realignment_holder, mutex_realignment_holder         
        )
      ); 
      boost::threadpool::schedule(realignment_pool, boost::bind(&RealignmentFunctor::run, job));
      cout << ">contig_" << counter << endl << it->assembled_sequence << endl;
      ++ counter;
    }
    realignment_pool.wait();
    double realignment_end_time = mytime();
    printElapsed(start_time, realignment_end_time, "path realignment");
    
    realignment_holder.sort(_cmp_realignment);
    //cout << "Search of sequence:  " << query_sequence_info_[i] << endl << endl;
    //SelectRecruitedReads(i, query_sequence_[i], realignment_holder, mapped_reads);
    WriteAlignments(q_index, query_sequence_[q_index], realignment_holder);
    delete mutex_realignment_holder;   
  }
  delete score_scheme_;
  WriteLog();
  return;
}

void GuidedAssemble::WriteLog(void)  {
  boost::filesystem::path abs_wdir = boost::filesystem::absolute(working_directory_);
  boost::filesystem::path abs_rdir = boost::filesystem::absolute(result_directory_);
  boost::filesystem::path abs_query = boost::filesystem::absolute(query_sequence_file_);
  boost::filesystem::path abs_target = boost::filesystem::absolute(read_mfasta_file_);
  string file_stem = GetFileStem(query_sequence_file_);
  string reads_file = result_directory_ + "/" + file_stem + ".log"; 
  ofstream out_info(reads_file.c_str(), ios_base::out);
  out_info << "Parameters used for the search task:" << endl << endl;
  out_info << "Work directory:  " << abs_wdir.string() << endl;
  out_info << "Result directory:  " << abs_rdir.string() << endl;
  out_info << "Query sequence:  " << abs_query.string() << endl;
  out_info << "Target reads:  " << abs_target.string() << endl;
  // setting default values
  out_info << "Number of threads:  " << num_threads_ << endl;
  out_info << "Seed length:  " << seed_len_ << endl;
  //alphabet_ = set_parameters.alphabet;
  out_info << "Alphabet in use:  ";
  switch(alphabet_)  {
    case 0:
      out_info << "ALL20" << endl;
      break;
    case 1:
      out_info << "DSSP5" << endl;
      break;
    case 2:
      out_info << "DSSP10" << endl;
      break;
    case 3:
      out_info << "GBMR4" << endl;
      break;
    case 4:
      out_info << "GBMR10" << endl;
      break;
    case 5:
      out_info << "HSDM5" << endl;
      break;
    case 6:
      out_info << "SDM6" << endl;
      break;
    case 7:
      out_info << "MURPHY5" << endl;
      break;
    case 8:
      out_info << "MURPHY10" << endl;
      break;
    case 9:
      out_info << "TD5" << endl;
      break;
    case 10:
      out_info << "TD10" << endl;
      break;
    default:
      out_info << "GBMR10" << endl;
  }
  out_info << "Amino-acids substitution matrix in use:  ";
  switch(scoring_matrix_)  {
    case 0:
      out_info << "BLOSUM62" << endl;
      break;
    case 1:
      out_info << "BLOSUM80" << endl;
      break;
    case 2:
      out_info << "BLOSUM90" << endl;
      break;
    case 3:
      out_info << "BLOSUM50" << endl;
      break;
    case 4:
      out_info << "BLOSUM45" << endl;
      break;
    case 5:
      out_info << "PAM250" << endl;
      break;
    case 6:
      out_info << "PAM70" << endl;
      break;
    case 7:
      out_info << "PAM30" << endl;
      break;
    default:
      out_info << "BLOSUM62" << endl;
  }
  out_info << "Gap open penalty:  " << gap_open_ << endl;
  out_info << "Gap extension penalty:  " << gap_extension_ << endl;
  out_info << "Alignment band size:  " << down_band_ * 2 << endl;
  out_info << "E-value cutoff:  " << evalue_cutoff_ << endl;
  out_info << "N_back size (path-read overlap):  " << n_back_check_ << endl;
  out_info << "Write Certificate:  " << write_certificate_ << endl;
  out_info.close();
  return;
}

int GuidedAssemble::NumGaps(const std::string& seq) {
  int count = 0;
  for(unsigned int i = 0; i < seq.length(); ++ i) {
    if(seq[i] == '-')  {
      ++ count;
    }
  }
  return count;
}

void GuidedAssemble::WriteRecruitedReads(
    const int& query_ID, std::unordered_map<RIDType, ReadMapType> mapped_reads
)  {
  string file_stem = GetFileStem(query_sequence_file_);
  ostringstream convert;
  convert << query_ID;
  string reads_file = result_directory_ + "/" + file_stem + ".reads";
  //cout << alignment_file << endl;
  
  if(!boost::filesystem::exists(reads_file) || query_ID == 0)  {
    ofstream out_alignment(reads_file.c_str(), ios_base::out);
    out_alignment << "# Fields: query_sequence_ID, read_ID, best_Evalue_aligned" << endl;
    for(auto it = mapped_reads.begin(); it != mapped_reads.end(); ++ it) {
      out_alignment << query_ID << "  " << it->first << " " << it->second.best_evalue << endl;
    }
    out_alignment.close();
  }  else  {
    ofstream out_alignment(reads_file.c_str(), ios_base::app);
    for(auto it = mapped_reads.begin(); it != mapped_reads.end(); ++ it) {
      out_alignment << query_ID << "  " << it->first << " " << it->second.best_evalue << endl;
    }
    out_alignment.close();
  }
  return;
}

void GuidedAssemble::WriteAlignments(
    const int& query_ID, const std::string& query_sequence, 
    const std::list<CertificateType<AlignScoreType> >& certificates
) {
  string file_stem = GetFileStem(query_sequence_file_);
  ostringstream convert;
  convert << query_ID;
  string alignment_file = result_directory_ + "/" + file_stem + ".aln";
  //cout << alignment_file << endl;
  ofstream *out_alignment;
  if(!boost::filesystem::exists(alignment_file) || query_ID == 0)  {
    out_alignment = new ofstream(alignment_file.c_str(), ios_base::out);
  }  else  {
    out_alignment = new ofstream(alignment_file.c_str(), ios_base::app);
  }
  for(auto it_cer = certificates.begin(); it_cer != certificates.end(); ++ it_cer) {
    
    double bit_score = score_scheme_->ComputeBitScore(it_cer->alignment_score);
    double e_value = score_scheme_->ComputeEValue(query_sequence.length(), sample_size_ * 1000000, it_cer->alignment_score);
    if(e_value > evalue_cutoff_)  {
      break;
    }
    *out_alignment << "****************** ALIGNMENT BEGINS ******************" << endl << endl;
    *out_alignment << "Score: " << bit_score << " (" << it_cer->alignment_score << ")  " << "Evalue: " << e_value << endl << endl;
    int query_begin = it_cer->alignment_regions[0] + 1; // use the '1'-based system
    for(unsigned int i = 0; i < it_cer->alignment_query.length(); i += 60)	{
      string query_frag = it_cer->alignment_query.substr(i, 60);
      int query_end = query_begin + query_frag.length() - NumGaps(query_frag) - 1;
        *out_alignment << "  " << query_begin << "\t" << it_cer->alignment_query.substr(i, 60) << "  " << query_end << std::endl;
        *out_alignment << "\t" << it_cer->alignment_symbol.substr(i, 60) << std::endl;
        *out_alignment << "\t" << it_cer->alignment_target.substr(i, 60) << std::endl;
        *out_alignment << std::endl;
        query_begin = query_end + 1;
      }	
    *out_alignment << "included reads: (read_ID,begin_location_in_assembled_sequence;)" << endl;
    for(auto it = it_cer->mapped_locations.begin(); it != it_cer->mapped_locations.end(); ++ it) {
      *out_alignment << (unsigned int) it->first << "," << (unsigned int) it->second << ";";
    }
    *out_alignment << std::endl;
    *out_alignment << endl;
  }
  out_alignment->close();
  return;
}

void GuidedAssemble::SelectRecruitedReads(
    const int& query_ID, const string& query_sequence, const std::list<CertificateType<AlignScoreType> >& certificates,
    std::unordered_map<RIDType, ReadMapType>& mapped_reads
)  {
  for(auto it_cer = certificates.begin(); it_cer != certificates.end(); ++ it_cer) {
    double e_value = score_scheme_->ComputeEValue(query_sequence.length(), sample_size_ * 1000000, it_cer->alignment_score);	
    if(e_value <= evalue_cutoff_)  {
	    for(auto it = it_cer->mapped_locations.begin(); it != it_cer->mapped_locations.end(); ++ it) {
	      auto it_read = mapped_reads.find(it->first);
	      if(it_read != mapped_reads.end() && it_read->second.best_evalue > e_value)  {
	        it_read->second.query_seq_ID = query_ID;
	        it_read->second.best_evalue = e_value;
	      } else  {
	        ReadMapType r_info;
	        r_info.query_seq_ID = query_ID;
	        r_info.best_evalue = e_value;
	        mapped_reads[it->first] = r_info;
	      }
	    }
    }
  }
  return;
}

void GuidedAssemble::WriteAlignment(const string& query_sequence, const CertificateType<AlignScoreType>& certificate) {
  cout << "****************** ALIGNMENT BEGINS ******************" << endl << endl;
  double bit_score = score_scheme_->ComputeBitScore(certificate.alignment_score);
  double e_value = score_scheme_->ComputeEValue(query_sequence.length(), sample_size_ * 1000000, certificate.alignment_score);
  cout << "Score: " << bit_score << " (" << certificate.alignment_score << ")  " << "Evalue: " << e_value << endl << endl;
  int query_begin = certificate.alignment_regions[0] + 1; // use the '1'-based system
  for(unsigned int i = 0; i < certificate.alignment_query.length(); i += 60)	{
    string query_frag = certificate.alignment_query.substr(i, 60);
    int query_end = query_begin + query_frag.length() - NumGaps(query_frag) - 1;
		std::cout << "  " << query_begin << "\t" << certificate.alignment_query.substr(i, 60) << "  " << query_end << std::endl;
		std::cout << "\t" << certificate.alignment_symbol.substr(i, 60) << std::endl;
		std::cout << "\t" << certificate.alignment_target.substr(i, 60) << std::endl;
		std::cout << std::endl;
		query_begin = query_end + 1;
	}
	
	
	std::cout << "included reads: " << endl;
	for(auto it = certificate.mapped_locations.begin(); it != certificate.mapped_locations.end(); ++ it) {
	  std::cout << (unsigned int) it->first << "," << (unsigned int) it->second << ";";
	}
	std::cout << std::endl;
  cout << endl;
	/*
	// module for checking sequence compatibility with the reads
	string all_seq = "";
	for(unsigned int i = 0; i < certificate.alignment_target.length(); ++ i) {
	  if(certificate.alignment_target[i] != '-')  {
	    all_seq += certificate.alignment_target[i];
	  }
	}  

	for(auto it = certificate.mapped_locations.begin(); it != certificate.mapped_locations.end(); ++ it) {
	  //std::cout << (unsigned int) it->first << "," << (unsigned int) it->second << ";";
	  string r_seq = string(sample_seqs_[it->first]);
	  string a_seq = all_seq.substr(it->second, r_seq.length());
	  if(r_seq != a_seq)  {
	    std::cerr << "*******************" << std::endl;
	    cerr << r_seq << endl;
	    cerr << a_seq << endl;
	  }
	}
	*/
	
  return;
}

/*
void GuidedAssemble::GetReadTags(std::vector<std::string>& tags)  {
  tags.resize(sample_num_reads_);
  for(unsigned int i = 0; i < sample_num_reads_; ++ i) {
    tags[i] = string(sample_tags_[i]);
  }
  return;
}
*/

int GuidedAssemble::ReconcileReads(
    const std::string& query_sequence,
    const AlignScoreType& cutoff,
    const AlignScoreType& opposite_max_score,
    const AlignScoreType& opposite_min_score,
    std::unordered_map<RIDType, ReadVertexType>& recorded_reads,
    std::unordered_map<RIDType, ReadMapType>& mapped_reads 
) {
  // adding all reads to the hash table
  int num_new_reads = 0;
  for(auto it = recorded_reads.begin(); it != recorded_reads.end(); ++ it) {
    //cout << "scores:  " << cutoff << "  " << it->second.max_score << "  " << opposite_max_score << endl;
    if(it->second.max_score + opposite_max_score >= cutoff)  {
      // record the read
      double best_evalue = score_scheme_->ComputeEValue(
          query_sequence.length(), sample_size_ * 1000000, 
          it->second.max_score + opposite_max_score
      );
      auto mapped_read_iter = mapped_reads.find(it->first);
      if(mapped_read_iter == mapped_reads.end())  {
        ReadMapType score_info;
        //score_info.is_bridging_read = false;
        //score_info.max_vertex_source = (BoostVertex) it->second.max_vertex_source;
        //score_info.max_vertex_target = (BoostVertex) it->second.max_vertex_target;
        score_info.best_evalue = best_evalue;
        mapped_reads[it->first] = score_info;
        ++ num_new_reads;
      } else  {
        if(best_evalue < mapped_read_iter->second.best_evalue)  {
          //mapped_read_iter->second.is_bridging_read = false;
          //mapped_read_iter->second.max_vertex_source = (BoostVertex) it->second.max_vertex_source;
          //mapped_read_iter->second.max_vertex_target = (BoostVertex) it->second.max_vertex_target;
          mapped_read_iter->second.best_evalue = best_evalue;
        }
      }
    }
  }
  return num_new_reads;
}

int GuidedAssemble::ReconcileBridgingReads(
    const std::string& query_sequence, 
    const AlignScoreType& cutoff,
    const AlignScoreType& left_max_score,
    const AlignScoreType& left_min_score,
    const AlignScoreType& right_max_score,
    const AlignScoreType& right_min_score,
    VertexPairType seed_pair,
    std::unordered_map<RIDType, ReadMapType>& mapped_reads 
) {
  int num_new_reads = 0;
  AlignScoreType seed_score = seed_pair.seed_match_score + seed_pair.left_score + seed_pair.right_score;
  for(auto it = seed_pair.bridging_reads.begin(); it != seed_pair.bridging_reads.end(); ++ it) {
    if(seed_score + left_max_score + right_max_score >= cutoff)  {
      // record the read
      double best_evalue = score_scheme_->ComputeEValue(
          query_sequence.length(), sample_size_ * 1000000, 
          seed_score + left_max_score + right_max_score
      );
      auto mapped_read_iter = mapped_reads.find(it->read_ID);
      if(mapped_read_iter == mapped_reads.end())  {
        ReadMapType score_info;
        //score_info.is_bridging_read = true; // do not need the best edge information
        score_info.best_evalue = best_evalue;
        mapped_reads[it->read_ID] = score_info;
        ++ num_new_reads;
      } else  {
        if(best_evalue < mapped_read_iter->second.best_evalue)  {
          //mapped_read_iter->second.is_bridging_read = true;
          mapped_read_iter->second.best_evalue = best_evalue;
        }
      }
    }
  }
  return num_new_reads;
}

void GuidedAssemble::TraceBack(
    const std::string& query_sequence,
    AssemblyGraph& graph, std::list<VertexPairType>& source_vertices, 
    unordered_map<RIDType, ReadMapType>& mapped_reads
)  {
  //cout << "***********TraceBack*************************" << endl;
  AlignScoreType score_cutoff = score_scheme_->ComputeRawScore(
      query_sequence.length(), sample_size_ * 1000000, evalue_cutoff_
  );
  for(auto it = source_vertices.begin(); it != source_vertices.end(); ++ it) {
    //list<EdgeProperty> assembled_seqs_left, assembled_seqs_right;
    unordered_map<RIDType, ReadVertexType> recorded_reads_left, recorded_reads_right;
    AlignScoreType left_max_score = 0, left_min_score = 0, right_max_score = 0, right_min_score = 0;
    AlignScoreType seed_score = it->seed_match_score + it->left_score +it->right_score;
    if(it->has_initialized_right)  {
      graph.TraverseGraphUnder(it->vertex_right, recorded_reads_right);
      right_max_score = graph.AccessVertex(it->vertex_right).max_score_below + seed_score;
    }
    if(it->has_initialized_left)  {
      graph.TraverseGraphUnder(it->vertex_left, recorded_reads_left);
      left_max_score = graph.AccessVertex(it->vertex_left).max_score_below + seed_score;
    }
    //cout << "Traversed reads left: " << recorded_reads_left.size() << endl;
    //cout << "Traversed reads right: " << recorded_reads_right.size() << endl;
    
    ReconcileBridgingReads(
        query_sequence, score_cutoff, 
        left_max_score, left_min_score, right_max_score, right_min_score, 
        *it, mapped_reads
    );
    ReconcileReads(
        query_sequence, score_cutoff, 
        left_max_score, left_min_score, 
        recorded_reads_right, mapped_reads
    );
    ReconcileReads(
        query_sequence, score_cutoff, 
        right_max_score, right_min_score, 
        recorded_reads_left, mapped_reads
    );
  }
  //cout << "Total included reads:  " << mapped_reads.size() << endl;
  return;
}

void GuidedAssemble::FetchSeeds(const std::string& sequence, std::set<std::string>& matched_seeds)  {
  IndexSample index_caller(seed_len_, alphabet_);
  unordered_map<RIDType, bool> seed_reads;
  for(unsigned int i = 0; i <= sequence.length() - seed_len_; ++ i) {
    // getting the k-mer sequence in the query
    string current_kmer = sequence.substr(i, seed_len_);
    KmerType encoded_kmer = index_caller.encode_kmer(current_kmer);
    // getting positions of k-mer sequences in the sample
    list<PositionType> seed_candidates = kmer_positions_[encoded_kmer];
    unordered_map<string, int> translated_kmers;
    for(auto it = seed_candidates.begin(); it != seed_candidates.end(); ++ it) {
      string seed_seq = sample_seqs_[(int) it->rid];
      string sample_kmer = seed_seq.substr(it->pos, seed_len_);
      translated_kmers[sample_kmer] = 1;
    }
    for(auto it = translated_kmers.begin(); it != translated_kmers.end(); ++ it) {
      AlignScoreType est_score = EvalHamming(current_kmer, it->first);
      // if the alignment score is high enough
      if(est_score >= (AlignScoreType) (0.6 * score_scheme_->GetAveMatch() * seed_len_))  {
        matched_seeds.insert(it->first);  
      }
    }
  }
  return;
}

void GuidedAssemble::PreFilterAll(void)  {
  set<string> matched_seeds;
  FetchSeeds(query_sequence_[0], matched_seeds);
  set<RIDType> candidate_reads;
  for(auto it = matched_seeds.begin(); it != matched_seeds.end(); ++ it) {
    set<RIDType> reached_reads;
    // progressively extend to neighboring reads
    ProgressiveSearch(*it, 3, 10, reached_reads);
    // copy the reached reads as candidate reads
    for(auto it_r = reached_reads.begin(); it_r != reached_reads.end(); ++ it_r) {
      candidate_reads.insert(*it_r);
    }
  }
  for(auto it = candidate_reads.begin(); it != candidate_reads.end(); ++ it) {
    cout << *it << endl;
  }
  return;
}

void GuidedAssemble::ProgressiveSearch(
    const std::string& seed, int step_forward, int branch_cutoff, 
    std::set<RIDType>& reached_reads
)  {
  // forward search
  set<string> fw_seeds;
  fw_seeds.insert(seed);
  for(int i = 0; i < step_forward; ++ i) {
    // search
    set<RIDType> c_reads;
    SearchSingleFw(fw_seeds, c_reads);
    // record reads
    for(auto it = c_reads.begin(); it != c_reads.end(); ++ it) {
      reached_reads.insert(*it);
    }
    // update
    fw_seeds.clear();
    DefineNewSearchSeqFw(c_reads, fw_seeds);
  }
  // reverse search
  set<string> re_seeds;
  string r_seed = string(seed.rbegin(), seed.rend());  
  re_seeds.insert(r_seed);
  for(int i = 0; i < step_forward; ++ i) {
    // search
    set<RIDType> c_reads;
    SearchSingleRe(re_seeds, c_reads);
    // record reads
    for(auto it = c_reads.begin(); it != c_reads.end(); ++ it) {
      reached_reads.insert(*it);
    }
    // update
    re_seeds.clear();
    DefineNewSearchSeqRe(c_reads, re_seeds);
  }
  return;
}

void GuidedAssemble::SearchSingleFw(std::set<std::string>& search_seed, std::set<RIDType>& containing_reads)  {
  for(auto it = search_seed.begin(); it != search_seed.end(); ++ it) {
    pair<SfaType, SfaType> range = suffix_array_->searchWithLCPs((SfaChar*) it->c_str(), it->length());
    for(int i = range.first; i <= range.second; ++ i) {
      containing_reads.insert(suffix_array_->getId(i));
    }
  }
  return;
}

void GuidedAssemble::SearchSingleRe(std::set<std::string>& search_seed, std::set<RIDType>& containing_reads)  {
  for(auto it = search_seed.begin(); it != search_seed.end(); ++ it) {
    pair<SfaType, SfaType> range = reverse_suffix_array_->searchWithLCPs((SfaChar*) it->c_str(), it->length());
    for(int i = range.first; i <= range.second; ++ i) {
      containing_reads.insert(reverse_suffix_array_->getId(i));
    }
  }
  return;
}

void GuidedAssemble::DefineNewSearchSeqFw(std::set<RIDType>& reads, std::set<std::string>& end_sequence)  {
  for(auto it = reads.begin(); it != reads.end(); ++ it) {
    string es = sample_seqs_[*it];
    end_sequence.insert(es.substr(es.length() - n_back_check_, n_back_check_));
  }
  return;
}

void GuidedAssemble::DefineNewSearchSeqRe(std::set<RIDType>& reads, std::set<std::string>& end_sequence)  {
  for(auto it = reads.begin(); it != reads.end(); ++ it) {
    string es = reversed_sample_seqs_[*it];
    end_sequence.insert(es.substr(es.length() - n_back_check_, n_back_check_));
  }
  return;
}
