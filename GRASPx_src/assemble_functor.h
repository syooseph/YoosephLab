#include "sequence_build.h"
#include "reduced_alphabet.h"
#include "database_index.h"
#include "reachable_reads.h"
#include "scoring_function.h"
#include "read_alignment.h"
#include "greedy_assembly.h"
#include "assemble_extend.h"
#include "contig_refinement.h"

#include "sequence.h"
#include "timer.h"

#include <boost/thread.hpp>
#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/smart_ptr.hpp>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <cassert>
#include <map>
#include <unordered_map>
#include <stack>

typedef boost::shared_mutex LockType;
typedef boost::shared_lock<LockType> SLock;
typedef boost::unique_lock<LockType> ULock;

class AssembleFunctor {
 public:
  AssembleFunctor(
      int thread_id, bool verbose, 
      LockType *mutex_out, LockType *mutex_seq, LockType *mutex_score, LockType *mutex_link,  
      std::string query_tag, std::string& query_seq,
      SequenceBuild* seq_obj, ScoringFunction<int>* score_obj, ReachableReads* link_obj,
      double seed_score_scale, int assembly_depth, int band_size, double e_value,
      std::list<ContigType>* contigs
  );
  ~AssembleFunctor(void);
  void Run(void);
 private:
  int thread_id_;
  bool verbose_;
  LockType* mutex_out_;
  LockType* mutex_seq_;
  LockType* mutex_score_;
  LockType* mutex_link_;
  std::string query_tag_;
  std::string query_seq_;
  SequenceBuild* seq_obj_;
  ScoringFunction<int>* score_obj_;
  ReachableReads* link_obj_;
  double seed_score_scale_;
  int assembly_depth_;
  int band_size_;
  double e_value_;
  std::list<ContigType>* contigs_; 
};

AssembleFunctor::AssembleFunctor(
    int thread_id, bool verbose, 
    LockType *mutex_out, LockType *mutex_seq, LockType *mutex_score, LockType *mutex_link,  
    std::string query_tag, std::string& query_seq,
    SequenceBuild* seq_obj, ScoringFunction<int>* score_obj, ReachableReads* link_obj,
    double seed_score_scale, int assembly_depth, int band_size, double e_value,
    std::list<ContigType>* contigs
) {
  thread_id_ = thread_id;
  verbose_ = verbose;
  mutex_out_ = mutex_out;
  mutex_seq_ = mutex_seq;
  mutex_score_ = mutex_score;
  mutex_link_ = mutex_link;
  query_tag_ = query_tag;
  query_seq_ = query_seq;
  seq_obj_ = seq_obj;
  score_obj_ = score_obj;
  link_obj_ = link_obj;
  seed_score_scale_ = seed_score_scale;
  assembly_depth_ = assembly_depth;
  band_size_ = band_size;
  e_value_ = e_value;
  contigs_ = contigs;
  return;
}

AssembleFunctor::~AssembleFunctor(void) {return;}

void AssembleFunctor::Run(void) {
  double start_time = mytime();
  double check_time;
  
  if(verbose_)  {
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPx::Assemble info: [" << query_tag_ << "] scheduled." << std::endl;
    w_lock.unlock();
  }
  SLock seq_lock(*mutex_seq_);
  SLock score_lock(*mutex_score_);
  SLock link_lock(*mutex_link_);
  std::map<int, std::list<SeedType> > candidate_seeds;
  link_obj_->SelectSeeds(seed_score_scale_, *score_obj_, query_seq_, candidate_seeds);  
  //std::cout << "Select seeds done: " << candidate_seeds.size() << std::endl;
  std::map<int, std::list<ReadPairType> > seed_reads;
  link_obj_->GetSeedReads(*seq_obj_, candidate_seeds, seed_reads);  
  //std::cout << "Get seed reads done:  " << seed_reads.size() << std::endl;
  seq_lock.unlock();
  score_lock.unlock();
  link_lock.unlock();
  /// recruit neighbouring reads
  std::unordered_map<RIDType, bool> candidate_reads;
  //link_obj_->CollectNeighbourReads(collect_depth_, seed_reads, candidate_reads);
  //std::cout << "Collect neighbour reads done: " << candidate_reads.size() << std::endl;
  
  if(verbose_)  {
    check_time = mytime();
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPx::Assemble info: [" << query_tag_ << "] preprocessing done. ";
    printElapsed(start_time, check_time, "");
    w_lock.unlock();
    start_time = mytime();
  }
  
  std::list<ContigType> pre_contigs;
  AssembleExtend guided_assembly;
  guided_assembly.AssembleAllContigs(
      query_seq_, *seq_obj_, *link_obj_, *score_obj_, band_size_, 
      seed_reads, assembly_depth_, e_value_ * 100, pre_contigs
  );
  
  if(verbose_)  {
    check_time = mytime();
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPx::Assemble info: [" << query_tag_ << "] assembly done. ";
    printElapsed(start_time, check_time, "");
    w_lock.unlock();
    start_time = mytime();
  }
  //std::cout << "Assembly done:  " << pre_contigs.size() << std::endl;
  //return 0;
  // refine the contigs
  int overlap_len = link_obj_->GetOverlapLen();
  ContigRefinement recalibrate_contigs;
  //refined_contigs = contigs;
  recalibrate_contigs.RefineContigs(
      query_seq_, *seq_obj_, *score_obj_, band_size_, 
      e_value_, overlap_len, pre_contigs, *contigs_
  );
  
  if(verbose_)  {
    check_time = mytime();
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPx::Assemble info: [" << query_tag_ << "] recalibration done. ";
    printElapsed(start_time, check_time, "");
    std::cout <<  "GRASPx::Assemble info: [" << query_tag_ << "] finished." << std::endl;
    w_lock.unlock();
  }
  //std::cout << "Recalibration done: " << contigs_->size() << std::endl;
  return;
}


