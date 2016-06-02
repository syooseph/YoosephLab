#include "hmm_profile.h"
#include "sequence_build.h"
#include "reduced_alphabet.h"
#include "database_index.h"
#include "reachable_reads.h"
#include "assemble_extend_p.h"
#include "contig_refinement_p.h"

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

class AssembleFunctorP {
 public:
  AssembleFunctorP(
      int thread_id, bool verbose, 
      LockType *mutex_out, LockType *mutex_seq, LockType *mutex_score, LockType *mutex_link,  
      HMMProfile &query_model, SequenceBuild* seq_obj, ReachableReads* link_obj,
      double seed_score_scale, int progressive_extend, bool is_dup_seeds, 
      int assembly_depth, int band_size, double p_value,
      std::string &output
  );
  ~AssembleFunctorP(void);
  void Run(void);
 private:
  int thread_id_;
  bool verbose_;
  LockType* mutex_out_;
  LockType* mutex_seq_;
  LockType* mutex_score_;
  LockType* mutex_link_;
  HMMProfile query_model_;
  SequenceBuild* seq_obj_;
  ReachableReads* link_obj_;
  double seed_score_scale_;
  int progressive_extend_;
  bool is_dup_seeds_;
  int assembly_depth_;
  int band_size_;
  double p_value_;
  std::string output_;
  
};

AssembleFunctorP::AssembleFunctorP(
    int thread_id, bool verbose, 
    LockType *mutex_out, LockType *mutex_seq, LockType *mutex_score, LockType *mutex_link,  
    HMMProfile &query_model, SequenceBuild* seq_obj, ReachableReads* link_obj,
    double seed_score_scale, int progressive_extend, bool is_dup_seeds, 
    int assembly_depth, int band_size, double p_value,
    std::string &output
) {
  thread_id_ = thread_id;
  verbose_ = verbose;
  mutex_out_ = mutex_out;
  mutex_seq_ = mutex_seq;
  mutex_score_ = mutex_score;
  mutex_link_ = mutex_link;
  query_model_ = query_model;
  seq_obj_ = seq_obj;
  link_obj_ = link_obj;
  seed_score_scale_ = seed_score_scale;
  progressive_extend_ = progressive_extend;
  is_dup_seeds_ = is_dup_seeds;
  assembly_depth_ = assembly_depth;
  band_size_ = band_size;
  p_value_ = p_value;
  output_ = output;
  return;
}

AssembleFunctorP ::~AssembleFunctorP(void) {return;}

void AssembleFunctorP::Run(void) {
  double start_time = mytime();
  double check_time;
  
  if(verbose_)  {
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPxp::Assemble info: [" << query_model_.GetName() << "] scheduled." << std::endl;
    w_lock.unlock();
  }
  SLock seq_lock(*mutex_seq_);
  SLock link_lock(*mutex_link_);
  
  std::map<int, std::list<SeedMatchType> > *seed_reads = new std::map<int, std::list<SeedMatchType> >;
  query_model_.GenHighScoreSeeds(-6.0, seed_score_scale_, *link_obj_, *seed_reads);
  //cout << "Finish getting high-score mers...  " << candidate_seeds.size() << endl;
  /*
  for(auto it = seed_reads->rbegin(); it != seed_reads->rend(); ++ it)  {
    std::cout << "SCORE:  " << it->first << std::endl;
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r)  {
      std::cout << "******" << std::endl;
      std::cout << it_r->seed << std::endl;
      std::cout << "******" << std::endl;
    }
  } 
  */ 
  
  seq_lock.unlock();
  link_lock.unlock();
  
  if(verbose_)  {
    check_time = mytime();
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPxp::Assemble info: [" << query_model_.GetName() << "] preprocessing done. ";
    printElapsed(start_time, check_time, "");
    w_lock.unlock();
    start_time = mytime();
  }
  
  std::list<ContigType> *pre_contigs = new std::list<ContigType>;
  std::list<ContigType> *post_contigs = new std::list<ContigType>;
  
  AssembleExtendP guided_assembly;
  guided_assembly.AssembleAllContigs(
      query_model_, *seq_obj_, *link_obj_, band_size_,
      *seed_reads, is_dup_seeds_, assembly_depth_,  
      p_value_, *pre_contigs
  );
  delete seed_reads;
  guided_assembly.ProgressiveExtension(*seq_obj_, *link_obj_, progressive_extend_, *pre_contigs);
  
  if(verbose_)  {
    check_time = mytime();
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPxp::Assemble info: [" << query_model_.GetName() << "] assembly done. ";
    printElapsed(start_time, check_time, "");
    w_lock.unlock();
    start_time = mytime();
  }
  
  

  //std::cout << "Assembly done:  " << pre_contigs.size() << std::endl;
  //return 0;
  // refine the contigs
  /*
  for(auto it_c = pre_contigs->begin(); it_c != pre_contigs->end(); ++ it_c) {
    post_contigs->push_back(*it_c);
  }
  */
  
  ContigRefinementP refine_worker;
  refine_worker.RefineContigsNoAln(
      *seq_obj_, link_obj_->GetOverlapLen(), 
      *pre_contigs, *post_contigs
  );
  delete pre_contigs;
  
  if(verbose_)  {
    check_time = mytime();
    ULock w_lock(*mutex_out_);
    std::cout <<  "GRASPxp::Assemble info: [" << query_model_.GetName() << "] recalibration done. ";
    printElapsed(start_time, check_time, "");
    std::cout <<  "GRASPxp::Assemble info: [" << query_model_.GetName() << "] finished." << std::endl;
    w_lock.unlock();
  }
  
  // aquire output MUTEX  
  ULock w_lock(*mutex_out_);
  std::ofstream out_fh;
  out_fh.open(output_.c_str(), std::ios_base::app);
  if(!out_fh.good())  {
    std::cout << output_ << std::endl;
    std::cout << "Error: graspxp-assemble: cannot write contig_out file." << std::endl;
    std::cout << "Please use \'--help\' for more details." << std::endl;
    exit(0);
  }
  int cid = 0;
  for(auto it_c = post_contigs->begin(); it_c != post_contigs->end(); ++ it_c) {
    out_fh << ">contig_" << cid << "||" << query_model_.GetName() << std::endl << it_c->sequence << std::endl;
    ++ cid;
  }
  out_fh.close();
  w_lock.unlock();
  
  
  delete post_contigs;
  //std::cout << "Recalibration done: " << contigs_->size() << std::endl;
  return;
}


