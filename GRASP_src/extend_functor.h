#include "guided_assemble.h"
#include "extend_and_assemble.h"


// A Extension working functor, defined in the class form
class ExtendFunctor {
 public:
  ExtendFunctor(const int& in_phase_counter, void* in_saa_obj, const int in_query_seq_ID, 
      const int& in_query_kmer_begin, const std::string& in_target_kmer, 
      std::map<int, int>* in_clump_map,
      AssemblyGraph* in_assembly_graph_phase, std::list<VertexPairType>* in_seed_vertices,
      std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* in_shared_read_mapping, 
      LockType* in_mutex_hash,
      std::unordered_map<RIDType, ReadMapType>* in_shared_output,
      LockType* in_mutex_output,
      std::vector<bool> *in_space_indicator,
      LockType* in_mutex_indicator
  ) :  phase_counter_(in_phase_counter),
      saa_obj_(in_saa_obj), 
      query_seq_ID_(in_query_seq_ID),
      query_kmer_begin_(in_query_kmer_begin), 
      target_kmer_(in_target_kmer),
      clump_map_(in_clump_map),
      assembly_graph_phase_(in_assembly_graph_phase),
      seed_vertices_(in_seed_vertices), 
      shared_read_mapping_(in_shared_read_mapping),
      mutex_hash_(in_mutex_hash),
      shared_output_(in_shared_output),
      mutex_output_(in_mutex_output),
      space_indicator_(in_space_indicator),
      mutex_indicator_(in_mutex_indicator)
  {
    //std::cout << "quer_seq_ID: " << query_seq_ID_ << std::endl; 
    return;
  }
  ~ExtendFunctor()  {
    return; 
  }
  int FinalizeReadRecruitment(
      LockType* mutex_shared_pool,
      std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* shared_recruited_pool,
      std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* phase_recruited_pool
  );
  void FinalizeOutput(
      LockType* mutex_shared_output_pool,
      std::unordered_map<RIDType, ReadMapType>* shared_output_pool,
      std::unordered_map<RIDType, ReadMapType>* phase_output_pool
  );
  void run();
  
 private:
  int phase_counter_;
  void* saa_obj_; 
  const int& query_seq_ID_; 
  const int& query_kmer_begin_; 
  const std::string& target_kmer_; 
  std::map<int, int>* clump_map_;
  AssemblyGraph *assembly_graph_phase_; 
  std::list<VertexPairType>* seed_vertices_;
  std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* shared_read_mapping_; 
  LockType* mutex_hash_;
  std::unordered_map<RIDType, ReadMapType>* shared_output_;
  LockType* mutex_output_;
  std::vector<bool> *space_indicator_;
  LockType* mutex_indicator_;
};

void ExtendFunctor::run()  {
  //std::cout << "Extension called" << std::endl;
  GuidedAssemble* func_caller = (GuidedAssemble*) saa_obj_;
  //std::cout << query_seq_ID_ << std::endl;
  //std::cout << func_caller->query_sequence_[query_seq_ID_] << std::endl;
  std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> > *reads_to_record 
      = new std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >;
  ExtendAndAssemble guided_assemble_phase(
      func_caller->num_threads_, func_caller->query_sequence_[query_seq_ID_], 
      query_kmer_begin_, func_caller->seed_len_, target_kmer_, 
      *(func_caller->suffix_array_), *(func_caller->reverse_suffix_array_), 
      func_caller->n_back_check_, func_caller->dropoff_, *clump_map_,
      func_caller->score_scheme_, func_caller->right_band_, func_caller->down_band_, 
      0, reads_to_record, shared_read_mapping_, mutex_hash_, *assembly_graph_phase_, *seed_vertices_
  );
  //cerr << "object initialized" << endl;
  guided_assemble_phase.ExtendToBothDirections();
  //cerr << "key function finished" << endl;
  int num_new_recruitment = FinalizeReadRecruitment(mutex_hash_, shared_read_mapping_, reads_to_record); 
  if(num_new_recruitment <= 0)  {
    //delete assembly_graph_phase_;
    ULock w_lock(*mutex_indicator_);
    (*space_indicator_)[phase_counter_] = false;
    w_lock.unlock();
    delete assembly_graph_phase_;
    delete seed_vertices_;
    return;
  }
  (*space_indicator_)[phase_counter_] = true;
  delete reads_to_record;
  // get the high-score reads
  
  //std::unordered_map<RIDType, ReadMapType>* mapped_reads_phase 
  //    = new std::unordered_map<RIDType, ReadMapType>;
  //func_caller->TraceBack(
  //    func_caller->query_sequence_[query_seq_ID_], 
  //    *assembly_graph_phase_, *seed_vertices_, *mapped_reads_phase
  //);
  //FinalizeOutput(mutex_output_, shared_output_, mapped_reads_phase);
  // collecting the memory
  
  //delete mapped_reads_phase;
  //std::cout << "Extension end" << std::endl;
  return;
}

void ExtendFunctor::FinalizeOutput(
    LockType* mutex_shared_output_pool,
    std::unordered_map<RIDType, ReadMapType>* shared_output_pool,
    std::unordered_map<RIDType, ReadMapType>* phase_output_pool
) {
  ULock w_lock(*mutex_shared_output_pool);
  for(auto it = phase_output_pool->begin(); it != phase_output_pool->end(); ++ it) {
    it->second.phase_counter = phase_counter_;
    it->second.query_seq_ID = query_seq_ID_;
    auto it_shared = shared_output_pool->find(it->first);
    if(it_shared == shared_output_pool->end())  {
      shared_output_pool->insert({it->first, it->second});
    } else if(it->second.best_evalue < it_shared->second.best_evalue) {
      it_shared->second = it->second;
    }
    
  }
  w_lock.unlock();
  return;
}

int ExtendFunctor::FinalizeReadRecruitment(
    LockType* mutex_shared_pool,
    std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* shared_recruited_pool,
    std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* phase_recruited_pool
)  {
  int update_counter = 0;
  // setting the phase counter
  for(auto it = phase_recruited_pool->begin(); it != phase_recruited_pool->end(); ++ it) {
    for(auto it_map_phase = it->second.begin(); it_map_phase != it->second.end(); ++ it_map_phase) {
      it_map_phase->second.phase_counter_left = it_map_phase->second.phase_counter_right = phase_counter_;
    }
  }  
  ULock w_lock(*mutex_shared_pool);
  for(auto it = phase_recruited_pool->begin(); it != phase_recruited_pool->end(); ++ it) {
    auto it_assigned = shared_recruited_pool->find(it->first);
    if(it_assigned == shared_recruited_pool->end())  {
      // insert the entire read
      shared_recruited_pool->insert({it->first, it->second});
      ++ update_counter;
    } else  {
      for(auto it_map_phase = it->second.begin(); it_map_phase != it->second.end(); ++ it_map_phase) {
        auto it_map_shared = it_assigned->second.find(it_map_phase->first);
        if(it_map_shared == it_assigned->second.end())  {
          // if the entire clump is not mapped
          it_assigned->second[it_map_phase->first] = it_map_phase->second;
          ++ update_counter;
          continue;
        } 
        // if the read has presented, see if we need to update it
        if(it_map_phase->second.left_extended)  {
          // update the left extension information
          if((it_map_shared->second.left_extended && 
              it_map_phase->second.assembled_score_left > it_map_shared->second.assembled_score_left) ||
              !it_map_shared->second.left_extended)  {
            ++ update_counter;
            it_map_shared->second.left_extended = true;
            it_map_shared->second.phase_counter_left = it_map_phase->second.phase_counter_left;
            it_map_shared->second.assembled_score_left = it_map_phase->second.assembled_score_left;
            it_map_shared->second.phase_counter_left = it_map_phase->second.phase_counter_left;
            //it_map_shared->second.end_vertex_left = it_map_phase->second.end_vertex_left;
            //it_map_shared->second.end_edge_left = it_map_phase->second.end_edge_left;
          }
        } else  {
          // update the right extension information
          if((it_map_shared->second.right_extended && 
              it_map_phase->second.assembled_score_right > it_map_shared->second.assembled_score_right) ||
              !it_map_shared->second.right_extended)  {
            ++ update_counter;
            it_map_shared->second.right_extended = true;
            it_map_shared->second.phase_counter_right = it_map_phase->second.phase_counter_right;
            it_map_shared->second.assembled_score_right = it_map_phase->second.assembled_score_right;
            it_map_shared->second.phase_counter_right = it_map_phase->second.phase_counter_right;
            //it_map_shared->second.end_vertex_right = it_map_phase->second.end_vertex_right;
            //it_map_shared->second.end_edge_right = it_map_phase->second.end_edge_right;
          }
        }
      }
    }
  }
  w_lock.unlock();
  //std::cerr << "Update counter:  " << update_counter << std::endl;
  return update_counter;
}
