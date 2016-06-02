#include "guided_assemble.h"
#include "extend_and_assemble.h"
#include "assembly_graph.h"

typedef int AlignScoreType;

class CertificateFunctor {
 public:
  CertificateFunctor(
      void *in_saa_obj, const int& in_query_ID,
      AssemblyGraph* in_assembly_graph, std::list<VertexPairType>* in_source_vertices,
      std::set<RIDType>* in_reads_to_visit, std::list<AssembledPathType>* in_path_holder, LockType* in_mutex_holder,
      std::unordered_map<RIDType, ReadMapType>* in_shared_output,
      LockType* in_mutex_output
  ) :  saa_obj_(in_saa_obj),
      query_ID_(in_query_ID),
      assembly_graph_(in_assembly_graph),
      source_vertices_(in_source_vertices),
      reads_to_visit_(in_reads_to_visit),
      path_holder_(in_path_holder),
      mutex_holder_(in_mutex_holder),
      shared_output_(in_shared_output),
      mutex_output_(in_mutex_output)  
  {
    return; 
  } 
  ~CertificateFunctor(void) {return;}
  void run(void);
  
 private:
  void *saa_obj_;
  int query_ID_;
  AssemblyGraph* assembly_graph_;
  std::list<VertexPairType>* source_vertices_;
  std::set<RIDType>* reads_to_visit_;
  std::list<AssembledPathType>* path_holder_;
  LockType* mutex_holder_;
  std::unordered_map<RIDType, ReadMapType>* shared_output_;
  LockType* mutex_output_;
  std::unordered_map<RIDType, int> best_path_for_read_;
  void ReverseLeftPath(AssembledPathType& left_path);
  void UpdateReadLocations(int adjustment, AssembledPathType& path);
  int RecordReads(const AssembledPathType& path, const int& path_ID);
};

int CertificateFunctor::RecordReads(const AssembledPathType& path, const int& path_ID)  {
  int num_better_hits = 0;
  GuidedAssemble* func_caller = (GuidedAssemble*) saa_obj_;
  ULock w_lock(*mutex_output_);
  for(auto it = path.mapped_locations.begin(); it != path.mapped_locations.end(); ++ it) {
    double evalue = func_caller->score_scheme_->ComputeEValue(
      func_caller->query_sequence_[query_ID_].length(), func_caller->sample_size_ * 1000000, path.align_score
    );
    auto it_r = shared_output_->find(it->first);
    if(it_r == shared_output_->end())  {
      ReadMapType info;
      info.query_seq_ID = query_ID_;
      info.best_evalue = evalue;
      (*shared_output_)[it->first] = info;
      best_path_for_read_[it->first] = path_ID;
      ++ num_better_hits;
    } else  {
      if(it_r->second.best_evalue > evalue)  {
        it_r->second.best_evalue = evalue;
        it_r->second.query_seq_ID = query_ID_;
        best_path_for_read_[it->first] = path_ID;
        ++ num_better_hits;
      }
    }
  }
  w_lock.unlock();
  return num_better_hits;
}

void CertificateFunctor::ReverseLeftPath(AssembledPathType& left_path)  {
  GuidedAssemble* func_caller = (GuidedAssemble*) saa_obj_;  
  left_path.assembled_sequence = ReverseShortSequence(left_path.assembled_sequence);
  for(auto it_lr = left_path.mapped_locations.begin(); it_lr != left_path.mapped_locations.end(); ++ it_lr) {
    it_lr->second = left_path.assembled_sequence.length()
        - it_lr->second - strlen(func_caller->sample_seqs_[(unsigned int) it_lr->first]);
  }
  return;
}

void CertificateFunctor::UpdateReadLocations(int adjustment, AssembledPathType& path)  {
  if(adjustment == 0)  {
    return;
  }
  for(auto it = path.mapped_locations.begin(); it != path.mapped_locations.end(); ++ it) {
    it->second += adjustment;
  }
  return;
}

void CertificateFunctor::run()  {
  GuidedAssemble* func_caller = (GuidedAssemble*) saa_obj_;
  AlignScoreType raw_cutoff = func_caller->score_scheme_->ComputeRawScore(
    func_caller->query_sequence_[query_ID_].length(), func_caller->sample_size_ * 1000000, func_caller->evalue_cutoff_
  );
  // for each source vertex
  for(auto it = source_vertices_->begin(); it != source_vertices_->end(); ++ it) {
    // each path in the current direction is concaternated only with the best path in the opposite direction,
    // but not each one of them
    std::list<AssembledPathType> paths_left, paths_right;
    AssembledPathType best_path_left, best_path_right;
    if(it->has_initialized_left)  {
      assembly_graph_->SpellAllPaths(
          *func_caller->suffix_array_, *func_caller->reverse_suffix_array_, it->vertex_left, paths_left
      );
    }
    if(it->has_initialized_right)  {
      assembly_graph_->SpellAllPaths(
          *func_caller->suffix_array_, *func_caller->reverse_suffix_array_, it->vertex_right, paths_right
      );
    }
    
    // transform all left paths 
    std::list<AssembledPathType> paths_to_align;
    std::list<AssembledPathType> formatted_left_paths;
    if(paths_left.size() > 0)  {
      for(auto it_lp = paths_left.begin(); it_lp != paths_left.end(); ++ it_lp) {
        AssembledPathType reconstructed_path;
        ReverseLeftPath(*it_lp);
        reconstructed_path = *it_lp;
        // inferring the aligned position for the sequence
        reconstructed_path.mapped_begin = it->seed_start - it->seed_in_gap_seq 
           - (int) reconstructed_path.assembled_sequence.length() - (int) func_caller->down_band_;
        reconstructed_path.mapped_end = it->seed_start - it->seed_in_gap_seq + (int) func_caller->down_band_;
        reconstructed_path.align_score = it_lp->align_score;
        formatted_left_paths.push_back(reconstructed_path);
      }
    }
    // append the bridging region sequence and score
    if(formatted_left_paths.size() > 0)  {
      // append
      for(auto it_f = formatted_left_paths.begin(); it_f != formatted_left_paths.end(); ++ it_f) {
        for(auto it_gr = it->bridging_reads.begin(); it_gr != it->bridging_reads.end(); ++ it_gr) {
          it_f->mapped_locations.push_back(
              {it_gr->read_ID, it_f->assembled_sequence.length() + it_gr->offset_to_start}
          );
        }
        it_f->assembled_sequence += it->gap_sequence;
        it_f->mapped_end += it->gap_sequence.length();
        it_f->align_score += it->seed_match_score + it->left_score + it->right_score;
      }
    } else  {
      // directly put the gap sequence as the formatted left path
      AssembledPathType reconstructed_path;
      for(auto it_gr = it->bridging_reads.begin(); it_gr != it->bridging_reads.end(); ++ it_gr) {
        reconstructed_path.mapped_locations.push_back({it_gr->read_ID, it_gr->offset_to_start});
      }
      reconstructed_path.mapped_begin = it->seed_start - it->seed_in_gap_seq - (int) func_caller->down_band_;
      reconstructed_path.mapped_end = reconstructed_path.mapped_begin + it->gap_sequence.length() + (int) func_caller->down_band_;
      reconstructed_path.assembled_sequence = it->gap_sequence;
      reconstructed_path.align_score = it->seed_match_score + it->left_score + it->right_score;
      formatted_left_paths.push_back(reconstructed_path);
    }
    // append the right reads
    for(auto it_lp = formatted_left_paths.begin(); it_lp != formatted_left_paths.end(); ++ it_lp) {
      for(auto it_rp = paths_right.begin(); it_rp != paths_right.end(); ++ it_rp) {
        AssembledPathType reconstructed_path;
        reconstructed_path.mapped_begin = it_lp->mapped_begin;
        reconstructed_path.mapped_end = it_lp->mapped_end + it_rp->assembled_sequence.length();
        reconstructed_path.assembled_sequence = it_lp->assembled_sequence + it_rp->assembled_sequence;
        reconstructed_path.align_score = it_lp->align_score + it_rp->align_score;
        reconstructed_path.mapped_locations = it_lp->mapped_locations;
        for(auto it_rr = it_rp->mapped_locations.begin(); it_rr != it_rp->mapped_locations.end(); ++ it_rr) {
          reconstructed_path.mapped_locations.push_back(
              {it_rr->first, it_lp->assembled_sequence.length() + it_rr->second}
          );
        }
        if(reconstructed_path.align_score > raw_cutoff)  {
          
          // ******* delay the recruitment of individual reads ******
          if(RecordReads(reconstructed_path, paths_to_align.size())) {
            paths_to_align.push_back(reconstructed_path);
          }
          
        }
      }
    }
    
    // writing the paths to the shared path holder
    std::set<int> best_path_IDs;
    for(auto it_rm = best_path_for_read_.begin(); it_rm != best_path_for_read_.end(); ++ it_rm) {
      best_path_IDs.insert(it_rm->second);
    }
    ULock w_lock(*mutex_holder_);
    int p_id = 0;
    for(auto it_hp = paths_to_align.begin(); it_hp != paths_to_align.end(); ++ it_hp) {
      if(best_path_IDs.find(p_id) != best_path_IDs.end())  {
        path_holder_->push_back(*it_hp);
      }
      ++ p_id;
    }
    w_lock.unlock(); 
  }
  return;
}
