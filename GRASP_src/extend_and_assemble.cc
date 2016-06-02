#include "extend_and_assemble.h"

using namespace std;

ExtendAndAssemble::ExtendAndAssemble(void)  {
  return;
}

ExtendAndAssemble::~ExtendAndAssemble(void) {
  //delete source_vertices_;
  return;
}

ExtendAndAssemble::ExtendAndAssemble(
    const unsigned int& in_num_threads,
    const std::string& in_query_seq, 
    const unsigned int in_start,
    const unsigned int in_seed_length,
    const std::string& in_seed_target,
    GSA& in_suffix_array,
    GSA& in_reverse_suffix_array,
    const unsigned int in_n_back_check,
    const unsigned int in_dropoff,
    std::map<int, int>& in_clump_map,
    ScoringFunction<AlignScoreType>* in_score_scheme,
    const unsigned int in_right_band, const unsigned int in_down_band,
    const double in_bit_score_cutoff,
    std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* in_reads_to_record,
    std::unordered_map<RIDType, std::unordered_map<int, ReadTerminateType> >* in_assigned_reads,
    LockType* in_mutex_assigned_reads,
    AssemblyGraph& in_assembly_graph,
    std::list<VertexPairType>& in_source_vertices
)  {
  // TODO: add left and right tree references for storing the paths
  // copy in inputs
  num_threads_ = in_num_threads;
  query_seq_ = &in_query_seq;
  //cout << in_query_seq << endl;
  //cout << *query_seq_ << endl;
  start_ = in_start;
  seed_length_ = in_seed_length;
  seed_target_ = in_seed_target;
  suffix_array_ = &in_suffix_array;
  reverse_suffix_array_ = &in_reverse_suffix_array;
  n_back_check_ = in_n_back_check;
  dropoff_ = in_dropoff;
  clump_map_ = &in_clump_map;
  score_scheme_ = in_score_scheme;
  bit_score_cutoff_ = in_bit_score_cutoff;
  right_band_ = in_right_band;
  down_band_ = in_down_band;
  assigned_reads_ = in_assigned_reads;
  reads_to_record_ = in_reads_to_record;
  mutex_assigned_reads_ = in_mutex_assigned_reads;
  assembly_graph_ = &in_assembly_graph;
  source_vertices_ = &in_source_vertices;
  //source_vertices_ = new list<VertexPairType>;
  // setting other values
  alignment_counter_ = 0;
  num_concurrent_threads_ = 0;
  return;
}

void ExtendAndAssemble::FillExtendSeq(
    const std::string& out_path,
    AlignmentInfoType<AlignScoreType>& current_alignment
) {
  unsigned int extend_len = n_back_check_;
  unsigned int assembled_len = current_alignment.current_assembled_seq.length();
  assert(extend_len > 0);
  assert(assembled_len > 0);
  if(extend_len >= assembled_len)  {
    current_alignment.seed_extend_seq = current_alignment.current_assembled_seq;
    current_alignment.is_seed_complete = false;
  } else  {
    current_alignment.seed_extend_seq = 
        current_alignment.current_assembled_seq.substr(assembled_len - extend_len, extend_len);
    current_alignment.is_seed_complete = true;
  }
  return;
}

void ExtendAndAssemble::ExtendToBothDirections() {
  //cout << "ExtendToBothDirections Called" << endl;
  assert(high_score_alignments_right_.empty());
  assert(high_score_alignments_left_.empty());
  //cout << "Good here!!!" << endl;
  //cout << *query_seq_ << endl;
  //cout << start_ << endl;
  //cout << seed_length_ << endl;
  string query_seed = query_seq_->substr(start_, seed_length_);
  string rev_query_seed = ReverseShortSequence(query_seed);
  //cout << "Good here!!!" << endl;

  ////////////////// prepare right extension ////////////////////////
  AlignmentInfoType<AlignScoreType> seed_align_info_right;
  //seed_align_info_right.current_query_seq = query_seed;
  seed_align_info_right.current_query_length = query_seed.length();
  seed_align_info_right.current_assembled_seq = seed_target_;
  //seed_align_info_right.seed_extend_seq = it->substr(1, it->length() - 1);
  FillExtendSeq(seed_target_, seed_align_info_right);
  seed_align_info_right.seed_coverage = GetKmerCoverageRight(seed_target_);
  seed_align_info_right.alignment_ID = alignment_counter_ ++;
  SeqAlign<AlignScoreType> align_init_right(
      query_seed, seed_target_, score_scheme_, GLOBAL, 
      down_band_, right_band_, false
  );
  align_init_right.Align();
  seed_align_info_right.max_score = align_init_right.GetBestGlobalScore();
  seed_align_info_right.best_bit_score = score_scheme_->ComputeBitScore(seed_align_info_right.max_score);
  seed_align_info_right.prev_score = 0;
  ComputeAlignmentPriority(&seed_align_info_right); 
  //cout << "Good here!!!" << endl;
  ////////////////// prepare left extension ////////////////////////  
  string rev_loaded_seed = ReverseShortSequence(seed_target_);
  AlignmentInfoType<AlignScoreType> seed_align_info_left;
  //seed_align_info_left.current_query_seq = rev_query_seed;
  seed_align_info_left.current_query_length = rev_query_seed.length();
  seed_align_info_left.current_assembled_seq = rev_loaded_seed;
  //seed_align_info_left.seed_extend_seq = rev_loaded_seed.substr(1, rev_loaded_seed.length() - 1);
  FillExtendSeq(rev_loaded_seed, seed_align_info_left);
  seed_align_info_left.seed_coverage = GetKmerCoverageLeft(rev_loaded_seed);
  seed_align_info_left.alignment_ID = alignment_counter_ ++;
  SeqAlign<AlignScoreType> align_init_left(
      rev_query_seed, rev_loaded_seed, score_scheme_, GLOBAL, 
      down_band_, right_band_, false
  );
  align_init_left.Align();
  seed_align_info_left.max_score = align_init_left.GetBestGlobalScore();
  seed_align_info_left.best_bit_score = score_scheme_->ComputeBitScore(seed_align_info_left.max_score);
  seed_align_info_left.prev_score = 0;
  ComputeAlignmentPriority(&seed_align_info_left);
  // TODO: find an elegant way to prune the possible seeds
  if(seed_align_info_right.priority_score > 3 && seed_align_info_left.priority_score > 3)  { 
    // record the left and right alignment scores
    align_init_right.FillEdgeScore(seed_align_info_right.edge_score);
    align_init_left.FillEdgeScore(seed_align_info_left.edge_score);
    // initialize the non_inclusion_count
    seed_align_info_right.non_increasing_count = seed_align_info_right.non_inclusion_count = 0;
    seed_align_info_left.non_increasing_count = seed_align_info_left.non_inclusion_count = 0;
    // initialize the seed pair in the assembly graph  
    InitSeedPair(seed_align_info_left, seed_align_info_right);
    // insert the alignment to the extension queue
    InsertToAlignmentHolderRight(&seed_align_info_right);    
    InsertToAlignmentHolderLeft(&seed_align_info_left);
    // Do the extension left and right
    while(high_score_alignments_right_.size() > 0) {    
      ExtendToRightPhase();
    }
    while(high_score_alignments_left_.size() > 0) {
      ExtendToLeftPhase();
    } 
    
    //cout << "Num of vertices in the graph after extension:  " << assembly_graph_->GetNumVertices() << endl;
    
    // first try to connect the two extensions
    ProcessBridgingReads(seed_vertices_);
    // if cannot, take them as separate ones
    if(source_vertices_->size() <= 0)  { 
      ProcessIsolatedReads(seed_vertices_);
    }
  }
  
  return;
}

void ExtendAndAssemble::InitSeedPair(
    AlignmentInfoType<AlignScoreType>& left_alignment, 
    AlignmentInfoType<AlignScoreType>& right_alignment
)  {
    //cout << "Begin of insert seed pair" << endl;
  // hold the submission of the seed vertices to the graph, because most of the seeds do not lead to successful extensions      
    seed_vertices_.seed_match_score = left_alignment.max_score;
    seed_vertices_.seed_sequence = seed_target_;
    seed_vertices_.has_initialized_left = seed_vertices_.has_initialized_right = true;
   
    // initiating right vertex
    right_alignment.extend_direction = EXT_RIGHT;
    right_alignment.has_source_initialized = true;
    //right_alignment.seed_iterator = std::prev(seed_vertices_.end());
    VertexProperty sv_right;
    //sv_right.mutex_vertex = new LockType;
    sv_right.query_seed_begin = start_;
    sv_right.query_end = start_ + seed_length_ - 1;
    sv_right.touched = sv_right.visited = false;
    sv_right.max_score_below = sv_right.max_score_above = -_MAX_VALUE;
    right_alignment.source_vertex = assembly_graph_->AddVertex(sv_right);
    //right_alignment.seed_iterator->seed_vertex_right = right_alignment.source_vertex;
    seed_vertices_.seed_vertex_right = right_alignment.source_vertex;
    // initiating left vertex
    left_alignment.extend_direction = EXT_LEFT;
    left_alignment.has_source_initialized = true;
    //left_alignment.seed_iterator = std::prev(seed_vertices_.end());
    VertexProperty sv_left;
    //sv_left.mutex_vertex = new LockType;
    sv_left.query_seed_begin = start_ + seed_length_ - 1;
    sv_left.query_end = start_;
    sv_left.touched = sv_left.visited = false;
    sv_left.max_score_below = sv_left.max_score_above = -_MAX_VALUE;
    left_alignment.source_vertex = assembly_graph_->AddVertex(sv_left);
    seed_vertices_.seed_vertex_left = left_alignment.source_vertex;
    //left_alignment.seed_iterator->seed_vertex_left = left_alignment.source_vertex;
    //cout << "End of insert seed pair" << endl;    
    return;
}

void ExtendAndAssemble::ProcessIsolatedReads(const SeedVertexType& sv) {
  if(!sv.has_initialized_left && !sv.has_initialized_right)  {
    return;
  }
  VertexPairType left_partial;
  left_partial.vertex_left = sv.seed_vertex_left;
  left_partial.has_initialized_left = true;
  left_partial.has_initialized_right = false;
  left_partial.seed_match_score = sv.seed_match_score;
  left_partial.left_score = left_partial.right_score = 0;
  left_partial.gap_sequence = sv.seed_sequence;
  //left_partial.seed_sequence = sv.seed_sequence;
  left_partial.seed_in_gap_seq = 0; // the begin of the seed in the gap_sequence
  left_partial.seed_start = start_; // the begin of the seed in the query sequence
  source_vertices_->push_back(left_partial);
        
  VertexPairType right_partial;
  right_partial.vertex_right = sv.seed_vertex_right;
  right_partial.has_initialized_left = false;
  right_partial.has_initialized_right = true;
  right_partial.seed_match_score = sv.seed_match_score;
  right_partial.left_score = right_partial.right_score = 0;
  right_partial.gap_sequence = sv.seed_sequence;
  //right_partial.seed_sequence = sv.seed_sequence;
  right_partial.seed_in_gap_seq = 0;
  right_partial.seed_start = start_;
  source_vertices_->push_back(right_partial);
  
  return;
}

void ExtendAndAssemble::ProcessBridgingReads(const SeedVertexType& sv) {
  if(!sv.has_initialized_left && !sv.has_initialized_right)  {
    return;
  }
  list<BoostEdge> ie_list, oe_list;
  assembly_graph_->GetOutEdges(sv.seed_vertex_left, ie_list);
  assembly_graph_->GetOutEdges(sv.seed_vertex_right, oe_list);
  for(auto it_i = ie_list.begin(); it_i != ie_list.end(); ++ it_i) {
    BoostEdge ie_single = *it_i;
    BoostVertex iv_single = assembly_graph_->GetTarget(ie_single);
    EdgeProperty ie = assembly_graph_->AccessEdge(ie_single); 
    for(auto it_j = oe_list.begin(); it_j != oe_list.end(); ++ it_j) {
      BoostEdge oe_single = *it_j;
      BoostVertex ov_single = assembly_graph_->GetTarget(oe_single);
      EdgeProperty oe = assembly_graph_->AccessEdge(oe_single);
      list<BridgingReadType> valid_reads;  // the second field indicates the estimated mapped end
      string ie_outgoing_seq = string(reverse_suffix_array_->getSuffix_explicit(ie.outgoing_seq_ID, ie.outgoing_seq_Pos));
      string oe_outgoing_seq = string(suffix_array_->getSuffix_explicit(oe.outgoing_seq_ID, oe.outgoing_seq_Pos));
      string gap_sequence = ReverseShortSequence(ie_outgoing_seq) + sv.seed_sequence + oe_outgoing_seq;
      // check reads spanning the left vertex
      for(auto it_b = ie.candidate_bridging_reads.begin(); it_b != ie.candidate_bridging_reads.end(); ++ it_b) {
        string full_sequence(suffix_array_->getSequence_explicit((unsigned int) it_b->read_ID));
        unsigned int px_len = reverse_suffix_array_->getPos(static_cast<size_t>(it_b->rank_in_sfa));
        string cmp_seq = full_sequence.substr(full_sequence.length() - px_len, px_len); 
        if(oe_outgoing_seq.length() > px_len && oe_outgoing_seq.compare(0, px_len, cmp_seq) == 0)  {
          // TODO: check correctness
          it_b->offset_to_start = ie_outgoing_seq.length() + sv.seed_sequence.length() 
              - reverse_suffix_array_->getSuffixLength(static_cast<size_t>(it_b->rank_in_sfa));
          valid_reads.push_back(*it_b);
        }
      }
      // check reads spanning the right vertex
      for(auto it_b = oe.candidate_bridging_reads.begin(); it_b != oe.candidate_bridging_reads.end(); ++ it_b) {
        string full_sequence(reverse_suffix_array_->getSequence_explicit((unsigned int) it_b->read_ID));
        unsigned int px_len = suffix_array_->getPos(static_cast<size_t>(it_b->rank_in_sfa));
        string cmp_seq = full_sequence.substr(full_sequence.length() - px_len, px_len);        
        if(ie_outgoing_seq.length() > px_len && ie_outgoing_seq.compare(0, px_len, cmp_seq) == 0)  {
          // TODO: check correctness
          it_b->offset_to_start = ie_outgoing_seq.length() - it_b->position;
          valid_reads.push_back(*it_b);
        }
      }
      // TODO: adding briding reads to the hash table!!!
      // add and edge between the source and the target node
      
      if(valid_reads.size() > 3)  {
        // record the vertex pair 
        VertexPairType bridge_link;
        bridge_link.vertex_left = iv_single;
        bridge_link.vertex_right = ov_single;
        bridge_link.seed_start = start_;  // the beginning of the seed sequence in the query sequence 
        bridge_link.gap_sequence = gap_sequence; 
        //bridge_link.seed_sequence = sv.seed_sequence;
        bridge_link.seed_in_gap_seq = ie_outgoing_seq.length(); // the beginning of the seed sequence in gap_sequence
        bridge_link.bridging_reads = valid_reads;
        bridge_link.seed_match_score = sv.seed_match_score;
        bridge_link.left_score = ie.increased_score;
        bridge_link.right_score = oe.increased_score;
        bridge_link.has_initialized_left = bridge_link.has_initialized_right = true;
        source_vertices_->push_back(bridge_link);  
        // record the read mappings
        for(auto it_list = valid_reads.begin(); it_list != valid_reads.end(); ++ it_list) {
          // record the read
          // record the mapping information as an assigned reads
          ReadTerminateType recorder;
          recorder.right_extended = recorder.left_extended = false;
          recorder.assembled_score_left = recorder.assembled_score_right = -_MAX_VALUE;
          recorder.alignment_position = it_list->alignment_position;
          if(it_list->extend_direction == EXT_RIGHT)  {
            recorder.right_extended = true;
            recorder.assembled_score_right = bridge_link.seed_match_score + bridge_link.left_score + bridge_link.right_score;
            //recorder.rank_in_sfa_right = it_list->rank_in_sfa;
            //recorder.end_vertex_right = ov_single;
            //recorder.end_edge_right = oe_single;
            //recorder.assembled_score = oe.current_score;
          } else  {
            recorder.left_extended = true;
            recorder.assembled_score_left = bridge_link.seed_match_score + bridge_link.left_score + bridge_link.right_score;
            //recorder.rank_in_sfa_left = it_list->rank_in_sfa;
            //recorder.end_vertex_left = iv_single;
            //recorder.end_edge_left = ie_single;
            //recorder.assembled_score = ie.current_score;
          }
          UpdateAssignedReads(it_list->read_ID, recorder);
        }
      }
    }
  }
  return;
}


int ExtendAndAssemble::ComputeBridgingReadsRightBound(const SfaType& rank, const DirectionType& direction) {
  unsigned int sx_len;
  if(direction == EXT_RIGHT)  {
    sx_len = (unsigned int) suffix_array_->getSuffixLength(static_cast<size_t>(rank));
  } else  {
    sx_len = (unsigned int) reverse_suffix_array_->getPos(static_cast<size_t>(rank)) + seed_length_;
  }
  int q_new_len = sx_len + down_band_ - 1;
  int end = start_ + q_new_len - 1;
  return end;
}

void ExtendAndAssemble::DoExtendAlignment(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    const std::string& query_seq_phase, const std::string& out_seq_phase,
    const unsigned int& read_ID, const unsigned int& position,
    const std::list<ReadAssignmentType>& included_reads,
    const std::list<BridgingReadType>& candidate_reads,
    AlignmentInfoType<AlignScoreType>& current_alignment
) {
  // this function does the actual alignment and record the alignment output in "current_alignment"
  string current_query_seq;
  if(source_alignment.extend_direction == EXT_RIGHT)  {
    current_query_seq = query_seq_->substr(start_, source_alignment.current_query_length);
  } else  {
    current_query_seq = query_seq_->substr(
        start_ + seed_length_ - source_alignment.current_query_length, 
        source_alignment.current_query_length
    );
    current_query_seq = ReverseShortSequence(current_query_seq);
  }
  //cout << "sequence used for alignment: " << endl;
  //cout << current_query_seq <<  " " << source_alignment.current_assembled_seq << endl;
  //cout << query_seq_phase << "  " << out_seq_phase << endl;
  //cout << "input score: " << source_alignment.edge_score.row[0] << "  " << source_alignment.edge_score.column[0] << endl;
  //cout << "input score: " << source_alignment.edge_score.row_ins[0] << "  " << source_alignment.edge_score.column_ins[0] << endl;
  //cout << "input score: " << source_alignment.edge_score.row_del[0] << "  " << source_alignment.edge_score.column_del[0] << endl;
  //cout << "end" << endl;

  SeqAlignExtend<AlignScoreType> align_phase = SeqAlignExtend<AlignScoreType>(
      current_query_seq, source_alignment.current_assembled_seq,
      query_seq_phase, out_seq_phase, source_alignment.edge_score, score_scheme_,
      GLOBAL, right_band_, down_band_, false
  );

  align_phase.Align();
  AlignScoreType align_score = align_phase.GetBestEdgeScore();  
  align_phase.FillEdgeScore(current_alignment.edge_score);

  //delete align_phase;
  //cout << "alignment score: " << align_score << endl; 
  // record this alignment extention and add it to high_score_alignments_
  current_alignment.alignment_ID = alignment_counter_ ++;
  // update the accumulated sequence
  //current_alignment.current_query_seq = source_alignment.current_query_seq + query_seq_phase;
  current_alignment.current_query_length = source_alignment.current_query_length + query_seq_phase.length();
  current_alignment.current_assembled_seq = source_alignment.current_assembled_seq + out_seq_phase;
  // record extension direction
  current_alignment.extend_direction = source_alignment.extend_direction;
  // update extension sequence by referring the the length set by n_back_check   
  FillExtendSeq(out_seq_phase, current_alignment);
  // record the alignment score
  current_alignment.max_score = align_score;
  
  current_alignment.best_bit_score = source_alignment.best_bit_score;
  double bit_score = score_scheme_->ComputeBitScore(align_score);
  if(bit_score > current_alignment.best_bit_score)  {
    current_alignment.best_bit_score = bit_score;
  }
  current_alignment.prev_score = source_alignment.max_score;
  if(current_alignment.max_score <= current_alignment.prev_score)  {
    current_alignment.non_increasing_count = source_alignment.non_increasing_count + 1;
  } else  {
    current_alignment.non_increasing_count = 0;
  }
  
  // count num reads that support the extension and record the recruited reads
  // and also based on the coverage compute the alignment priority
  current_alignment.seed_coverage = included_reads.size();
  ComputeAlignmentPriority(&current_alignment);
  // record the seed iterator
  current_alignment.has_source_initialized = source_alignment.has_source_initialized;
  //current_alignment.seed_iterator = source_alignment.seed_iterator;
  // record the current alignment information as vertex and edge information
  current_alignment.source_vertex = source_alignment.source_vertex;
  // the edge
  EdgeProperty alignment_edge;
  // TODO: write a function similar to RecordRecruitedReadsRight to fill included reads into the edge
  //alignment_edge.mutex_edge = new LockType;
  alignment_edge.included_reads = included_reads;
  //alignment_edge.outgoing_seq = out_seq_phase;
  alignment_edge.outgoing_seq_ID = read_ID;
  alignment_edge.outgoing_seq_Pos = position;
  if(source_alignment.extend_direction == EXT_RIGHT)  {
    //int index = start_ + current_alignment.current_query_seq.length() - 1;
    int index = start_ + current_alignment.current_query_length - 1;
    if(index > (int) query_seq_->length() - 1)  {
      index = query_seq_->length() - 1;
    }
    alignment_edge.query_end = index;
  } else  {
    //int index = start_ - current_alignment.current_query_seq.length() + seed_length_;
    int index = start_ - current_alignment.current_query_length + seed_length_;
    if(index < 0)  {
      index = 0;
    }
    alignment_edge.query_end = index;
  }
  alignment_edge.current_score = align_score;
  alignment_edge.increased_score = current_alignment.max_score - current_alignment.prev_score;
  alignment_edge.extend_direction = source_alignment.extend_direction;
  alignment_edge.candidate_bridging_reads = candidate_reads;
  // the vertex
  VertexProperty alignment_vertex;
  //alignment_vertex.mutex_vertex = new LockType;
  alignment_vertex.query_seed_begin = start_;
  alignment_vertex.query_end = alignment_edge.query_end;
  alignment_vertex.touched = alignment_vertex.visited = false;
  alignment_vertex.max_score_below = alignment_vertex.max_score_above = -_MAX_VALUE;
  // record them into alignment
  current_alignment.edge_holder = source_alignment.edge_holder;
  current_alignment.vertex_holder = source_alignment.vertex_holder;
  current_alignment.non_increasing_count = source_alignment.non_increasing_count;
  current_alignment.non_inclusion_count = source_alignment.non_inclusion_count;
#if DEBUG
  cout << "DoExtendAlignment: size vertex before adding:  " << current_alignment.vertex_holder.size() << endl;
#endif
  current_alignment.edge_holder.push_back(alignment_edge);
  current_alignment.vertex_holder.push_back(alignment_vertex); 
#if DEBUG
  cout << "DoExtendAlignment: size vertex after adding:  " << current_alignment.vertex_holder.size() << endl;
#endif
  return;
}

bool ExtendAndAssemble::LookupExtendedPaths(
    const OutgoingPathType& extend_path, 
    int& num_new_reads, BoostEdge& end_edge, BoostVertex& end_vertex
) {
  // checks all reads included in the "extend_path"

  //unordered_map<BoostVertex, int> vertex_occurrence;
  //BoostVertex max_vertex;
  //BoostEdge max_edge;
  //int max_count = 0;
  int num_redundant_reads = 0;
  num_new_reads = 0;
  for(auto it = extend_path.included_reads.begin(); it != extend_path.included_reads.end(); ++ it) {
    //BoostVertex v;
    //BoostEdge e;
    if(LookupAssignedVertex(*it))  {
      ++ num_redundant_reads;
      /*
      if(vertex_occurrence.find(v) != vertex_occurrence.end())  {
        ++ vertex_occurrence[v];
      } else  {
        vertex_occurrence[v] = 1;
      }
      // record the end_vertex that occur most frequently
      if(vertex_occurrence[v] > max_count)  {
        max_count = vertex_occurrence[v];
        //max_vertex = v;
        //max_edge = e;
      }
      */
    } else  {
      num_new_reads ++;
    }
  }
  for(auto it = extend_path.candidate_bridging_reads.begin(); it != extend_path.candidate_bridging_reads.end(); ++ it) {
    //BoostVertex v;
    //BoostEdge e;
    if(LookupAssignedVertex(*it))  {
      ++ num_redundant_reads;
      /*
      if(vertex_occurrence.find(v) != vertex_occurrence.end())  {
        ++ vertex_occurrence[v];
      } else  {
        vertex_occurrence[v] = 1;
      }
      // record the end_vertex that occur most frequently
      if(vertex_occurrence[v] > max_count)  {
        max_count = vertex_occurrence[v];
        //max_vertex = v;
        //max_edge = e;
      }
      */
    } else  {
      num_new_reads ++;
    }
  }
  if(num_redundant_reads == (int) (extend_path.included_reads.size() + extend_path.candidate_bridging_reads.size()))  { 
    // a more stringent criteiron to solve the non-deterministic problem for multi-threaded
    //cerr << "Redundant detected" << endl;
    //end_vertex = max_vertex;
    //end_edge = max_edge;
    return true;
  } else  {
    return false;
  }
}

bool ExtendAndAssemble::DetectPathMerge(
    AlignmentInfoType<AlignScoreType>& source_alignment,
    OutgoingPathType& out_path
) {
  BoostEdge end_edge;
  BoostVertex end_vertex;
  EdgeProperty e;
  int num_new_reads;
  if(LookupExtendedPaths(out_path, num_new_reads, end_edge, end_vertex))  {   
    //cerr << "Redundant extension found" << endl; 
    return true;
    /*
    if(ConsolidateMergedPaths(source_alignment, out_path, end_edge, end_vertex, e)) {
      if(source_alignment.has_source_initialized)  { 
        assembly_graph_->AddEdge(source_alignment.source_vertex, end_vertex, e);
      }
    }
    return true;
    */
  }
  return false;
}

unsigned int ExtendAndAssemble::GetNumConcurrentThreads(void) {
  //SLock r_lock(mutex_thread_counter_);
  //int n = num_concurrent_threads_;
  //r_lock.unlock();
  return 0;
}

void ExtendAndAssemble::IncreaseThreadCounter(void) {
  //ULock w_lock(mutex_thread_counter_);
  //++ num_concurrent_threads_;
  //w_lock.unlock();
  return;
}

void ExtendAndAssemble::DecreaseThreadCounter(void) {
  //ULock w_lock(mutex_thread_counter_);
  //-- num_concurrent_threads_;
  //w_lock.unlock();
  return;
}

void ExtensionWorkerWrapper(bool is_right, void* ext_obj, void* alignment_info, void* path_info, BoostEdge& ext_e, BoostVertex& ext_v) {
  ExtendAndAssemble* func_caller = (ExtendAndAssemble*) ext_obj;
  //cout << "ExtensionWorkerWrapper::all " << ((AlignmentInfoType<AlignScoreType>*) alignment_info)->edge_score.row[0] << endl;
  // increase the thread counter
  func_caller->IncreaseThreadCounter();
  if(is_right)  {
    func_caller->RightExtensionWorker(alignment_info, path_info, ext_e, ext_v);
  } else  {
    func_caller->LeftExtensionWorker(alignment_info, path_info, ext_e, ext_v);
  }
  // decrease the thread counter
  func_caller->DecreaseThreadCounter();
  return;
}

void ExtendAndAssemble::LeftExtensionWorker(void* alignment_info, void* path_info, BoostEdge& ext_e, BoostVertex& ext_v) {

  //cout << "Left extend worker begin" << endl;
  AlignmentInfoType<AlignScoreType>* alignment_to_extend = (AlignmentInfoType<AlignScoreType>*) alignment_info;
  OutgoingPathType* out_path = (OutgoingPathType*) path_info;

  unsigned int r_id = reverse_suffix_array_->getId(static_cast<size_t>(out_path->key_sequence_position));
  unsigned int r_pos = 
      reverse_suffix_array_->getPos(static_cast<size_t>(out_path->key_sequence_position)) + 
      out_path->len_search_seed;
  string out_sequence = string(reverse_suffix_array_->getSuffix_explicit(r_id, r_pos));
  if(out_sequence.length() <= 0)  {
    return;
  }
    
  string query_extend_seq = GetQueryExtendSeqLeft(*alignment_to_extend, out_path->key_sequence_position);
  AlignmentInfoType<AlignScoreType> current_alignment;
  DoExtendAlignment(
      *alignment_to_extend, query_extend_seq, out_sequence, r_id, r_pos,
      out_path->included_reads, out_path->candidate_bridging_reads, current_alignment
  );
  double align_bit_score = score_scheme_->ComputeBitScore(current_alignment.max_score);
  int num_new_reads = SubmitEdgeReads(current_alignment, ext_e, ext_v);
  if(num_new_reads > 0)  {
    current_alignment.non_inclusion_count = 0;
  } else  {
    ++ current_alignment.non_inclusion_count;
  }
  
  if(align_bit_score > bit_score_cutoff_ && align_bit_score + dropoff_ >= current_alignment.best_bit_score && 
      current_alignment.non_increasing_count < 3 && current_alignment.non_inclusion_count < 3
  )  {  
    // if the current alignment is good enough, insert it to the extension queue
    ComputeAlignmentPriority(&current_alignment); 
    InsertToAlignmentHolderLeft(&current_alignment);
  } 
  return;
}

void ExtendAndAssemble::ExtendToLeftPhase(void) {

  // takes the first alignment path which has the highest prioirty
  //ULock w_lock(mutex_holder_queue_left_);
  AlignmentInfoType<AlignScoreType> alignment_to_extend = 
      *high_score_alignments_left_.begin();
  high_score_alignments_left_.pop_front();
  //w_lock.unlock();
  
  if(alignment_to_extend.edge_score.begin_column >= alignment_to_extend.current_query_length)  {
    return;
  }
  IndexSample decode_caller;
  // computes the possible outgoing paths
  list<OutgoingPathType> outgoing_paths;
  ComputeOutgoingPathsLeft(alignment_to_extend, outgoing_paths);
  // holding paths that share no redundant read recruitment
  list<OutgoingPathType> paths_to_extend;
  for(auto it = outgoing_paths.begin(); it != outgoing_paths.end(); ++ it) {
    //LeftExtensionWorker(&alignment_to_extend, &(*it));
    if(!DetectPathMerge(alignment_to_extend, *it))  {
      paths_to_extend.push_back(*it);
    }
  }
  boost::thread_group extension_group;
  list<AlignmentInfoType<AlignScoreType>* > alignment_pointer_holder;
  //int num_concurrent_threads = 0;
  // for each outgoing path, do the alignment and insert the alignment to list
  //cout << "num paths to extend: " << paths_to_extend.size() << endl;
  auto it = paths_to_extend.begin();
  while(it != paths_to_extend.end()) {
    AlignmentInfoType<AlignScoreType>* alignment_copy = new AlignmentInfoType<AlignScoreType>;
    *alignment_copy = alignment_to_extend;
    alignment_pointer_holder.push_back(alignment_copy);
    // add place-holder vertex
    VertexProperty foo_vertex;
    //foo_vertex.mutex_vertex = new LockType;
    BoostVertex extend_vertex = assembly_graph_->AddVertex(foo_vertex);
    // add place-holder edge
    EdgeProperty foo_edge;
    //foo_edge.mutex_edge = new LockType;
    BoostEdge extend_edge = assembly_graph_->AddEdge(alignment_to_extend.source_vertex, extend_vertex, foo_edge);
    // try to obtain execution slot
    /*
    do  {
      num_concurrent_threads = GetNumConcurrentThreads();
    } while(num_concurrent_threads >= num_threads_);
    boost::thread *ext_phase = new boost::thread(
        ExtensionWorkerWrapper, 0, this, (void *) alignment_copy, (void *) &(*it), extend_edge, extend_vertex
    );
    extension_group.add_thread(ext_phase);
    */
    
    ExtensionWorkerWrapper(0, this, (void *) alignment_copy, (void *) &(*it), extend_edge, extend_vertex);
    ++ it;
  }  
  extension_group.join_all();
  if(!alignment_pointer_holder.empty())  {
    auto it_holder = alignment_pointer_holder.begin();
    while(it_holder != alignment_pointer_holder.end()) {
      delete *it_holder;
      ++ it_holder;
    }
  }
  return;
}

void ExtendAndAssemble::RightExtensionWorker(void* alignment_info, void* path_info, BoostEdge& ext_e, BoostVertex& ext_v)  {

  //cout << "Right extend worker begin" << endl;
  AlignmentInfoType<AlignScoreType>* alignment_to_extend = (AlignmentInfoType<AlignScoreType>*) alignment_info;
  OutgoingPathType* out_path = (OutgoingPathType*) path_info;
  
  unsigned int r_id = suffix_array_->getId(static_cast<size_t>(out_path->key_sequence_position));
  unsigned int r_pos = 
      suffix_array_->getPos(static_cast<size_t>(out_path->key_sequence_position)) + 
      out_path->len_search_seed;
  string out_sequence = string(suffix_array_->getSuffix_explicit(r_id, r_pos));
   
  if(out_sequence.length() <= 0)  {
    return;
  }
   
  string query_extend_seq = GetQueryExtendSeqRight(*alignment_to_extend, out_path->key_sequence_position);

  // construct and do the alignment
  AlignmentInfoType<AlignScoreType> current_alignment;
  //cout << query_extend_seq << " " << out_sequence << endl;
  DoExtendAlignment(
      *alignment_to_extend, query_extend_seq, out_sequence, r_id, r_pos,
      out_path->included_reads, out_path->candidate_bridging_reads, current_alignment
  );    
  double align_bit_score = score_scheme_->ComputeBitScore(current_alignment.max_score);
  int num_new_reads = SubmitEdgeReads(current_alignment, ext_e, ext_v);
  if(num_new_reads > 0)  {
    current_alignment.non_inclusion_count = 0;
  } else  {
    ++ current_alignment.non_inclusion_count;
  }

  if(align_bit_score > bit_score_cutoff_ && align_bit_score + dropoff_ >= current_alignment.best_bit_score &&
      current_alignment.non_increasing_count < 3 && current_alignment.non_inclusion_count < 3
  )  {
    ComputeAlignmentPriority(&current_alignment); 
    InsertToAlignmentHolderRight(&current_alignment);
  }
  return; 
}

void ExtendAndAssemble::ExtendToRightPhase(void)  {

  //cout << "ExtendToRight Begin" << endl;
  //ULock w_lock(mutex_holder_queue_right_);
  AlignmentInfoType<AlignScoreType> alignment_to_extend = 
      *high_score_alignments_right_.begin();
  high_score_alignments_right_.pop_front();
  //w_lock.unlock();

  if(alignment_to_extend.edge_score.begin_column >= alignment_to_extend.current_query_length)  {
    return;
  }
  IndexSample decode_caller;
  // computes the possible outgoing paths
  list<OutgoingPathType> outgoing_paths;
  ComputeOutgoingPathsRight(alignment_to_extend, outgoing_paths);
  list<OutgoingPathType> paths_to_extend;
  for(auto it = outgoing_paths.begin(); it != outgoing_paths.end(); ++ it) {
    //LeftExtensionWorker(&alignment_to_extend, &(*it));
    if(!DetectPathMerge(alignment_to_extend, *it))  {
      paths_to_extend.push_back(*it);
    }
  }
  //cout << "num paths to extend: " << paths_to_extend.size() << endl;
  // for each outgoing path, do the alignment and insert the alignment to list
  boost::thread_group extension_group;
  list<AlignmentInfoType<AlignScoreType>* > alignment_pointer_holder;
  //int num_concurrent_threads = 0;
  auto it = paths_to_extend.begin();
  while(it != paths_to_extend.end()) {
    AlignmentInfoType<AlignScoreType>* alignment_copy = new AlignmentInfoType<AlignScoreType>;
    *alignment_copy = alignment_to_extend;
    alignment_pointer_holder.push_back(alignment_copy);
    // add place-holder vertex
    VertexProperty foo_vertex;
    //foo_vertex.mutex_vertex = new LockType;
    BoostVertex extend_vertex = assembly_graph_->AddVertex(foo_vertex);
    // add place-holder edge
    EdgeProperty foo_edge;
    //foo_edge.mutex_edge = new LockType;
    BoostEdge extend_edge = assembly_graph_->AddEdge(alignment_to_extend.source_vertex, extend_vertex, foo_edge);
    // try to obtain execution slot
    /*
    do  {
      num_concurrent_threads = GetNumConcurrentThreads();
    } while(num_concurrent_threads >= num_threads_);
    boost::thread *ext_phase = new boost::thread(
        ExtensionWorkerWrapper, 1, this, (void *) alignment_copy, (void *) &(*it), extend_edge, extend_vertex
    );
    extension_group.add_thread(ext_phase); 
    */
    ExtensionWorkerWrapper(1, this, (void *) alignment_copy, (void *) &(*it), extend_edge, extend_vertex);
    ++ it;
  }
  extension_group.join_all();
  if(!alignment_pointer_holder.empty())  {
    auto it_holder = alignment_pointer_holder.begin();
    while(it_holder != alignment_pointer_holder.end()) {
      delete *it_holder;
      ++ it_holder;
    }
  }
  //cout << "###########################################ExtensionRightPhase end: " << endl;
  return;
}

int ExtendAndAssemble::ComputeOutSeqDifference(const ReadAssignmentType& extend_read, const ReadAssignmentType& edge_read) {
  int extend_pos, edge_pos;
  if(extend_read.extend_direction == EXT_RIGHT && edge_read.extend_direction == EXT_RIGHT)  {
    extend_pos = suffix_array_->getPos(static_cast<size_t>(extend_read.rank_in_sfa));
    edge_pos = suffix_array_->getPos(static_cast<size_t>(edge_read.rank_in_sfa));
    return edge_pos - extend_pos;
  } else if(extend_read.extend_direction == EXT_LEFT && edge_read.extend_direction == EXT_LEFT) {
    extend_pos = reverse_suffix_array_->getPos(static_cast<size_t>(extend_read.rank_in_sfa));
    edge_pos = reverse_suffix_array_->getPos(static_cast<size_t>(edge_read.rank_in_sfa));
    return edge_pos - extend_pos;
  } 
  else  {
    return _MAX_VALUE;
  }
}

/*
bool ExtendAndAssemble::ConsolidateMergedPaths(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    OutgoingPathType& extend_path, const BoostEdge& end_edge, 
    const BoostVertex& end_vertex, EdgeProperty& consolidated_edge
)  {
  EdgeProperty end_edge_content = assembly_graph_->AccessEdge(end_edge);
#if DEBUG
  cout << "ConsolidateMergedPaths: begin of function" << endl;
  cout << "ConsolidateMergedPaths: edge sequence: " << end_edge_content.outgoing_seq << endl;
#endif
  // get reads (in both included_reads and candidate_bridging_reads) that spells the edge
  map<RIDType, ReadAssignmentType> reads_in_edge;
  list<ReadAssignmentType> all_reads_edge = end_edge_content.included_reads;
  all_reads_edge.splice(all_reads_edge.end(), end_edge_content.candidate_bridging_reads);
  for(auto it = all_reads_edge.begin(); it != all_reads_edge.end(); ++ it) {
#if DEBUG
    cout << "ConsolidateMergedPaths: included reads" << endl;
    if(it->extend_direction == EXT_RIGHT)  {
      cout << "ConsolidateMergedPaths: included reads:  " << (unsigned int) it->read_ID << " extends right" << endl;
    } else  {
      cout << "ConsolidateMergedPaths: included reads:  " << (unsigned int) it->read_ID << " extends left" << endl;
    }
#endif
    reads_in_edge[it->read_ID] = *it;
  }
  bool has_initialized = false;
  DirectionType extend_direction = EXT_RIGHT, edge_direction = EXT_RIGHT;
  int distance = -_MAX_VALUE;
  list<ReadAssignmentType> all_reads_extend = extend_path.included_reads;
  all_reads_extend.splice(all_reads_extend.end(), extend_path.candidate_bridging_reads);
  for(auto it = all_reads_extend.begin(); it != all_reads_extend.end(); ++ it) {
    auto edge_it = reads_in_edge.find(it->read_ID);
    if(edge_it == reads_in_edge.end())  {
      continue;
    }
    if(!has_initialized)  {
      has_initialized = true;
      extend_direction = it->extend_direction;
      edge_direction = edge_it->second.extend_direction;
      distance = ComputeOutSeqDifference(*it, edge_it->second);
    }
    if(extend_direction != it->extend_direction || 
      edge_direction != edge_it->second.extend_direction || 
      extend_direction != edge_direction ||
      distance != ComputeOutSeqDifference(*it, edge_it->second)
    )  {
#if DEBUG
      cout << "ConsolidateMergedPaths:  failed" << endl;
#endif
      return false;
    }
  }
  if(has_initialized)  {
    BuildNewEdge(source_alignment, extend_path, end_edge_content, extend_direction, edge_direction, distance, consolidated_edge);
  }
  if(consolidated_edge.included_reads.size() > 0)  {
#if DEBUG
    cout << "ConsolidateMergedPaths:  success" << endl;
#endif
    return true;
  } else  {
#if DEBUG
    cout << "ConsolidateMergedPaths:  fail" << endl;
#endif
    return false;
  }
}
*/

bool ExtendAndAssemble::IsSeqEndMatch(std::string seqA, std::string seqB) {
  if(seqA == "" || seqB == "")  {
    return false;
  }
  if(seqA.length() >= seqB.length())  {
    if(seqA.compare(seqA.length() - seqB.length(), seqB.length(), seqB) == 0)  {
      return true;
    } else  {
      return false;
    }
  } else  {
    if(seqB.compare(seqB.length() - seqA.length(), seqA.length(), seqA) == 0)  {
      return true;
    } else  {
      return false;
    }
  }
}

/*
void ExtendAndAssemble::BuildNewEdge(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    const OutgoingPathType& extend_path, const EdgeProperty& edge_content, 
    const DirectionType& extend_direction, const DirectionType& read_direction, int distance,
    EdgeProperty& new_edge
) {
 
#if DEBUG
  cout << "BuildNewEdge:  begin of function" << endl;
#endif
  assert(extend_direction == read_direction);
  //new_edge.mutex_edge = new LockType;
  new_edge.query_end = edge_content.query_end;
  new_edge.extend_direction = extend_direction;
  // get the query sequence to be aligned
  int query_end;
  if(source_alignment.edge_holder.size() > 0)  { // when there are alignments that are not submitted, get the last position
    query_end = source_alignment.edge_holder.rbegin()->query_end;
  } else if(source_alignment.has_source_initialized) {
    query_end = (assembly_graph_->AccessVertex(source_alignment.source_vertex)).query_end;
  } else  {
    return;
  }
  string query_aln_seq = "";
  // fill up the outgoing sequence and the included reads
  if(extend_direction == EXT_RIGHT)  {
#if DEBUG
    cout << "BuildNewEdge:  *********************************************" << endl;
    cout << "BuildNewEdge:  extending right" << endl;
    cout << "BuildNewEdge:  distance: " << distance << endl;
    cout << "BuildNewEdge:  edge sequence:  " << edge_content.outgoing_seq << endl;
    cout << "BuildNewEdge:  extend sequence:  " << suffix_array_->getSuffix(static_cast<size_t>(extend_path.key_sequence_position)) << endl;
#endif
    
    string extend_seq = suffix_array_->getSuffix(static_cast<size_t>(extend_path.key_sequence_position));
    // construct the new assembled sequence
    int exp_len_extend = edge_content.outgoing_seq.length() + distance;
    if((int) (extend_seq.length() - extend_path.len_search_seed) >= exp_len_extend)  {
      // if the extend sequence contain full sequence required, then take full new edge sequence from extend sequence
      new_edge.outgoing_seq = extend_seq.substr(extend_path.len_search_seed, exp_len_extend);
    } else  {
      // otherwise we have to return, or we will create holes in the assembled sequence
      return;
    }
    
    if(new_edge.outgoing_seq.empty() || !IsSeqEndMatch(new_edge.outgoing_seq, edge_content.outgoing_seq))  {
      return;
    }
#if DEBUG
    cout << "BuildNewEdge:  right new seq: dist:  " << distance << endl;
    cout << "BuildNewEdge:  edge outgoing sequence: " << edge_content.outgoing_seq << endl;
    cout << "BuildNewEdge:  new outgoing sequence:  " << new_edge.outgoing_seq << endl;
#endif

    
    // construct new query sequence
    if(query_end < (int) edge_content.query_end)  {
      query_aln_seq = query_seq_->substr(query_end + 1, edge_content.query_end - query_end);  
    }
    // include the reads in the path to the edge
    for(auto it = extend_path.included_reads.begin(); it != extend_path.included_reads.end(); ++ it) {
      if((unsigned int) suffix_array_->getSuffixLength(static_cast<size_t>(it->rank_in_sfa)) 
          - extend_path.len_search_seed <= new_edge.outgoing_seq.length())  {
        new_edge.included_reads.push_back(*it);
#if DEBUG
        int begin = (unsigned int) it->position + source_alignment.seed_extend_seq.length();
        cout << "BuildNewEdge right:  included seq: " << it->full_seq.substr(begin, it->full_seq.length() - begin) << endl;
#endif
      }
    }

  } else  {
#if DEBUG
    cout << "BuildNewEdge:  *********************************************" << endl;
    cout << "BuildNewEdge:  extending left" << endl;
    cout << "BuildNewEdge:  distance: " << distance << endl;
    cout << "BuildNewEdge:  edge sequence:  " << edge_content.outgoing_seq << endl;
    cout << "BuildNewEdge:  extend sequence:  " << reverse_suffix_array_->getSuffix(static_cast<size_t>(extend_path.key_sequence_position)) << endl;
#endif
    // construct the new assembled sequence
    string extend_seq = reverse_suffix_array_->getSuffix(static_cast<size_t>(extend_path.key_sequence_position));
    int exp_len_extend = edge_content.outgoing_seq.length() + distance;
    if((int) (extend_seq.length() - extend_path.len_search_seed) >= exp_len_extend)  {
      // if the extend sequence contain full sequence required, then take full new edge sequence from extend sequence
      new_edge.outgoing_seq = extend_seq.substr(extend_path.len_search_seed, exp_len_extend);
    } else  {
      // otherwise we need to return, otherwise we will create holes in the assembled sequences
      return;
    }
    if(new_edge.outgoing_seq.empty() || !IsSeqEndMatch(new_edge.outgoing_seq, edge_content.outgoing_seq))  {
      return;
    }
#if DEBUG
    cout << "BuildNewEdge:  left new seq: dist:  " << distance << endl;
    cout << "BuildNewEdge:  edge outgoing sequence: " << edge_content.outgoing_seq << endl;
    cout << "BuildNewEdge:  new outgoing sequence:  " << new_edge.outgoing_seq << endl;
#endif

    
    // construct the new query sequence
    if(query_end > (int) edge_content.query_end)  {
      query_aln_seq = query_seq_->substr(edge_content.query_end, query_end - edge_content.query_end);
      query_aln_seq = ReverseShortSequence(query_aln_seq);
    }
    // include the reads
    for(auto it = extend_path.included_reads.begin(); it != extend_path.included_reads.end(); ++ it) {
      if((unsigned int) reverse_suffix_array_->getSuffixLength(static_cast<size_t>(it->rank_in_sfa)) 
          - extend_path.len_search_seed <= new_edge.outgoing_seq.length())  {
        new_edge.included_reads.push_back(*it);
#if DEBUG 
        int begin = (unsigned int) it->position + source_alignment.seed_extend_seq.length();
        cout << "BuildNewEdge left:  included seq: " << it->full_seq.substr(begin, it->full_seq.length() - begin) << endl;
#endif
      }
    }

  }
  // align the new sequence to estimate the alignment score
  string current_query_seq;
  if(source_alignment.extend_direction == EXT_RIGHT)  {
    current_query_seq = query_seq_->substr(start_, source_alignment.current_query_length);
  } else  {
    current_query_seq = query_seq_->substr(
        start_ + seed_length_ - source_alignment.current_query_length, 
        source_alignment.current_query_length
    );
    current_query_seq = ReverseShortSequence(current_query_seq);
  }
  SeqAlignExtend<AlignScoreType> align_phase(
    current_query_seq, source_alignment.current_assembled_seq,
    query_aln_seq, new_edge.outgoing_seq, source_alignment.edge_score, score_scheme_,
    GLOBAL, right_band_, down_band_, false
  );
  align_phase.Align();
  new_edge.current_score = align_phase.GetBestEdgeScore();
  new_edge.increased_score = new_edge.current_score - source_alignment.max_score;
  new_edge.query_end = edge_content.query_end;
  return;
}
*/

int ExtendAndAssemble::SubmitEdgeReads(
    AlignmentInfoType<AlignScoreType>& current_alignment,
    BoostEdge& current_edge, BoostVertex& current_vertex    
) {
  //cout << "begin of SubmitEdgeReads" << endl;
  assert(current_alignment.vertex_holder.size() == current_alignment.edge_holder.size());
  assert(current_alignment.extend_direction == EXT_RIGHT  || current_alignment.extend_direction == EXT_LEFT);
  
  // if the source of the current extension has not been initialized
  // if the source has not been initialized by its partner who extends the alignment to the opposite direction, initialize it
  
  // add to the assembly graph vertex by vertex
  auto it_vertex = current_alignment.vertex_holder.begin();
  int num_valid_updates = 0;
  for(auto it = current_alignment.edge_holder.begin(); it != current_alignment.edge_holder.end(); ++ it) {
    //cout << "here A" << endl;
    // updating the edge in the graph
    assembly_graph_->UpdateEdge(current_edge, *it);

    // updating the vertex in the graph
    assembly_graph_->UpdateVertex(current_vertex, *it_vertex);

    // update the current alignment information
    current_alignment.source_vertex = current_vertex;
    // go to the next vertex
    ++ it_vertex;
    //cout << "here 0" << endl;
    // submitting the included reads
    
    for(auto it_r = it->included_reads.begin(); it_r != it->included_reads.end(); ++ it_r)  {
      RIDType read_ID = it_r->read_ID;
      ReadTerminateType recorder;
      recorder.right_extended = recorder.left_extended = false;
      recorder.assembled_score_left = recorder.assembled_score_right = -_MAX_VALUE;
      recorder.alignment_position = it_r->alignment_position;
      //recorder.assembled_score = current_alignment.max_score;
      if(it_r->extend_direction == EXT_RIGHT)  {
        recorder.right_extended = true;
        recorder.assembled_score_right = current_alignment.max_score;
        //recorder.end_vertex_right = current_vertex;
        //recorder.end_edge_right = current_edge;
        //recorder.rank_in_sfa_right = it_r->rank_in_sfa;
      } else  {
        recorder.left_extended = true;
        recorder.assembled_score_left = current_alignment.max_score;
        //recorder.end_vertex_left = current_vertex;
        //recorder.end_edge_left = current_edge;
        //recorder.rank_in_sfa_left = it_r->rank_in_sfa;
      }
      //cout << "here C" << endl;
      num_valid_updates += (int) UpdateAssignedReads(read_ID, recorder);
      //cout << "here D" << endl;
    }
    
    //cout << "here B" << endl;
  }
  //cout << "here1" << endl;
  // reinitialize the non_increasing count
  current_alignment.non_increasing_count = 0;
  // clear the lists that hold the vertices and edges
  current_alignment.vertex_holder.clear();
  current_alignment.edge_holder.clear();
  //cout << "End of SubmitEdgeReads" << endl;
  return num_valid_updates;
}

void ExtendAndAssemble::RecordRecruitedReadsLeft(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    const std::list<SfaType>& reads_positions,
    std::list<ReadAssignmentType>& reads_assignments
) {
  for(auto it = reads_positions.begin(); it != reads_positions.end(); ++ it) {
    int mapped_end = ComputeQueryRightBoundLeft(source_alignment, *it);
    RIDType r_id = (RIDType) reverse_suffix_array_->getId(static_cast<size_t>(*it)); 
    int s_len = (int) reverse_suffix_array_->getSuffixLength(static_cast<size_t>(*it));
    //cout << " AddRecruitedReadsLeft:: updating hash table list:  " << r_id << "  " << mapped_end << endl;
    ReadAssignmentType mapping;
    mapping.read_ID = r_id;
    mapping.position = (POSType) reverse_suffix_array_->getPos(static_cast<size_t>(*it));
    //mapping.full_seq_len = reverse_suffix_array_->getFullSequenceLength(static_cast<size_t>(*it));
    mapping.len_search_seed = source_alignment.seed_extend_seq.length();
    mapping.alignment_position = mapped_end;
    mapping.rank_in_sfa = *it;
    mapping.offset_to_start = 
        -(source_alignment.current_assembled_seq.length() - seed_length_ + s_len - source_alignment.seed_extend_seq.length());
        
    /************************************************************************************
    string read_seq = reverse_suffix_array_->getSequence_explicit((unsigned int) r_id);
    int seg_len = read_seq.length() - s_len + source_alignment.seed_extend_seq.length();
    string seq_read = read_seq.substr(0, seg_len);
    string seq_assemble = source_alignment.current_assembled_seq.substr(seed_length_ - mapping.offset_to_start - read_seq.length(), seg_len);
    if(seq_read != seq_assemble)  {
      cout << "@@@: fetched seq:  " << seq_read << "  " << seq_assemble << endl;
      cout << "@@@: assembled:  " << source_alignment.current_assembled_seq << endl;
      cout << "@@@: seed extend:  " << source_alignment.seed_extend_seq << endl;
      cout << "@@@: suffix length: " << (int) s_len << endl; 
      cout << "@@@: full read:  " << reverse_suffix_array_->getSequence_explicit((unsigned int) r_id) << endl;
      cout << "@@@: computed: " << mapping.offset_to_start << endl;
    }
    ************************************************************************************/
      
    mapping.extend_direction = source_alignment.extend_direction;
    //if(!HasAssignedLeft(source_alignment, *it))  {
      //cout << " AddRecruitedReadsLeft:: adding new element:" << endl;
    reads_assignments.push_back(mapping);
    //}
  }
  return;
}

void ExtendAndAssemble::RecordCandidateReadsLeft(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    const std::list<SfaType>& candidate_positions,
    std::list<BridgingReadType>& candidate_bridging_reads
) {
  for(auto it = candidate_positions.begin(); it != candidate_positions.end(); ++ it) {
    RIDType r_id = (RIDType) reverse_suffix_array_->getId(static_cast<size_t>(*it));
    POSType s_len = (int) reverse_suffix_array_->getSuffixLength(static_cast<size_t>(*it));
    //POSType pos = (POSType) reverse_suffix_array_->getPos(static_cast<size_t>(*it)); 
    BridgingReadType candidate;
    candidate.read_ID = r_id;
    candidate.position = (POSType) reverse_suffix_array_->getPos(static_cast<size_t>(*it));
    //candidate.full_seq_len = reverse_suffix_array_->getFullSequenceLength(static_cast<size_t>(*it));
    candidate.len_search_seed = source_alignment.seed_extend_seq.length();
    candidate.rank_in_sfa = *it;
    candidate.offset_to_start = 
        -(source_alignment.current_assembled_seq.length() - seed_length_ + s_len - source_alignment.seed_extend_seq.length());
    
    /**********************************************************************  
    int seg_len = source_alignment.current_assembled_seq.length();
    string read_seq = reverse_suffix_array_->getSequence_explicit((unsigned int) r_id);
    string seq_read = read_seq.substr(read_seq.length() + candidate.offset_to_start - seg_len, seg_len);
    if(seq_read != source_alignment.current_assembled_seq)  {
      cout << "@@@: fetched seq:  " << seq_read << "  " << source_alignment.current_assembled_seq << endl;
      cout << "@@@: assembled:  " << source_alignment.current_assembled_seq << endl;
      cout << "@@@: seed extend:  " << source_alignment.seed_extend_seq << endl;
      cout << "@@@: suffix length: " << (int) s_len << endl; 
      cout << "@@@: full read:  " << reverse_suffix_array_->getSequence_explicit((unsigned int) r_id) << endl;
      cout << "@@@: computed: " << candidate.offset_to_start << endl;
    }
    **********************************************************************/    
    candidate.extend_direction = source_alignment.extend_direction;
    candidate.alignment_position = ComputeBridgingReadsRightBound(candidate.rank_in_sfa, candidate.extend_direction);
    candidate_bridging_reads.push_back(candidate);
  }
  return;
}

void ExtendAndAssemble::RecordRecruitedReadsRight(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    const std::list<SfaType>& reads_positions,
    std::list<ReadAssignmentType>& reads_assignments
) {
  // compute the index of the ending in the query sequence, note that the query_end
  // may be shorter than the sequence length
  for(auto it = reads_positions.begin(); it != reads_positions.end(); ++ it) {
    int mapped_end = ComputeQueryRightBoundRight(source_alignment, *it);
    RIDType r_id = (RIDType) suffix_array_->getId(static_cast<size_t>(*it));
    POSType pos = (POSType) suffix_array_->getPos(static_cast<size_t>(*it));
    //cout << " AddRecruitedReadsRight:: updating hash table list:  " << r_id << "  " << mapped_end << endl;
    ReadAssignmentType mapping;
    mapping.read_ID = r_id;
    mapping.position = pos;
    mapping.alignment_position = mapped_end;
    //mapping.full_seq_len = suffix_array_->getFullSequenceLength(static_cast<size_t>(*it));
    mapping.len_search_seed = source_alignment.seed_extend_seq.length();   
    mapping.rank_in_sfa = *it;
    mapping.offset_to_start = source_alignment.current_assembled_seq.length() - source_alignment.seed_extend_seq.length() - (int) pos;
    
    
    mapping.extend_direction = source_alignment.extend_direction;
    //if(!HasAssignedRight(source_alignment, *it))  {
    reads_assignments.push_back(mapping);
    //}
  }
  return;
}

void ExtendAndAssemble::RecordCandidateReadsRight(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    const std::list<SfaType>& candidate_positions,
    std::list<BridgingReadType>& candidate_bridging_reads
) {
  for(auto it = candidate_positions.begin(); it != candidate_positions.end(); ++ it) {
    RIDType r_id = (RIDType) suffix_array_->getId(static_cast<size_t>(*it));
    POSType pos = (POSType) suffix_array_->getPos(static_cast<size_t>(*it));
    BridgingReadType candidate;
    candidate.read_ID = r_id;
    candidate.position = pos;
    //candidate.full_seq_len = suffix_array_->getFullSequenceLength(static_cast<size_t>(*it));
    candidate.len_search_seed = source_alignment.seed_extend_seq.length();
    candidate.rank_in_sfa = *it;
    candidate.offset_to_start = source_alignment.current_assembled_seq.length() - source_alignment.seed_extend_seq.length() - (int) pos;
    
    candidate.extend_direction = source_alignment.extend_direction;
    candidate.alignment_position = ComputeBridgingReadsRightBound(candidate.rank_in_sfa, candidate.extend_direction);
    candidate_bridging_reads.push_back(candidate);
  }
  return;
}

bool ExtendAndAssemble::HasAssignedLeft(
  const AlignmentInfoType<AlignScoreType>& source_alignment,
  SfaType extend_position
) {
	int q_end = ComputeQueryRightBoundLeft(source_alignment, extend_position);
	if(q_end > (int) query_seq_->length() - 1)  {
	  q_end = query_seq_->length() - 1;
	}
	
	RIDType read_ID = (RIDType) reverse_suffix_array_->getId(static_cast<size_t>(extend_position));
	
	SLock r_lock(*mutex_assigned_reads_);
	auto it = assigned_reads_->find(read_ID);
	if(it != assigned_reads_->end())  {
	  auto it_h = clump_map_->find(q_end); 
	  assert(it_h != clump_map_->end());
	  int idx = it_h->second;
	  auto it_map = it->second.find(idx);
	  if(it_map != it->second.end())  {
	    r_lock.unlock();
	    return true;  // ID found from hash table and is aligned nearby
	  }
	  r_lock.unlock();
	  return false; // ID found from hash table but not aligned nearby
	} 
	r_lock.unlock();
	return false;
}

bool ExtendAndAssemble::HasAssignedRight(
  const AlignmentInfoType<AlignScoreType>& source_alignment,
  SfaType extend_position
) {
	int q_end = ComputeQueryRightBoundRight(source_alignment, extend_position);
	if(q_end > (int) query_seq_->length() - 1)  {
	  q_end = query_seq_->length() - 1;
	}
	RIDType read_ID = (RIDType) suffix_array_->getId(static_cast<size_t>(extend_position));
	
	SLock r_lock(*mutex_assigned_reads_);
	auto it = assigned_reads_->find(read_ID);
	if(it != assigned_reads_->end())  {
	  auto it_h = clump_map_->find(q_end); 
	  assert(it_h != clump_map_->end());
	  int idx = it_h->second;
	  auto it_map = it->second.find(idx);
	  if(it_map != it->second.end())  {
	    r_lock.unlock();
	    return true;  // ID found from hash table and is aligned nearby
	  }
	  r_lock.unlock();
	  return false; // ID found from hash table but not aligned nearby
	}
	r_lock.unlock();
  return false;
}

string ExtendAndAssemble::GetQueryExtendSeqLeft(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    SfaType extend_position
) {
  int str_start = ComputeQueryLeftBoundLeft(source_alignment, extend_position);
  if(str_start < 0)  {
    str_start = 0;
  }
  //int str_end = (start_ + seed_length_ - 1) - source_alignment.current_query_seq.length();
	int str_end = (start_ + seed_length_ - 1) - source_alignment.current_query_length;
	//cout << "GetQueryExtendSeqLeft: start and len:  " << str_start << " " << str_len << endl;
	string seq_holder =  query_seq_->substr(str_start, str_end - str_start + 1);
	return ReverseShortSequence(seq_holder);
}


string ExtendAndAssemble::GetQueryExtendSeqRight(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    SfaType extend_position
) {
  int str_end = ComputeQueryRightBoundRight(source_alignment, extend_position);
  //int str_start = start_ + source_alignment.current_query_seq.length();
  int str_start = start_ + source_alignment.current_query_length; 
  return query_seq_->substr(str_start, str_end - str_start + 1);
}

void ExtendAndAssemble::ComputeOutgoingPathsLeft(
    AlignmentInfoType<AlignScoreType>& source_alignment,
    std::list<OutgoingPathType>& out_paths
) {
  int i;  
  // search the seed sequence from the suffix array
  string seed_seq = source_alignment.seed_extend_seq;
  string assembled_seq = source_alignment.current_assembled_seq;
#if DEBUG  
  cout << " ComputeOutgoingpathsLeft:: seed sequence: " << seed_seq << endl;
#endif  
  BoundType range = reverse_suffix_array_->searchWithLCPs((SfaChar*) seed_seq.c_str(), seed_seq.length());
  while(range.first <= range.second &&
      (unsigned int) reverse_suffix_array_->getSuffixLength(static_cast<size_t>(range.first)) <= seed_seq.length()) {
    ++ range.first;
  }
  if(range.first > range.second)  {
    return;
  }
  //cout << " ComputeOutgoingpathsLeft:: range 1: " << range.first << "  " << range.second << endl;
  // note that the suffix array search will return the seed sequence itself, 
  // need to jump forward to the first extension sequence
  
  unordered_map<unsigned int, bool> suffix_valid_indicator;
  // if the seed is incomplete, then it means the assembled sequence is not long enough
  // to carry out the n_back_check, and therefore we cannot confidently include any reads
  // In this case, mark all reads invalid
  if(source_alignment.is_seed_complete)  {
    SelectValidSuffixLeft(assembled_seq, seed_seq, range, suffix_valid_indicator);
  } else  {
    for(i = range.first; i <= range.second; ++ i) {
      suffix_valid_indicator[i] = true;
    }
  }
 
  unsigned int position;
  position = reverse_suffix_array_->getPos(static_cast<size_t>(range.first));
  while(range.first <= range.second && !suffix_valid_indicator[range.first]) {
#if DEBUG
    cout << "ComputeOutgoingPathsLeft: position: " << position << endl;
    cout << "ComputeOutgoingPathsLeft: seed_len: " << seed_seq.length() << endl;
    cout << "ComputeOutgoingPathsLeft: assembled_len:  " << assembled_seq.length() << endl;
#endif
    ++ range.first;
    position = reverse_suffix_array_->getPos(static_cast<size_t>(range.first));
  }
  if(range.first > range.second)  {
    return;
  }
#if DEBUG
  cout << " ComputeOutgoingpathsLeft:: range final: " << range.first << "  " << range.second << endl;
#endif  
  
  IndexSample index_caller;
  
  list<SfaType> reads_recorder;     // reads that are put in this list indicate reads that can be confidently placed
  list<SfaType> candidate_recorder; // reads that are put in this list indicate reads that need to be checked when combining both extension directions (bridging reads)
  LcpType prev_lcp = reverse_suffix_array_->getSuffixLength(static_cast<size_t>(range.first));
  position = (unsigned int) reverse_suffix_array_->getPos(static_cast<size_t>(range.first));
  LcpType min_lcp_between_valid = prev_lcp;
  SfaType last_valid_seq = range.first;
  if(suffix_valid_indicator[range.first] && position + seed_seq.length() <= assembled_seq.length())  {
    reads_recorder.push_back(range.first);
  } else if (suffix_valid_indicator[range.first] && position + seed_seq.length() > assembled_seq.length()) {
    candidate_recorder.push_back(range.first);
  }
  for(i = range.first + 1; i <= range.second; ++ i) {
    LcpType adjacent_lcp = reverse_suffix_array_->getLcp(static_cast<size_t>(i));
    min_lcp_between_valid = min_lcp_between_valid < adjacent_lcp ? min_lcp_between_valid : adjacent_lcp;
    position = reverse_suffix_array_->getPos(static_cast<size_t>(i));
    // check if suffix i is a valid sequence to look at
    
    if(!suffix_valid_indicator[i])  {
      continue;
    }
    
    //cout << " ComputeOutgoingpathsLeft:: range passed and checked: " << i << endl;
    // if suffix i is a valid sequence to look at
    LcpType current_lcp = min_lcp_between_valid;
    if(current_lcp < reverse_suffix_array_->getSuffixLength(static_cast<size_t>(last_valid_seq)) || current_lcp < prev_lcp)  {
      OutgoingPathType single_outgoing_path;
      single_outgoing_path.extend_direction = EXT_LEFT;
      single_outgoing_path.key_sequence_position = last_valid_seq;
      single_outgoing_path.len_search_seed = seed_seq.length();
      /* TODO: obsolete part, should be removed after testing
      for(auto it_reads_recorder = reads_recorder.begin(); it_reads_recorder != reads_recorder.end(); ++ it_reads_recorder) {
        //single_outgoing_path.included_reads.push_back(suffix_array_->getId(static_cast<size_t>(*it_reads_recorder)));
        single_outgoing_path.included_reads.push_back(*it_reads_recorder);
      }
      */
      RecordRecruitedReadsLeft(source_alignment, reads_recorder, single_outgoing_path.included_reads);
      RecordCandidateReadsLeft(source_alignment, candidate_recorder, single_outgoing_path.candidate_bridging_reads);
      out_paths.push_back(single_outgoing_path);
      // We should not clear the entire reads_recorder, we need to maintain
      // the reads that are shorter than the current_lcp, as they might also
      // contribute to the current outgoing path
      //cout << " ComputeOutgoingPathsRight:: before pop:  " << reads_recorder.size() << endl;
      while(!reads_recorder.empty()) {
        SfaType last_pos = *(std::prev(reads_recorder.end()));
        unsigned int r_id = reverse_suffix_array_->getId(static_cast<size_t>(last_pos));
        unsigned int r_pos = reverse_suffix_array_->getPos(static_cast<size_t>(last_pos));
        if(reverse_suffix_array_->getSeqLength_RID(r_id) - r_pos > current_lcp)  {
          reads_recorder.pop_back();
        } else  {
          break;
        }
      }
      while(!candidate_recorder.empty()) {
        SfaType last_pos = *(std::prev(candidate_recorder.end()));
        unsigned int r_id = reverse_suffix_array_->getId(static_cast<size_t>(last_pos));
        unsigned int r_pos = reverse_suffix_array_->getPos(static_cast<size_t>(last_pos));
        if(reverse_suffix_array_->getSeqLength_RID(r_id) - r_pos > current_lcp)  {
          candidate_recorder.pop_back();
        } else  {
          break;
        }
      }
      //cout << " ComputeOutgoingPathsRight:: after pop: " << reads_recorder.size() << endl;
    } 
    // pushing the current suffix to the stack
    //ReadMapType read_to_be_added;
    //read_to_be_added.read_position = index_caller.encode_read_position(
    //    static_cast<RIDType>(suffix_array_->getId(i)), 
    //    static_cast<POSType>(suffix_array_->getPos(i))
    //);
    //read_to_be_added.aligned_position_assembled = 0;
    if(suffix_valid_indicator[i] && position + seed_seq.length() <= assembled_seq.length())  {
      reads_recorder.push_back(i);
    } else if(suffix_valid_indicator[i] && position + seed_seq.length() > assembled_seq.length())  {
      candidate_recorder.push_back(i);
    }
    prev_lcp = current_lcp;
    min_lcp_between_valid = reverse_suffix_array_->getSuffixLength(static_cast<size_t>(i));
    last_valid_seq = i;
  }
  // record last_valid_seq into the recruited reads
  OutgoingPathType last_outgoing_path;
  last_outgoing_path.extend_direction = EXT_LEFT;
  last_outgoing_path.key_sequence_position = last_valid_seq;
  last_outgoing_path.len_search_seed = seed_seq.length();
  RecordRecruitedReadsLeft(source_alignment, reads_recorder, last_outgoing_path.included_reads);
  RecordCandidateReadsLeft(source_alignment, candidate_recorder, last_outgoing_path.candidate_bridging_reads);
  /* TODO: the following segment of code is obsolete. remove after testing
  for(auto it_reads_recorder = reads_recorder.begin(); it_reads_recorder != reads_recorder.end(); ++ it_reads_recorder) {
    //single_outgoing_path.included_reads.push_back(suffix_array_->getId(static_cast<size_t>(*it_reads_recorder)));
    last_outgoing_path.included_reads.push_back(*it_reads_recorder);
  }
  */
  out_paths.push_back(last_outgoing_path);
  
  //for(auto it = out_paths.begin(); it != out_paths.end(); ++ it) {
  //  unsigned int r_id = suffix_array_->getId(static_cast<size_t>(it->key_sequence_position));
  //  unsigned int r_pos = 
  //      suffix_array_->getPos(static_cast<size_t>(it->key_sequence_position)) + it->len_search_seed;
  //  string out_sequence = string(suffix_array_->getSuffix_explicit(r_id, r_pos)); 
  //  cout << " ComputeOutgoingPathsRight:: out sequence:  " << out_sequence << endl;
  //  for(auto it_list = it->included_reads.begin(); it_list != it->included_reads.end(); ++ it_list) {
  //    cout << *it_list << endl;
  //  }
  //}
  return;
}

void ExtendAndAssemble::ComputeOutgoingPathsRight(
    AlignmentInfoType<AlignScoreType>& source_alignment,
    std::list<OutgoingPathType>& out_paths
)  {
  int i;
  // search the seed sequence from the suffix array
  string seed_seq = source_alignment.seed_extend_seq;
  string assembled_seq = source_alignment.current_assembled_seq;
#if DEBUG  
  cout << "ComputeOutgoingpathsRight:: seed sequence: " << seed_seq << endl;
#endif  
  BoundType range = suffix_array_->searchWithLCPs((SfaChar*) seed_seq.c_str(), seed_seq.length());
  while(range.first <= range.second &&
      (unsigned int) suffix_array_->getSuffixLength(static_cast<size_t>(range.first)) <= seed_seq.length()) {
    ++ range.first;
  }
  if(range.first > range.second)  {
    return;
  }
  //cerr << " ComputeOutgoingpathsRight:: range 1: " << range.first << "  " << range.second << "  seq:  " << seed_seq << endl;
  // note that the suffix array search will return the seed sequence itself, 
  // need to jump forward to the first extension sequence
  
  unordered_map<unsigned int, bool> suffix_valid_indicator;
  // if the seed is incomplete, then it means the assembled sequence is not long enough
  // to carry out the n_back_check, and therefore we cannot confidently include any reads
  // In this case, mark all reads invalid
  if(source_alignment.is_seed_complete)  {
    SelectValidSuffixRight(assembled_seq, seed_seq, range, suffix_valid_indicator);
  } else  {
    for(i = range.first; i <= range.second; ++ i) {
      suffix_valid_indicator[i] = true;
    }
  }
  // need to jump to the first valid outgoing path
  //cout << " ComputeOutgoingPathsRight:: before the while loop" << endl;
  //cout << " ComputeOutgoingPathsRight:: printing the suffix_valid_indicator" << endl;
  /********************************************************************************************************
  for(auto it_v = suffix_valid_indicator.begin(); it_v != suffix_valid_indicator.end(); ++ it_v) {
    int non_overlap_suffix_begin = suffix_array_->getPos(it_v->first) + source_alignment.seed_extend_seq.length();
    int max_len = non_overlap_suffix_begin > source_alignment.current_assembled_seq.length() ? 
        source_alignment.current_assembled_seq.length() : non_overlap_suffix_begin;
    int asm_len = source_alignment.current_assembled_seq.length();
    string asm_seq = source_alignment.current_assembled_seq.substr(asm_len - max_len, max_len);
    string full_seq = suffix_array_->getSequence_explicit(suffix_array_->getId(it_v->first));
    string sfa_seq = full_seq.substr(non_overlap_suffix_begin - max_len, max_len);
    cout << "$$$  asm_seq: " << asm_seq << " sfa_seq:  " << sfa_seq << endl;
    cout << "$$$  valid:  " << it_v->second << endl;
    if((asm_seq != sfa_seq && it_v->second) || (asm_seq == sfa_seq && !it_v->second))  {
      cout << "$$$  Error" << endl;
    }
  }
  **********************************************************************************************************/
  unsigned int position;
  position = suffix_array_->getPos(static_cast<size_t>(range.first));
  while(range.first <= range.second && !suffix_valid_indicator[range.first]) {
#if DEBUG
    cout << "ComputeOutgonigPathsRight: position: " << position << endl;
    cout << "ComputeOutgoingPathsRight: seed_len: " << seed_seq.length() << endl;
    cout << "ComputeOutgoingPathsRight: assembled_len:  " << assembled_seq.length() << endl;
#endif
    ++ range.first;
    position = suffix_array_->getPos(static_cast<size_t>(range.first));
  }
  if(range.first > range.second)  {
    return;
  }
#if DEBUG
  cout << "ComputeOutgoingPathsRight:: range final: " << range.first << "  " << range.second << endl;
#endif  
  
  IndexSample index_caller;
  
  list<SfaType> reads_recorder;
  list<SfaType> candidate_recorder;
  LcpType prev_lcp = suffix_array_->getSuffixLength(static_cast<size_t>(range.first));
  LcpType min_lcp_between_valid = prev_lcp;
  position = (unsigned int) suffix_array_->getPos(static_cast<size_t>(range.first));
  SfaType last_valid_seq = range.first;
  //cerr << "ComputeOutgoingPathsRight:: last_valid_seq " << last_valid_seq << endl;
  if(suffix_valid_indicator[range.first] && position + seed_seq.length() <= assembled_seq.length())  {
#if DEBUG
    cout << "ComputeOutgoingPathsRight:: first position " << range.first << " is valid" << endl;
#endif
    reads_recorder.push_back(range.first);
  } else if (suffix_valid_indicator[range.first] && position + seed_seq.length() > assembled_seq.length()) {
#if DEBUG
    cout << "ComputeOutgoingPathsRight:: first position " << range.first << " is Invalid" << endl;
#endif
    candidate_recorder.push_back(range.first);
  }
  for(i = range.first + 1; i <= range.second; ++ i) {
    LcpType adjacent_lcp = suffix_array_->getLcp(static_cast<size_t>(i));
    min_lcp_between_valid = min_lcp_between_valid < adjacent_lcp ? min_lcp_between_valid : adjacent_lcp;
    position = suffix_array_->getPos(static_cast<size_t>(i));
    // check if suffix i is a valid sequence to look at
    if(!suffix_valid_indicator[i])  {
      continue;
    }
    //cout << " ComputeOutgoingpathsRight:: range passed and checked: " << i << endl;
    // if suffix i is a valid sequence to look at
    LcpType current_lcp = min_lcp_between_valid;
    if(current_lcp < suffix_array_->getSuffixLength(static_cast<size_t>(last_valid_seq)) || current_lcp < prev_lcp)  {
      OutgoingPathType single_outgoing_path;
      single_outgoing_path.extend_direction = EXT_RIGHT;
      single_outgoing_path.key_sequence_position = last_valid_seq;
      single_outgoing_path.len_search_seed = seed_seq.length();
      /* TODO: obsolete version, should be removed after testing
      for(auto it_reads_recorder = reads_recorder.begin(); it_reads_recorder != reads_recorder.end(); ++ it_reads_recorder) {
        //single_outgoing_path.included_reads.push_back(suffix_array_->getId(static_cast<size_t>(*it_reads_recorder)));
        single_outgoing_path.included_reads.push_back(*it_reads_recorder);
      }
      */
      RecordRecruitedReadsRight(source_alignment, reads_recorder, single_outgoing_path.included_reads);
      RecordCandidateReadsRight(source_alignment, candidate_recorder, single_outgoing_path.candidate_bridging_reads);
      out_paths.push_back(single_outgoing_path);
      // We should not clear the entire reads_recorder, we need to maintain
      // the reads that are shorter than the current_lcp, as they might also
      // contribute to the current outgoing path
      //cout << " ComputeOutgoingPathsRight:: before pop:  " << reads_recorder.size() << endl;
      while(!reads_recorder.empty()) {
        SfaType last_pos = *(std::prev(reads_recorder.end()));
        unsigned int r_id = suffix_array_->getId(static_cast<size_t>(last_pos));
        unsigned int r_pos = suffix_array_->getPos(static_cast<size_t>(last_pos));
        if(suffix_array_->getSeqLength_RID(r_id) - r_pos > current_lcp)  {
          reads_recorder.pop_back();
        } else  {
          break;
        }
      }
      while(!candidate_recorder.empty()) {
        SfaType last_pos = *(std::prev(candidate_recorder.end()));
        unsigned int r_id = suffix_array_->getId(static_cast<size_t>(last_pos));
        unsigned int r_pos = suffix_array_->getPos(static_cast<size_t>(last_pos));
        if(suffix_array_->getSeqLength_RID(r_id) - r_pos > current_lcp)  {
          candidate_recorder.pop_back();
        } else  {
          break;
        }
       }
      //cout << " ComputeOutgoingPathsRight:: after pop: " << reads_recorder.size() << endl;
    } 
    // pushing the current suffix to the stack
    //ReadMapType read_to_be_added;
    //read_to_be_added.read_position = index_caller.encode_read_position(
    //    static_cast<RIDType>(suffix_array_->getId(i)), 
    //    static_cast<POSType>(suffix_array_->getPos(i))
    //);
    //read_to_be_added.aligned_position_assembled = 0;
    if(suffix_valid_indicator[i] && position + seed_seq.length() <= assembled_seq.length())  {
      reads_recorder.push_back(i);
    } else if(suffix_valid_indicator[i] && position + seed_seq.length() > assembled_seq.length())  {
      candidate_recorder.push_back(i);
      //cerr << "#############" << endl;
      //cerr << seed_seq << endl;
      //cerr << assembled_seq << endl;
      //unsigned int rid = (unsigned int) suffix_array_->getId(static_cast<size_t>(i));
      //cerr << suffix_array_->getSequence_explicit(rid) << endl;
    }
    prev_lcp = current_lcp;
    min_lcp_between_valid = suffix_array_->getSuffixLength(static_cast<size_t>(i));
    last_valid_seq = i;
  }
  OutgoingPathType last_outgoing_path;
  last_outgoing_path.extend_direction = EXT_RIGHT;
  last_outgoing_path.key_sequence_position = last_valid_seq;
  //cerr << "ComputeOutgoingPathsRight:: last_valid_seq2 " << last_valid_seq << endl;
  last_outgoing_path.len_search_seed = seed_seq.length();
  RecordRecruitedReadsRight(source_alignment, reads_recorder, last_outgoing_path.included_reads);
  RecordCandidateReadsRight(source_alignment, candidate_recorder, last_outgoing_path.candidate_bridging_reads);
  /* TODO: the following segment of code is obsolete. remove after test
  for(auto it_reads_recorder = reads_recorder.begin(); it_reads_recorder != reads_recorder.end(); ++ it_reads_recorder) {
    //single_outgoing_path.included_reads.push_back(suffix_array_->getId(static_cast<size_t>(*it_reads_recorder)));
    last_outgoing_path.included_reads.push_back(*it_reads_recorder);
  }
  */
  out_paths.push_back(last_outgoing_path);
  
  //for(auto it = out_paths.begin(); it != out_paths.end(); ++ it) {
  //  unsigned int r_id = suffix_array_->getId(static_cast<size_t>(it->key_sequence_position));
  //  unsigned int r_pos = 
  //      suffix_array_->getPos(static_cast<size_t>(it->key_sequence_position)) + it->len_search_seed;
  //  string out_sequence = string(suffix_array_->getSuffix_explicit(r_id, r_pos)); 
  //  cout << " ComputeOutgoingPathsRight:: out sequence:  " << out_sequence << endl;
  //  for(auto it_list = it->included_reads.begin(); it_list != it->included_reads.end(); ++ it_list) {
  //    cout << *it_list << endl;
  //  }
  //}
  return;
}

int ExtendAndAssemble::GetKmerCoverageLeft(std::string kmer) {
  assert(suffix_array_ != NULL);
  BoundType range = reverse_suffix_array_->searchWithLCPs((SfaChar*) kmer.c_str(), kmer.length());
  return (range.second - range.first + 1);
}

int ExtendAndAssemble::GetKmerCoverageRight(std::string kmer) {
  assert(suffix_array_ != NULL);
  BoundType range = suffix_array_->searchWithLCPs((SfaChar*) kmer.c_str(), kmer.length());
  return (range.second - range.first + 1);
}

void ExtendAndAssemble::ComputeAlignmentPriority(AlignmentInfoType<AlignScoreType> *a)  {
  a->priority_score = a->seed_coverage;
  return;
}

bool __cmp_align_priority(
    AlignmentInfoType<AlignScoreType> align_info, 
    PriorityType priority_score
) {
  if(align_info.priority_score > priority_score)  {
    return true;
  } else  {
    return false;
  }
}

void ExtendAndAssemble::InsertToAlignmentHolderRight(const AlignmentInfoType<AlignScoreType>* a)  {
  // finds the correct position to insert the alignment "a"
  //SLock r_lock(mutex_holder_queue_right_);
  //auto to_insert = lower_bound(
  //    high_score_alignments_right_.begin(), 
  //    high_score_alignments_right_.end(), 
  //    a->priority_score,
  //    __cmp_align_priority
  //);
  //r_lock.unlock();
  
  //ULock w_lock(mutex_holder_queue_right_);
  //high_score_alignments_right_.insert(to_insert, *a);
  high_score_alignments_right_.push_front(*a);
  //w_lock.unlock();
  
  return;
}

void ExtendAndAssemble::InsertToAlignmentHolderLeft(const AlignmentInfoType<AlignScoreType>* a)  {
  // finds the correct position to insert the alignment "a"
  //SLock r_lock(mutex_holder_queue_right_);
  //auto to_insert = lower_bound(
  //    high_score_alignments_left_.begin(), 
  //    high_score_alignments_left_.end(), 
  //    a->priority_score,
  //    __cmp_align_priority
  //);
  //r_lock.unlock();
  
  //ULock w_lock(mutex_holder_queue_left_);
  //high_score_alignments_left_.insert(to_insert, *a);
  high_score_alignments_left_.push_front(*a);
  //w_lock.unlock();
  
  return;
}

// Selecting valid reads. For each read to be considered as valid, it must be fully compatible with 
// the assembled sequence (meaning that their sequences are exactly the same). If the prefix of the read
// (taking the seed as the center) happen to be longer than the assembled sequence, it is not taken as
// a valid read, but as a candidate bridging read and will be processed in later processes

void ExtendAndAssemble::SelectValidSuffixLeft(
      const string& assembled_seq,
      const string& seed_seq,
      BoundType seed_range,
      unordered_map<unsigned int, bool>& valid_indicator
) {
#if DEBUG
  cout << "   SelectValidSuffixLeft:: selecting valid sequence..." << endl;
#endif
  if(seed_range.first < 0 || seed_range.second < 0 ||
      seed_range.first >= reverse_suffix_array_->getSize() ||
      seed_range.second >= reverse_suffix_array_->getSize())  {
    return;
  }
  // first check the length of all prefixs that ends at seed_alignment.seed_extend_seq
  typedef int SFAIndexType;
  typedef unsigned int SeqLengthType;
  unsigned int i;
  SeqLengthType max_prefix_len = 0;
  for(i = static_cast<unsigned int>(seed_range.first); i <= static_cast<unsigned int>(seed_range.second); ++ i) {
    SeqLengthType prefix_len = reverse_suffix_array_->getPos(static_cast<size_t>(i));
    max_prefix_len = prefix_len > max_prefix_len ? prefix_len : max_prefix_len;
  }
  // search the reverse suffix array with the reversed seed_seq to define upper and lower bounds
  string rev_seed_seq(seed_seq.rbegin(), seed_seq.rend());
  BoundType rev_seed_range = suffix_array_->searchWithLCPs(
      (SfaChar*) rev_seed_seq.c_str(), rev_seed_seq.length()
  );
  
  //while(rev_seed_range.first <= rev_seed_range.second &&
  //  (unsigned int)suffix_array_->getSuffixLength(static_cast<size_t>(rev_seed_range.first)) <= seed_seq.length()) {
  //  ++ rev_seed_range.first;
  //}
#if DEBUG  
  cout << "   SelectValidSuffixLeft:: suffix array range:  " << seed_range.first << "\t" << seed_range.second << endl;
  cout << "   SelectValidSuffixLeft:: prefix array range:  " << rev_seed_range.first << "\t" << rev_seed_range.second << endl;
  cout << "   SelectValidSuffixLeft:: max prefix length: " << max_prefix_len << endl;
#endif
  // identify the position to start searching upward
  SFAIndexType position_search_upward = rev_seed_range.second;
  SeqLengthType prefix_array_search_len = seed_seq.length();
  
  // note that the assembled_sequence should contain the seed_seq as its suffix
  if(max_prefix_len >= assembled_seq.length() - seed_seq.length())  {
    //cout << "   SelectValidSuffixLeft:: When prefix longer than or equal to assembled sequence" << endl;
    prefix_array_search_len = assembled_seq.length();
    string search_seq(assembled_seq.rbegin(), assembled_seq.rend());
    //cout << "   SelectValidSuffixLeft:: assembled_seq:  " << assembled_seq << endl;
    BoundType rev_range = suffix_array_->searchWithLCPs(
        (SfaChar*) search_seq.c_str(), search_seq.length()
    );
    //cout << "   SelectValidSuffixLeft:: final range:  " << rev_range.first << "  " << rev_range.second << endl;
    if(rev_range.second >= rev_range.first)  {
      position_search_upward = rev_range.second;
    } else  {
      position_search_upward = -1;
    }
    //cout << "   SelectValidSuffixLeft:: position determiend: " << position_search_upward << endl;  
  } else  {
    //cout << "   SelectValidSuffixLeft:: When prefix shorter than assembled sequence" << endl;
    SeqLengthType left_bound = seed_seq.length();
    SeqLengthType right_bound = seed_seq.length() + max_prefix_len;
    right_bound = right_bound <= assembled_seq.length() ? right_bound : assembled_seq.length();
    prefix_array_search_len = right_bound;
    BoundType rev_range = rev_seed_range;
    BoundType search_range = rev_seed_range;
    //SeqLengthType sfa_lower_seq_len = suffix_array_->getSuffixLength(rev_seed_range.second);
    while(left_bound <= right_bound)  { 
      // a binary search style recursion to figure out the upward search position
      string search_seq(assembled_seq.rbegin(), assembled_seq.rbegin() + prefix_array_search_len);
      //cout << "   SelectValidSuffixLeft:: sequence to search:  " << search_seq << endl; 
      rev_range = suffix_array_->searchWithLCPs(
          (SfaChar*) search_seq.c_str(), search_seq.length()
      );
      //sfa_lower_seq_len = suffix_array_->getSuffixLength(rev_range.second);
      if(rev_range.second >= rev_range.first)  {
        // the pattern is found but with more than 1 hit, need to make the length longer
        position_search_upward = rev_range.second;
        //cout << "   SelectValidSuffixLeft:: in binary search: pattern found, extending sequence" << endl;
        left_bound = prefix_array_search_len + 1;
        prefix_array_search_len = floor((left_bound + right_bound) / 2);  // record the last sucessful match
        search_range.first = rev_range.first + 1;
        // TODO: determine the suffix array search behavior when a search fails
      } else if(rev_range.second < rev_range.first) {
        // the pattern is not found, need to shorten the length to half
        //cout << "   SelectValidSuffixLeft:: in binary search: pattern not found, shortening sequence" << endl;
        right_bound = prefix_array_search_len - 1;
        prefix_array_search_len = floor((left_bound + right_bound) / 2);
        search_range.second = search_range.second > rev_range.second ? rev_range.second : search_range.second;
      }
    } 
    //cout << "   SelectValidSuffixLeft:: position determiend: " << position_search_upward << endl;  
  } 
  // intersect the read IDs in the defined boundaries
  BoundType rev_lookup_range = rev_seed_range;
  rev_lookup_range.second = position_search_upward;
  //cout << "   SelectValidSuffixLeft: rev_lookup_range:  " << position_search_upward << endl;
  IntersectBoundRanges(
      reverse_suffix_array_, seed_range, 
      suffix_array_, rev_lookup_range, 
      seed_seq.length(), n_back_check_, assembled_seq.length(),
      valid_indicator
  );
  return;
}

void ExtendAndAssemble::SelectValidSuffixRight(
      const string& assembled_seq,
      const string& seed_seq,
      BoundType seed_range,
      unordered_map<unsigned int, bool>& valid_indicator
) {
#if DEBUG  
  cout << "   SelectValidSuffixRight:: selecting valid sequence..." << endl;
  cout << "   SelectValidSuffixRight:: assembled_seq  " << assembled_seq << endl;
  cout << "   SelectValidSuffixRight:: seed_seq " << seed_seq << endl;
#endif  
  if(seed_range.first < 0 || seed_range.second < 0 ||
      seed_range.first >= suffix_array_->getSize() ||
      seed_range.second >= suffix_array_->getSize())  {
    return;
  }
  // first check the length of all prefixs that ends at seed_alignment.seed_extend_seq
  typedef int SFAIndexType;
  typedef unsigned int SeqLengthType;
  unsigned int i;
  SeqLengthType max_prefix_len = 0;
  for(i = static_cast<unsigned int>(seed_range.first); i <= static_cast<unsigned int>(seed_range.second); ++ i) {
    SeqLengthType prefix_len = suffix_array_->getPos(static_cast<size_t>(i));
    max_prefix_len = prefix_len > max_prefix_len ? prefix_len : max_prefix_len;
  }
  // search the reverse suffix array with the reversed seed_seq to define upper and lower bounds
  string rev_seed_seq(seed_seq.rbegin(), seed_seq.rend());
  BoundType rev_seed_range = reverse_suffix_array_->searchWithLCPs(
      (SfaChar*) rev_seed_seq.c_str(), rev_seed_seq.length()
  );
  
  //while(rev_seed_range.first <= rev_seed_range.second &&
  //  (unsigned int) reverse_suffix_array_->getSuffixLength(static_cast<size_t>(rev_seed_range.first)) <= seed_seq.length()) {
  //  ++ rev_seed_range.first;
  //}
#if DEBUG  
  cout << "   SelectValidSuffixRight:: suffix array range:  " << seed_range.first << "\t" << seed_range.second << endl;
  cout << "   SelectValidSuffixRight:: prefix array range:  " << rev_seed_range.first << "\t" << rev_seed_range.second << endl;
  cout << "   SelectValidSuffixRight:: max prefix length: " << max_prefix_len << endl;
#endif
  // identify the position to start searching upward
  SFAIndexType position_search_upward = rev_seed_range.second;
  SeqLengthType prefix_array_search_len = seed_seq.length();
  
  // note that the assembled_sequence should contain the seed_seq as its suffix
  if(max_prefix_len >= assembled_seq.length() - seed_seq.length())  {
    //cout << "   SelectValidSuffixRight:: When prefix longer than or equal to assembled sequence" << endl;
    prefix_array_search_len = assembled_seq.length();
    string search_seq(assembled_seq.rbegin(), assembled_seq.rend());
    BoundType rev_range = reverse_suffix_array_->searchWithLCPs(
        (SfaChar*) search_seq.c_str(), search_seq.length()
    );
    if(rev_range.second >= rev_range.first)  {
      position_search_upward = rev_range.second;
    } else  {
      position_search_upward = -1;
    }
    //cout << "   SelectValidSuffixRight:: position determiend: " << position_search_upward << endl;  
  } else  {
    //cout << "   SelectValidSuffixRight:: When prefix shorter than assembled sequence" << endl;
    SeqLengthType left_bound = seed_seq.length();
    SeqLengthType right_bound = seed_seq.length() + max_prefix_len;
    right_bound = right_bound <= assembled_seq.length() ? right_bound : assembled_seq.length();
    prefix_array_search_len = right_bound;
    BoundType rev_range = rev_seed_range;
    BoundType search_range = rev_seed_range;
    //SeqLengthType sfa_lower_seq_len = reverse_suffix_array_->getSuffixLength(rev_seed_range.second);
    while(left_bound <= right_bound)  { 
      // a binary search style recursion to figure out the upward search position
      string search_seq(assembled_seq.rbegin(), assembled_seq.rbegin() + prefix_array_search_len);
      //cout << "   SelectValidSuffixRight:: sequence to search:  " << search_seq << endl; 
      rev_range = reverse_suffix_array_->searchWithLCPs(
          (SfaChar*) search_seq.c_str(), search_seq.length()
      );
      //sfa_lower_seq_len = reverse_suffix_array_->getSuffixLength(rev_range.second);
      //cout << "   SelectValidSuffixRight:: lower sequence length: " << sfa_lower_seq_len << endl;
      if(rev_range.second >= rev_range.first)  {
        // the pattern is found but with more than 1 hit, need to make the length longer
        position_search_upward = rev_range.second;
        //cout << "   SelectValidSuffixRight:: in binary search: pattern found, extending sequence" << endl;
        left_bound = prefix_array_search_len + 1;
        prefix_array_search_len = floor((left_bound + right_bound) / 2);  // record the last sucessful match
        search_range.first = rev_range.first + 1;
        // TODO: determine the suffix array search behavior when a search fails
      } else if(rev_range.second < rev_range.first) {
        // the pattern is not found, need to shorten the length to half
        //cout << "   SelectValidSuffixRight:: in binary search: pattern not found, shortening sequence" << endl;
        right_bound = prefix_array_search_len - 1;
        prefix_array_search_len = floor((left_bound + right_bound) / 2);
        search_range.second = search_range.second > rev_range.second ? rev_range.second : search_range.second;
      }
    } 
    //cout << "   SelectValidSuffixRight:: position determiend: " << position_search_upward << endl;  
  } 
  // intersect the read IDs in the defined boundaries
  BoundType rev_lookup_range = rev_seed_range;
  rev_lookup_range.second = position_search_upward;
  IntersectBoundRanges(
      suffix_array_, seed_range, 
      reverse_suffix_array_, rev_lookup_range, 
      seed_seq.length(), n_back_check_, assembled_seq.length(), 
      valid_indicator
  );
  return;
}

void ExtendAndAssemble::IntersectBoundRanges(
    GSA* ref_array,             // the reference array whose elements are being tested
    const BoundType ref_array_bound,  // the range in the reference array
    GSA* rev_ref_array,         // the reverse of the ref_array to check the prefix
    const BoundType rev_ref_array_bound,  // the range in the reverse reference array
    int seed_len,  // the length of the seed that was used to search both suffix arrays
    int n_back_check,
    int max_rev_length,
    unordered_map<unsigned int, bool>& valid_indicator
) {
  // first building a has of the read ids in the reference array
  //cout << "     IntersectBoundRages:: in IntersectBoundRanges, ref_bound: " << ref_array_bound.first << "\t" << ref_array_bound.second << endl;
  //cout << "     IntersectBoundRages:: in IntersectBoundRanges, srh_bound: " << rev_ref_array_bound.first << "\t" << rev_ref_array_bound.second << endl;
  //assert(ref_array_bound.first >= 0 && ref_array_bound.first < ref_array->getSize());
  //assert(ref_array_bound.second >= 0 && ref_array_bound.second < ref_array->getSize());
  //assert(rev_ref_array_bound.first >= 0 && rev_ref_array_bound.first < rev_ref_array->getSize());
  //assert(rev_ref_array_bound.second >= 0 && rev_ref_array_bound.second < rev_ref_array->getSize());
  
  if( ref_array_bound.first < 0 || ref_array_bound.second < 0 ||
      rev_ref_array_bound.first < 0 || rev_ref_array_bound.second < 0 ||
      ref_array_bound.first >= ref_array->getSize() || 
      ref_array_bound.second >= ref_array->getSize() ||
      rev_ref_array_bound.first >= rev_ref_array->getSize() || 
      rev_ref_array_bound.second >= rev_ref_array->getSize() ||
      ref_array_bound.first > ref_array_bound.second || 
      rev_ref_array_bound.first > rev_ref_array_bound.second)  {
    return;
  }
  
  if(n_back_check > 20 || n_back_check < 0)  {
    cout << "Warning: n-back-check size is " << n_back_check << ", which is an unusual setting and likely to cause incorrect results.\n";
  }
  unsigned int i;
  unsigned int bound_begin = static_cast<unsigned int>(ref_array_bound.first);
  unsigned int bound_end = static_cast<unsigned int>(ref_array_bound.second);
  unsigned int rev_bound_begin = static_cast<unsigned int>(rev_ref_array_bound.first);
  unsigned int rev_bound_end = static_cast<unsigned int>(rev_ref_array_bound.second);
  unordered_map<unsigned int, unsigned int> rid_hash; // the first field is the read ID,
      // the second field is the corresponding suffix array location in ref_array_bound
  for(i = bound_begin; i <= bound_end; ++ i) {
    valid_indicator[i] = false;
  }
  
  // constructing hash
  for(i = bound_begin; i <= bound_end; ++ i) {
    rid_hash[ref_array->getId(static_cast<size_t>(i))] = i;
  }
  unsigned int min_lcp = rev_ref_array->getSuffixLength(static_cast<size_t>(rev_ref_array_bound.second));
  for(i = rev_bound_end; i >= rev_bound_begin; -- i) {
    //cout << "      IntersectBoundRages:*" << rev_ref_array->getId(i) << "*checked" << endl;
    unsigned int len = rev_ref_array->getSuffixLength(static_cast<size_t>(i));
    //cout << "max_prefix_len:  " << max_prefix_len << endl;
    RIDType rev_read_ID = (RIDType) rev_ref_array->getId(i);
    auto it = rid_hash.find(rev_read_ID);
    if(len <= min_lcp && it != rid_hash.end())  {
      int rev_pos = rev_ref_array->getPos(static_cast<size_t>(i));
      int full_len = ref_array->getFullSequenceLength(static_cast<size_t>(it->second));
      int pos = ref_array->getPos(static_cast<size_t>(it->second));
      //cout << "     IntersectBoundRanges: positions:  " << (int) i << " " << (int) it->second << endl;
      //cout << "     IntersectBoundRanges: pos:  " << rev_pos << " " << pos << " " << seed_len << "  " << full_len << endl;
      if(rev_pos + pos + seed_len == full_len && pos <= max_rev_length)  {
        //cout << "      IntersectBoundRanges:*" << rev_read_ID << "*taken" << endl;
        valid_indicator[rid_hash[rev_read_ID]] = true;
      }
    }
    unsigned int current_lcp = rev_ref_array->getLcp(static_cast<size_t>(i));
    min_lcp = min_lcp < current_lcp ? min_lcp : current_lcp;
  }
  //cout << "     IntersectBoundRages:: list of valid reads: " << endl;
  //for(auto it = valid_indicator.begin(); it != valid_indicator.end(); ++ it) {
  //  if(it->second)  {
  //    cout << it->first << endl;
  //  }
  //}
  return;
} 

int ExtendAndAssemble::ComputeQueryRightBoundRight(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    SfaType extend_position
) {
  int suffix_len = suffix_array_->getSuffixLength(static_cast<size_t>(extend_position));
  int q_new_len = 
      source_alignment.current_assembled_seq.length() + 
      suffix_len - source_alignment.seed_extend_seq.length() + down_band_ - 1;
  //cout << "ComputeQueryRightBoundRight::  " << suffix_array_->getSuffix(extend_position) << endl; 
  //cout << "ComputeQueryRightBoundRight::  " << suffix_len << "  " << q_new_len << endl;
  int end = start_ + q_new_len - 1;
  return end;
}

int ExtendAndAssemble::ComputeQueryLeftBoundRight(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    SfaType extend_position
) {
  int end = ComputeQueryRightBoundRight(source_alignment, extend_position);
  int full_len = suffix_array_->getFullSequenceLength(static_cast<size_t>(extend_position));
  int begin = end - (full_len + right_band_ - 1 + down_band_ - 1) + 1; 
  return begin;
}
  
int ExtendAndAssemble::ComputeQueryLeftBoundLeft(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    SfaType extend_position
) {
  int suffix_len = reverse_suffix_array_->getSuffixLength(static_cast<size_t>(extend_position));
  int q_new_len = 
      source_alignment.current_assembled_seq.length() + 
      suffix_len - source_alignment.seed_extend_seq.length() + down_band_ - 1;
  int begin = start_ + seed_length_ - q_new_len;
  return begin;
}

int ExtendAndAssemble::ComputeQueryRightBoundLeft(
    const AlignmentInfoType<AlignScoreType>& source_alignment,
    SfaType extend_position
) {
  int begin = ComputeQueryLeftBoundLeft(source_alignment, extend_position);
  int full_len = reverse_suffix_array_->getFullSequenceLength(static_cast<size_t>(extend_position));
  int end = begin + (full_len + right_band_ - 1 + down_band_ - 1) - 1;
  return end;
}

bool ExtendAndAssemble::HasAssigned(const RIDType read_ID, const AlignmentPositionType position) {
  int ck_position = position;
  if(ck_position > (int) query_seq_->length() - 1)  {
    ck_position = query_seq_->length() - 1;
  }
  SLock r_lock(*mutex_assigned_reads_);
  auto it = assigned_reads_->find(read_ID);
	if(it != assigned_reads_->end())  {
	  auto it_h = clump_map_->find(ck_position); 
	  assert(it_h != clump_map_->end());
	  int idx = it_h->second;
	  auto it_map = it->second.find(idx);
	  if(it_map != it->second.end())  {
	    r_lock.unlock();
	    return true;  // ID found from hash table and is aligned nearby
	  }
	  r_lock.unlock();
	  return false; // ID found from hash table but not aligned nearby
	} 
	r_lock.unlock();
	return false;
}

bool ExtendAndAssemble::LookupAssignedVertex(const ReadAssignmentType& lookup_read) {
  int ck_position = lookup_read.alignment_position;
  if(ck_position > (int) query_seq_->length() - 1)  {
    ck_position = query_seq_->length() - 1;
  }
  auto it_h = clump_map_->find(ck_position); 
	assert(it_h != clump_map_->end());
	int idx = it_h->second;
  
  // searching the private read hash
	auto it = reads_to_record_->find(lookup_read.read_ID);
	if(it != reads_to_record_->end())  { 
    auto it_map = it->second.find(idx);
	  if(it_map != it->second.end())  {
      if(lookup_read.extend_direction == EXT_RIGHT && it_map->second.right_extended)  {
	      //end_edge = it_map->second.end_edge_right;
	      //end_vertex = it_map->second.end_vertex_right;
	      return true;
	    } else if(lookup_read.extend_direction == EXT_LEFT && it_map->second.left_extended) {
	      //end_edge = it_map->second.end_edge_left;
	      //end_vertex = it_map->second.end_vertex_left;
	      return true;
  	  } 
	  }
	}
	
  // searching the shared read hash
  SLock r_lock(*mutex_assigned_reads_);
  it = assigned_reads_->find(lookup_read.read_ID);
	if(it != assigned_reads_->end())  { 
    auto it_map = it->second.find(idx);
	  if(it_map != it->second.end())  {
      if(lookup_read.extend_direction == EXT_RIGHT && it_map->second.right_extended)  {
	      //end_edge = it_map->second.end_edge_right;
	      //end_vertex = it_map->second.end_vertex_right;
	      return true;
	    } else if(lookup_read.extend_direction == EXT_LEFT && it_map->second.left_extended) {
	      //end_edge = it_map->second.end_edge_left;
	      //end_vertex = it_map->second.end_vertex_left;
	      return true;
  	  } 
	  }
	} 
	r_lock.unlock();
	
	// ID not found from the hash table
  return false;
}

bool ExtendAndAssemble::UpdateAssignedReads(const RIDType& read_ID, const ReadTerminateType& mapped_info)  {
  //cout << "Update begin" << endl;
  bool updated = false;
  int ck_position = mapped_info.alignment_position;
  if(ck_position > (int) query_seq_->length() - 1)  {
    ck_position = query_seq_->length() - 1;
  }
  auto it_h = clump_map_->find(ck_position); 
	assert(it_h != clump_map_->end());
	int idx = it_h->second;
  //ULock w_lock(*mutex_assigned_reads_);
  //int idx = (int) mapped_info.alignment_position / clump_range_;
  auto it = reads_to_record_->find(read_ID);
	if(it != reads_to_record_->end())  {  
	  //cout << "reads found" << endl;
	  assert(mapped_info.right_extended || mapped_info.left_extended);	  
	  auto it_map = it->second.find(idx);
	  if(it_map != it->second.end())  {    
	    if(mapped_info.right_extended)  {
	      if(!it_map->second.right_extended || 
	          (it_map->second.right_extended && mapped_info.assembled_score_right > it_map->second.assembled_score_right))  {
	        it_map->second.right_extended = true;
	        it_map->second.assembled_score_right = mapped_info.assembled_score_right;
	        //it_map->second.end_vertex_right = mapped_info.end_vertex_right;
	        //it_map->second.end_edge_right = mapped_info.end_edge_right;
	        updated = true;
	      }
	    } else if(mapped_info.left_extended) {
	      if(!it_map->second.left_extended || 
	          (it_map->second.left_extended && mapped_info.assembled_score_left > it_map->second.assembled_score_left))  {
	        it_map->second.left_extended = true;
	        it_map->second.assembled_score_left = mapped_info.assembled_score_left;
	        //it_map->second.end_vertex_left = mapped_info.end_vertex_left;
	        //it_map->second.end_edge_left = mapped_info.end_edge_left;
	        updated = true;
	      }
	    } 
	  }
	} 
	// ID not found from the hash table
	//cout << "reads not found" << endl;
	//cout << "aligned position:  " << mapped_info.alignment_position << endl;
	//cout << "query_len: " << query_seq_->length() << endl;
  //cout << "size to initialize:  " << ((unsigned int) (query_seq_->length() + down_band_) / clump_range_) + 1 << endl;
  //cout << "idx: " << idx << endl;
  (*reads_to_record_)[read_ID].insert({idx, mapped_info});
	
	//cout << "before adding the vector" << endl;
	
	updated = true;
	//w_lock.unlock();
	//cout << "Update end" << endl;
	return updated;
}


