#include "assembly_graph.h"

using namespace std;

AssemblyGraph::AssemblyGraph()  {
  return;
}

AssemblyGraph::~AssemblyGraph() {
  return;
}

BoostVertex AssemblyGraph::AddVertex(const VertexProperty& data) {
  //ULock w_lock(mutex_);
  BoostVertex v = boost::add_vertex(data, graph_);
  if(boost::num_vertices(graph_) <= 1)  {
    source_ = v;
  }
  //w_lock.unlock();
  return v;
}

BoostEdge AssemblyGraph::AddEdge(const BoostVertex& u, const BoostVertex& v, const EdgeProperty& data)  {
  BoostEdge e;
  bool is_success = false;
  //ULock w_lock(mutex_);
  tie(e, is_success) = boost::add_edge(u, v, data, graph_);
  //w_lock.unlock();
  if(!is_success)  {
    std::cerr << "AssemblyGraph::AddEdge:Error in adding edge (through calling Boost Graph Library)." << std::endl;
  }
  return e;
}

VertexProperty AssemblyGraph::AccessVertex(const BoostVertex& v) {
  //SLock r_lock(*(graph_[v].mutex_vertex));
  VertexProperty data = graph_[v];
  //r_lock.unlock();
  return data;
}

EdgeProperty AssemblyGraph::AccessEdge(const BoostEdge& e) {
  //SLock r_lock(*(graph_[e].mutex_edge));
  EdgeProperty data = graph_[e];
  //r_lock.unlock();
  return data;
}


void AssemblyGraph::UpdateVertex(const BoostVertex& v, const VertexProperty& n)  {
  //ULock w_lock(*(graph_[v].mutex_vertex));
  graph_[v] = n;
  //w_lock.unlock();
  return;
}

void AssemblyGraph::UpdateEdge(const BoostEdge& e, const EdgeProperty& n)  {
  //SLock w_lock(*(graph_[e].mutex_edge));
  graph_[e] = n;
  //w_lock.unlock();
  return;
}

int AssemblyGraph::GetNumVertices(void) {
  //SLock r_lock(mutex_);
  int n = boost::num_vertices(graph_);
  //r_lock.unlock();
  return n;
}

int AssemblyGraph::GetNumEdges(void)  {
  //SLock r_lock(mutex_);
  int n = boost::num_edges(graph_);
  //r_lock.unlock();
  return n;
} 

std::pair<BoostVertexIter, BoostVertexIter> AssemblyGraph::GetVertices(void) {
  //SLock r_lock(mutex_);
  std::pair<BoostVertexIter, BoostVertexIter> n;
  n = boost::vertices(graph_);
  //r_lock.unlock();
  return n;
}

std::pair<BoostEdgeIter, BoostEdgeIter> AssemblyGraph::GetEdges(void) {
  //SLock r_lock(mutex_);
  std::pair<BoostEdgeIter, BoostEdgeIter> n;
  n = boost::edges(graph_);
  //r_lock.unlock();
  return n;
}

void AssemblyGraph::GetInEdges(const BoostVertex& v, std::list<BoostEdge>& ie_list) {
  //SLock r_lock(*(graph_[v].mutex_vertex));
  if(boost::in_degree(v, graph_) <= 0)  {
    //r_lock.unlock();
    return;
  }
  auto edge_range = boost::in_edges(v, graph_);
  //r_lock.unlock();
  while(edge_range.first != edge_range.second)  {
    ie_list.push_back(*edge_range.first);
    ++ edge_range.first;
  }
  return;
}

void AssemblyGraph::GetOutEdges(const BoostVertex& v, std::list<BoostEdge>& oe_list)  {
  //SLock r_lock(*(graph_[v].mutex_vertex));
  if(boost::out_degree(v, graph_) <= 0)  {
    //r_lock.unlock();
    return;
  }
  auto edge_range = boost::out_edges(v, graph_);
  //r_lock.unlock();
  while(edge_range.first != edge_range.second)  {
    oe_list.push_back(*edge_range.first);
    ++ edge_range.first;
  }
  return;
}

int AssemblyGraph::GetInDegree(const BoostVertex& v) {
  //SLock r_lock(*(graph_[v].mutex_vertex));
  int n = boost::in_degree(v, graph_);  
  //r_lock.unlock();
  return n;
}

int AssemblyGraph::GetOutDegree(const BoostVertex& v)  {
  //SLock r_lock(*(graph_[v].mutex_vertex));
  int n = boost::out_degree(v, graph_);
  //r_lock.unlock();
  return n;
}

int AssemblyGraph::GetDegree(const BoostVertex& v) {
  //SLock r_lock(*(graph_[v].mutex_vertex));
  int n = boost::degree(v, graph_);
  //r_lock.unlock();
  return n;
}

BoostVertex AssemblyGraph::GetSource(const BoostEdge& e)  {
  //SLock r_lock(*(graph_[e].mutex_edge));
  BoostVertex n = boost::source(e, graph_);
  //r_lock.unlock();
  return n;
}

BoostVertex AssemblyGraph::GetTarget(const BoostEdge& e) {
  //SLock r_lock(*(graph_[e].mutex_edge));
  BoostVertex n = boost::target(e, graph_);
  //r_lock.unlock();
  return n;
}

void AssemblyGraph::RemoveEdge(const BoostEdge& e) {
  //ULock w_lock(mutex_);
  boost::remove_edge(e, graph_);
  //w_lock.unlock();
  return;
}

void AssemblyGraph::RemoveVertex(const BoostVertex& v) {
  //ULock w_lock(mutex_);
  boost::remove_vertex(v, graph_);
  //w_lock.unlock();
  return;
}

BoostEdge AssemblyGraph::GetEdge(const BoostVertex& u, const BoostVertex& v) {
  //cout << "GetEdge called" << endl;
  BoostEdge e;
  bool success;
  //SLock r_lock_u(*(graph_[u].mutex_vertex));
  //SLock r_lock_v(*(graph_[v].mutex_vertex));
  //SLock r_lock_g(mutex_);
  tie(e, success) = boost::edge(u, v, graph_);
  //r_lock_g.unlock();
  //r_lock_v.unlock();
  //r_lock_u.unlock();
  if(!success)  {
    cerr << "AssemblyGraph::GetEdge: error in finding edge between two vertices: no such edge exists." << endl;
    exit(1);
  }
  return e;
}

void AssemblyGraph::UpdateAbove(
    BoostVertex& source, BoostVertex& target, 
    AlignScoreType path_score_max
) {
  if(graph_[source].max_score_below < path_score_max)  {
    graph_[source].max_score_below = path_score_max;
    //cout << " updated max_vertex_above" << endl;
    graph_[source].max_vertex_below = target;
  }
  return;
}

void AssemblyGraph::UpdateBelow(
    BoostVertex& source, BoostVertex& target, 
    AlignScoreType path_score_max
) {
  if(graph_[target].max_score_above < path_score_max)  {
    graph_[target].max_score_above = path_score_max;
    //cout << " updated max_vertex_below" << endl;
    graph_[target].max_vertex_above = source;
  }
  return;
}

/*
void AssemblyGraph::TrackDownBestPath(
    const BoostVertex& pivot,
    EdgeProperty& e_current,
    EdgeProperty& e_complete
) {
  // tracks down the best path that follows the pivot vertex
  EdgeProperty e_exist, e_extend, e_phase;
  e_exist = e_current;
  BoostVertex vertex_phase = pivot;
  while(GetOutDegree(vertex_phase) > 0) {
    BoostEdge max_edge = GetEdge(vertex_phase, graph_[vertex_phase].max_vertex_below); 
    //cout << "score_below: " << graph_[vertex_phase].max_score_below << "  " << graph_[vertex_phase].min_score_below << endl;
    e_extend = graph_[max_edge];
    ExtendEdge(e_exist, e_extend, e_phase);
    e_exist = e_phase;
    vertex_phase = graph_[vertex_phase].max_vertex_below;
  }
  e_complete = e_phase;
  return;
}
*/

/*
void AssemblyGraph::ExtendEdge(
    EdgeProperty& edge_current, 
    EdgeProperty& edge_extend, 
    EdgeProperty& edge_next
) {
  
  if(edge_extend.extend_direction == EXT_RIGHT)  {
    edge_next.outgoing_seq = edge_current.outgoing_seq + edge_extend.outgoing_seq;
  } else  {
    edge_next.outgoing_seq = ReverseShortSequence(edge_extend.outgoing_seq) + edge_current.outgoing_seq;
  }
  edge_next.query_end = edge_extend.query_end;
  edge_next.extend_direction = edge_extend.extend_direction;
  edge_next.current_score = edge_extend.current_score;
  edge_next.increased_score = edge_current.increased_score + edge_extend.increased_score;
  edge_next.included_reads = edge_current.included_reads;
  for(auto it_nr = edge_extend.included_reads.begin(); it_nr != edge_extend.included_reads.end(); ++ it_nr) {
    it_nr->offset_to_start = edge_current.outgoing_seq.length() - (unsigned int) it_nr->position - it_nr->len_search_seed;
    edge_next.included_reads.push_back(*it_nr);     
  }
  return;
}
*/

void AssemblyGraph::MergeReadsWithParent(BoostVertex& v_parent, BoostVertex& v_child)  {
  // first add the accumulative reads in the child vertex to the parent
  BoostEdge e_between = GetEdge(v_parent, v_child);
  AlignScoreType inc_score = graph_[e_between].increased_score;
  for(auto it = graph_[v_child].reads_in_children.begin(); it != graph_[v_child].reads_in_children.end(); ++ it) {
    it->second.max_score += inc_score;
    it->second.min_score += inc_score;
    auto parent_read_iter = graph_[v_parent].reads_in_children.find(it->first);
    if(parent_read_iter == graph_[v_parent].reads_in_children.end())  {
      (graph_[v_parent].reads_in_children)[it->first] = it->second;
    } else  {
      if(it->second.max_score > parent_read_iter->second.max_score)  {
        parent_read_iter->second.max_score = it->second.max_score;
        parent_read_iter->second.max_vertex_source = it->second.max_vertex_source;
        parent_read_iter->second.max_vertex_target = it->second.max_vertex_target;
      }
      if(it->second.min_score < parent_read_iter->second.min_score)  {
        parent_read_iter->second.min_score = it->second.min_score;
        parent_read_iter->second.min_vertex_source = it->second.min_vertex_source;
        parent_read_iter->second.min_vertex_target = it->second.min_vertex_target;
      }
    }
  }
  return;
}

void AssemblyGraph::AddEdgeReadsToParent(
    BoostVertex& v_parent, BoostVertex& v_child,
    AlignScoreType max_score
) {
  BoostEdge e_between = GetEdge(v_parent, v_child);
  for(auto it_read = graph_[e_between].included_reads.begin(); it_read != graph_[e_between].included_reads.end(); ++ it_read) {
    // add single read to the read has in parent
    auto parent_read_iter = graph_[v_parent].reads_in_children.find(it_read->read_ID);
    if(parent_read_iter == graph_[v_parent].reads_in_children.end())  {
      ReadVertexType path_info;
      path_info.max_score = max_score;
      path_info.max_vertex_source = v_parent;
      path_info.max_vertex_target = v_child;
      (graph_[v_parent].reads_in_children)[it_read->read_ID] = path_info;
    } else  {
      if(max_score > parent_read_iter->second.max_score)  {
        parent_read_iter->second.max_score = max_score;
        parent_read_iter->second.max_vertex_source = v_parent;
        parent_read_iter->second.max_vertex_target = v_child;
      }
    }
  }
  return;
}

void AssemblyGraph::TraverseGraphUnder(
    const BoostVertex& source, 
    std::unordered_map<RIDType, ReadVertexType>& recorded_reads
)  {
  //cout << "NEW CALL OF TRAVERSEGRAPHUNDER *********************************" << endl;
  stack<BoostVertex> unvisited;
  stack<AlignScoreType> accumulate_score;
  // initialize the traversal (DFS)
  unvisited.push(source);
  graph_[source].max_score_above = 0;
  accumulate_score.push(0);
  while(!unvisited.empty()) {
    //cout << "Visiting a node: ###############################" << endl;
    BoostVertex vertex_phase = unvisited.top();
    AlignScoreType score_phase = accumulate_score.top();
    //cout << "score_below: " << graph_[vertex_phase].max_score_below << "  " << graph_[vertex_phase].min_score_below << endl;
    if(!graph_[vertex_phase].touched && !graph_[vertex_phase].visited)  {
      //cout << "reach a new node" << endl;
      // the vertex is neither touched or visited, update score for its children
      graph_[vertex_phase].touched = true;
      list<BoostEdge> out_edges;
      GetOutEdges(vertex_phase, out_edges);
      for(auto it_oe = out_edges.begin(); it_oe != out_edges.end(); ++ it_oe) {
        BoostVertex vertex_next = GetTarget(*it_oe);
        // update max and min above score
        AlignScoreType path_score_max = graph_[*it_oe].increased_score + graph_[vertex_phase].max_score_above;
        // update the score of the vertex
        UpdateBelow(vertex_phase, vertex_next, path_score_max);
        if(graph_[vertex_next].visited)  {
          // don't do anything if the next vertex has already been visited
          //cout << "visited node observed" << endl;
          continue;
        }
        // otherwise the vertex is untouched, push the vertices to the stack, also push the new edge property
        //cout << "new node pushed" << endl;
        unvisited.push(vertex_next);
        accumulate_score.push(score_phase + graph_[*it_oe].increased_score);
      }
    } else if(graph_[vertex_phase].touched && !graph_[vertex_phase].visited)  {
      // the vertex is touched but not visited, update score for its parents
      //cout << "revisiting a node" << endl;
      graph_[vertex_phase].visited = true;
      
      // if the node is a leaf node
      if(GetOutDegree(vertex_phase) == 0)  {
        //cout << "Reaching a leaf node" << endl;
        graph_[vertex_phase].max_score_below = 0;
        graph_[vertex_phase].max_vertex_below = NULL;
      }
      // if the node is the source
      if(unvisited.size() == 1)  {
        // record all the recruited reads
        recorded_reads = graph_[source].reads_in_children;
        unvisited.pop();
        return;
      }
      
      
      list<BoostEdge> in_edges;
      GetInEdges(vertex_phase, in_edges);
      //cout << "In degree: " << GetInDegree(vertex_phase) << endl;
      //cout << "Out degree:  " << GetOutDegree(vertex_phase) << endl;
      //cout << "size in edges: " << in_edges.size() << endl;
      for(auto it_ie = in_edges.begin(); it_ie != in_edges.end(); ++ it_ie) {
        BoostVertex vertex_parent = GetSource(*it_ie);

        AlignScoreType path_score_max = graph_[vertex_phase].max_score_below + graph_[*it_ie].increased_score; 
        
        UpdateAbove(vertex_parent, vertex_phase, path_score_max);
        MergeReadsWithParent(vertex_parent, vertex_phase);
        AddEdgeReadsToParent(vertex_parent, vertex_phase, path_score_max);
      }
      // delete all reads in the child node
      graph_[vertex_phase].reads_in_children.clear();
      unvisited.pop();
    } else if(graph_[vertex_phase].visited) {
      //cout << "visiting a resolved node" << endl;
      unvisited.pop();
    }
  }
  
  return;
}

bool _cmp_assembled_path(const AssembledPathType& item_a, const AssembledPathType& item_b)  {
  if(item_a.assembled_sequence.length() > item_b.assembled_sequence.length())  {
    return true;
  }
  return false;
}

void AssemblyGraph::SpellMinimumPaths(
    GSA& suffix_array, GSA& reverse_suffix_array,
    const BoostVertex& source, std::set<RIDType> reads_to_visit,
    std::list<AssembledPathType>& assembled_seqs,
    AssembledPathType& best_seq
) {

  set<RIDType> included_reads;
  list<AssembledPathType> assembled_holder;
  
  stack<BoostVertex> unvisited;
  stack<AssembledPathType> accumulate_seqs;
  stack<std::set<RIDType> > visited_in_path;
  stack<AlignScoreType> alignment_scores;
  // initialize the traversal (DFS)
  unvisited.push(source);
  AssembledPathType init_seq;
  init_seq.assembled_sequence = "";
  accumulate_seqs.push(init_seq);
  set<RIDType> init_set;
  visited_in_path.push(init_set);
  alignment_scores.push(0);
  // visit all vertices
  while(!unvisited.empty()) {
    // get the top node and edge
    BoostVertex v = unvisited.top();
    AssembledPathType acc_seq_phase = accumulate_seqs.top();
    set<RIDType> visited_phase = visited_in_path.top();
    AlignScoreType align_score_phase = alignment_scores.top();
    // remove the vertex and edge from the stack
    unvisited.pop();
    accumulate_seqs.pop();
    visited_in_path.pop();
    alignment_scores.pop();
    // for all the vertex's children
    list<BoostEdge> out_edges;
    GetOutEdges(v, out_edges);
    if(out_edges.size() == 0)  {
      assembled_holder.push_back(acc_seq_phase);
      continue;
    }

    for(auto it = out_edges.begin(); it != out_edges.end(); ++ it) {
      // extend the sequence and recruit the reads
      EdgeProperty e_property = AccessEdge(*it);
      AssembledPathType next_seq_phase = acc_seq_phase;
      set<RIDType> next_visited_reads = visited_phase;
      for(auto it_nr = e_property.included_reads.begin(); it_nr != e_property.included_reads.end(); ++ it_nr) {
        int end_position = next_seq_phase.assembled_sequence.length() - (unsigned int) it_nr->position - it_nr->len_search_seed;
        next_seq_phase.mapped_locations.push_back({it_nr->read_ID, end_position});
        if(reads_to_visit.find(it_nr->read_ID) != reads_to_visit.end())  {
          next_visited_reads.insert(it_nr->read_ID);
        }
      }
      //cout << "******************************" << endl;
      //next_seq_phase.assembled_sequence += e_property.outgoing_seq;
      //cout << e_property.outgoing_seq << endl;
      if(e_property.extend_direction == EXT_RIGHT)  {
        //cout << "right" << endl;
        //cout << suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos) << endl;
        next_seq_phase.assembled_sequence += 
            string(suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos));
      } else  {
        //cout << "left" << endl;
        //cout << reverse_suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos) << endl;
        next_seq_phase.assembled_sequence +=
            string(reverse_suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos));
      }
      
      accumulate_seqs.push(next_seq_phase);
      // record the next vertex to be visited
      unvisited.push(GetTarget(*it));
      // record vertices that have been visited
      visited_in_path.push(next_visited_reads);
      // computes the alignment score;
      alignment_scores.push(align_score_phase + e_property.increased_score);
    }
  }
  
  // finalize the path selection
  assembled_holder.sort(_cmp_assembled_path);
  for(auto it = assembled_holder.begin(); it != assembled_holder.end(); ++ it) {
    //if((double) included_reads.size() / reads_to_visit.size() >= 1.0 - 0.01)  {
    //  break;
    //}
    int num_covered = 0, num_critical = 0;
    set<RIDType> read_not_covered;
    for(auto it_r = it->mapped_locations.begin(); it_r != it->mapped_locations.end(); ++ it_r) {
      if(reads_to_visit.find(it_r->first) != reads_to_visit.end())  {
        ++ num_critical;
        if(included_reads.find(it_r->first) != included_reads.end())  {
          ++ num_covered;
        } else  {
          read_not_covered.insert(it_r->first);
        }
      }
    }
    //if(num_critical > 0 && (double) read_not_covered.size() / reads_to_visit.size() >= 0.01 
    //    && (double) num_covered / num_critical <= 1.0 - 0.01)  {
    assembled_seqs.push_back(*it);
    if(best_seq.assembled_sequence.length() <= 0)  {
      best_seq = *it;
    }
    for(auto it_id = read_not_covered.begin(); it_id != read_not_covered.end(); ++ it_id) {
      included_reads.insert(*it_id);
    }
    //}
  }

  return;
}

void AssemblyGraph::SpellAllPaths(
    GSA& suffix_array, GSA& reverse_suffix_array,
    const BoostVertex& source, std::list<AssembledPathType>& assembled_seqs
) {

  set<RIDType> included_reads;
  list<AssembledPathType> assembled_holder;
  
  stack<BoostVertex> unvisited;
  stack<AssembledPathType> accumulate_seqs;
  stack<std::set<RIDType> > visited_in_path;
  stack<AlignScoreType> alignment_scores;
  // initialize the traversal (DFS)
  unvisited.push(source);
  AssembledPathType init_seq;
  init_seq.assembled_sequence = "";
  init_seq.align_score = 0;
  accumulate_seqs.push(init_seq);
  set<RIDType> init_set;
  visited_in_path.push(init_set);
  alignment_scores.push(0);
  // visit all vertices
  while(!unvisited.empty()) {
    // get the top node and edge
    BoostVertex v = unvisited.top();
    AssembledPathType acc_seq_phase = accumulate_seqs.top();
    set<RIDType> visited_phase = visited_in_path.top();
    AlignScoreType align_score_phase = alignment_scores.top();
    // remove the vertex and edge from the stack
    unvisited.pop();
    accumulate_seqs.pop();
    visited_in_path.pop();
    alignment_scores.pop();
    // for all the vertex's children
    list<BoostEdge> out_edges;
    GetOutEdges(v, out_edges);
    if(out_edges.size() == 0)  {
      assembled_seqs.push_back(acc_seq_phase);
      continue;
    }

    for(auto it = out_edges.begin(); it != out_edges.end(); ++ it) {
      // extend the sequence and recruit the reads
      EdgeProperty e_property = AccessEdge(*it);
      AssembledPathType next_seq_phase = acc_seq_phase;
      set<RIDType> next_visited_reads = visited_phase;
      for(auto it_nr = e_property.included_reads.begin(); it_nr != e_property.included_reads.end(); ++ it_nr) {
        int end_position = next_seq_phase.assembled_sequence.length() - (unsigned int) it_nr->position - it_nr->len_search_seed;
        next_seq_phase.mapped_locations.push_back({it_nr->read_ID, end_position});
      }
      //cout << "******************************" << endl;
      //next_seq_phase.assembled_sequence += e_property.outgoing_seq;
      //cout << e_property.outgoing_seq << endl;
      if(e_property.extend_direction == EXT_RIGHT)  {
        //cout << "right" << endl;
        //cout << suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos) << endl;
        next_seq_phase.assembled_sequence += 
            string(suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos));
      } else  {
        //cout << "left" << endl;
        //cout << reverse_suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos) << endl;
        next_seq_phase.assembled_sequence +=
            string(reverse_suffix_array.getSuffix_explicit(e_property.outgoing_seq_ID, e_property.outgoing_seq_Pos));
      }
      next_seq_phase.align_score = align_score_phase + e_property.increased_score;
      accumulate_seqs.push(next_seq_phase);
      // record the next vertex to be visited
      unvisited.push(GetTarget(*it));
      // record vertices that have been visited
      visited_in_path.push(next_visited_reads);
      // computes the alignment score;
      alignment_scores.push(align_score_phase + e_property.increased_score);
    }
  }
 

  return;
}

