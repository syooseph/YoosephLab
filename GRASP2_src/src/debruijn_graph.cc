#include "../include/debruijn_graph.h"

using namespace std;

DeBruijnGraph::~DeBruijnGraph() {
  if(initialized_)  {
    delete p_graph_;
    delete p_info_map_;
  }
}

bool DeBruijnGraph::PassFreqFilter(
    const std::string &seq, 
    const FrequencyTable<KmerUnitType> &unit_freq, const int freq_cutoff
) {
  int ml = encoder_.GetMerLen();
  if(seq.length() < ml) return false;
  KmerUnitType kunit; 
  for(int i = 0; i < seq.length() - encoder_.mer_len_; ++ i) {
    if(i == 0) kunit = encoder_.Encode(&(seq.c_str()[0]));
    else kunit = encoder_.RightExt(kunit, seq[i + ml - 1]);
    if(unit_freq.GetFreq(kunit) < freq_cutoff)  {
      //cout << seq << "  " << i << " " << kunit << endl;
      return false;
    }
  }
  return true; 
}


void DeBruijnGraph::Construct(
    const int num_seq, char ** const seq, std::vector<bool> &read_mark,
    const FrequencyTable<KmerUnitType> &unit_freq, const int freq_cutoff
) {
  // check for parameters
  if(freq_cutoff < 0)  {
    cout << "DeBruijnGraph::Warning: minimum frequency cutoff set to negative value " << freq_cutoff << "! Reset to 0." << endl;
  }
  // initilize the object
  p_graph_ = new BoostGraph;
  p_info_map_ = new std::unordered_map<std::string, BoostVertex>;
  initialized_ = true;
  // construct the graph  
  int i, j;
  int kmer_id = 0;
  for(i = 0; i < num_seq; ++ i) {
    if(read_mark[i]) continue;
    string s = seq[i];
    if(s.length() < kmer_size_) continue;
    BoostVertex prev_vertex;
    bool is_prev_set = false;
    for(int j = 0; j < s.length() - kmer_size_ + 1; ++ j) {
      string mer = s.substr(j, kmer_size_);
      /*****/
      //cout << mer << endl;
      /*****/
      if(PassFreqFilter(mer, unit_freq, freq_cutoff))  {
        /*****/
        //cout << "Passed frequency filter" << endl;
        /*****/
      
        auto it = p_info_map_->find(mer);
        // add the vertex to the hash table and the graph
        if(it == p_info_map_->end())  {
          (*p_info_map_)[mer] = add_vertex(VertexType(mer), *p_graph_);
        } else  {
          (*p_graph_)[it->second].multiplicity_ ++;
        }
        // add the edge to the graph
        if(is_prev_set)  {
          std::pair<BoostEdge, bool> e_search = 
            boost::edge(prev_vertex, (*p_info_map_)[mer], *p_graph_);
          if(!e_search.second) add_edge(prev_vertex, (*p_info_map_)[mer], *p_graph_);
          else (*p_graph_)[e_search.first].multiplicity_ ++;
        }
        prev_vertex = (*p_info_map_)[mer];  is_prev_set = true;
      } else is_prev_set = false; // we need two consequtive good kmers to add an edge
    }
  }
  //cout << "Recorded " << p_info_map_->size() << " high-frequency kmers." << endl;
  return;
}


void DeBruijnGraph::RefineGraph(const int freq_cutoff) {
  list<BoostEdge> to_delete;
  for(auto it_v = boost::vertices(*p_graph_).first; 
    it_v != boost::vertices(*p_graph_).second; ++ it_v
  ) {
    if((*p_graph_)[*it_v].multiplicity_ < freq_cutoff)  {
      for(auto it_e = in_edges(*it_v, *p_graph_).first; 
        it_e != in_edges(*it_v, *p_graph_).second; ++ it_e
      ) {
        to_delete.push_back(*it_e);
      }
      for(auto it_e = out_edges(*it_v, *p_graph_).first; 
        it_e != out_edges(*it_v, *p_graph_).second; ++ it_e
      ) {
        to_delete.push_back(*it_e);
      }
    }
  }
  RemoveEdges(to_delete);
  RemoveOrphantVertices();
  return;
}

inline void DeBruijnGraph::ReplaceColor(std::list<BoostVertex> &vertex_list, const int color) {
  for(auto it = vertex_list.begin(); it != vertex_list.end(); ++ it) 
    (*p_graph_)[*it].color_ = color;
}

void DeBruijnGraph::FormMinSpanningTree(const int min_cov) {  
  list<BoostEdge> to_delete;
  unordered_map<int, list<BoostEdge> > multiplicity_map;
  for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++ it_e)
    multiplicity_map[(*p_graph_)[*it_e].multiplicity_].push_back(*it_e);
  // sort the edges based on multiplicity
  list<int> multiplicity_keys;
  for(auto it = multiplicity_map.begin(); it != multiplicity_map.end(); ++ it)
    multiplicity_keys.push_back(it->first);
  multiplicity_keys.sort();
  // visit edges one-by-one
  int color_id = 1;
  unordered_map<int, list<BoostVertex> > color_map;
  for(auto it = multiplicity_keys.rbegin(); it != multiplicity_keys.rend(); ++ it) {
    // DEBUG: cout << "multiplicity:  " << *it << " size: " << multiplicity_map[*it].size() << endl;
    // if the vertex has color of "0" that the vertex is not colored
    for(auto it_e = multiplicity_map[*it].begin(); it_e != multiplicity_map[*it].end(); ++ it_e) {
      if(*it < min_cov) {
        to_delete.push_back(*it_e); 
        continue;
      } 
      BoostVertex u = source(*it_e, *p_graph_), v = target(*it_e, *p_graph_);
      int uc = (*p_graph_)[u].color_, vc = (*p_graph_)[v].color_;
      if(uc == vc && uc == 0)  {
        // the two vertices are not connected before
        // create a new color, and connect the two vertices
        // DEBUG: cout << "Create New vertex" << endl;
        (*p_graph_)[u].color_ = (*p_graph_)[v].color_ = color_id;
        color_map[color_id].push_back(u);
        color_map[color_id].push_back(v);
        ++ color_id;
      } else if(uc == vc && uc != 0) {
        // the two vertices belong to trees that have previously connected
        // remove the edge
        // DEBUG: cout << "Redundant edge" << endl;
        to_delete.push_back(*it_e);
      } else if(uc == 0 && vc != 0) {
        // v belong to a sub-tree and u was not connected
        // connect u to v
        // DEBUG: cout << "Absort vertex" << endl;
        (*p_graph_)[u].color_ = vc;
        color_map[vc].push_back(u);
      } else if(uc != 0 && vc == 0) {
        // symmetirc to the previous case
        // DEBUG: cout << "Absort vertex" << endl;
        (*p_graph_)[v].color_ = uc;
        color_map[uc].push_back(v);
      } else  {
        // DEBUG: cout << "Merging trees vertex" << endl;
        // both vertices belong to different sub-trees; merge the two sub-trees
        if(color_map[vc].size() < color_map[uc].size()) {
          ReplaceColor(color_map[vc], uc);
          for(auto it = color_map[vc].begin(); it != color_map[vc].end(); ++ it)
            color_map[uc].push_back(*it);
          color_map.erase(vc);
        } else  {
          ReplaceColor(color_map[uc], vc);
          for(auto it = color_map[uc].begin(); it != color_map[uc].end(); ++ it)
            color_map[vc].push_back(*it);
          color_map.erase(uc);
        }
      }
    }
  }
  RemoveEdges(to_delete);
  RemoveOrphantVertices();
  return;
}

void DeBruijnGraph::RemoveTipsHead(const int length)  {
  list<BoostEdge> to_delete;
  for(auto it_v = vertices(*p_graph_).first; it_v != vertices(*p_graph_).second; ++ it_v) {
    // for edges that are heads
    if(in_degree(*it_v, *p_graph_) == 0)  {
      list<BoostEdge> traversed;
      BoostVertex current = *it_v;
      int depth = 0;
      while(depth < length && 
          in_degree(current, *p_graph_) <= 1 && out_degree(current, *p_graph_) == 1
      ) {
        BoostEdge out_e = *(out_edges(current, *p_graph_).first);
        traversed.push_back(out_e);
        current = target(out_e, *p_graph_);
        ++ depth;
      }
      if(depth < length)  {
        // delete all traversed edges
        for(auto it = traversed.begin(); it != traversed.end(); ++ it)
          to_delete.push_back(*it);
      }
    }
  }
  RemoveEdges(to_delete);
  RemoveOrphantVertices();
  return;
}

void DeBruijnGraph::RemoveTipsTail(const int length)  {
  list<BoostEdge> to_delete;
  for(auto it_v = vertices(*p_graph_).first; it_v != vertices(*p_graph_).second; ++ it_v) {
    // for edges that are tails
    if(out_degree(*it_v, *p_graph_) == 0)  {
      list<BoostEdge> traversed;
      BoostVertex current = *it_v;
      int depth = 0;
      while(depth < length && 
          out_degree(current, *p_graph_) <= 1 && in_degree(current, *p_graph_) == 1
      ) {
        BoostEdge out_e = *(in_edges(current, *p_graph_).first);
        traversed.push_back(out_e);
        current = source(out_e, *p_graph_);
        ++ depth;
      }
      if(depth < length)  {
        // delete all traversed edges
        for(auto it = traversed.begin(); it != traversed.end(); ++ it)
          to_delete.push_back(*it);
      }
    }
  }
  RemoveEdges(to_delete);
  RemoveOrphantVertices();
  return;
}

void DeBruijnGraph::BreakBadCovBranch(const float fraction) {
  list<BoostEdge> to_delete;
  for(auto it = vertices(*p_graph_).first; it != vertices(*p_graph_).second; ++ it) {
    // has multiple out degrees
    for(auto it_e = out_edges(*it, *p_graph_).first; 
        it_e != out_edges(*it, *p_graph_).second; ++ it_e
    ) {
      if((float) (*p_graph_)[*it_e].multiplicity_ < 
          (float) fraction * (*p_graph_)[*it].multiplicity_
      ) {
#ifdef DEBUG
        cout << "To delete out_edge:  " << (*p_graph_)[*it_e].multiplicity_ << "  " << (*p_graph_)[*it].multiplicity_ << endl; 
#endif
        to_delete.push_back(*it_e);     
      }
    }
    // has multiple in degrees
    for(auto it_e = in_edges(*it, *p_graph_).first; 
        it_e != in_edges(*it, *p_graph_).second; ++ it_e
    ) {
      if((float) (*p_graph_)[*it_e].multiplicity_ < 
          (float) fraction * (*p_graph_)[*it].multiplicity_
      ) {
#ifdef DEBUG
        cout << "To delete in_edge:  " << (*p_graph_)[*it_e].multiplicity_ << "  " << (*p_graph_)[*it].multiplicity_ << endl;
#endif
        to_delete.push_back(*it_e);
      }
    }
  }
  RemoveEdges(to_delete);
  RemoveOrphantVertices();
  return;
}

void DeBruijnGraph::RemoveEdges(std::list<BoostEdge> &to_delete)  {
  // remove edges
  // DEBUG: cout << "To remove " << to_delete.size() << " edges." << endl;
  for(auto it = to_delete.begin(); it != to_delete.end(); ++ it) {
    if(edge(source(*it, *p_graph_), target(*it, *p_graph_), *p_graph_).second)
      boost::remove_edge(*it, *p_graph_);
  }
  // DEBUG: cout << "Remain " << num_edges(*p_graph_) << " edges."  << endl;
  // DEBUG: cout << "Remain " << num_vertices(*p_graph_) << " vertices." << endl;
  return;
}

void DeBruijnGraph::RemoveOrphantVertices() {
  // remove orphant vertices
  list<BoostVertex> vertex_delete;
  auto it_v = boost::vertices(*p_graph_).first;
  while(it_v != boost::vertices(*p_graph_).second) {
    if(degree(*it_v, *p_graph_) <= 0) vertex_delete.push_back(*it_v);
    ++ it_v;
  }
  // DEBUG: cout << "To remove " << vertex_delete.size() << " vertices." << endl;
  for(auto it = vertex_delete.begin(); it != vertex_delete.end(); ++ it) {
    boost::remove_vertex(*it, *p_graph_);
  }
  return;
}


void DeBruijnGraph::TraverseUnbranched(
  std::list<std::string> &paths, const int min_length
) {
  std::list<BoostVertex> sources;
  auto it = boost::edges(*p_graph_).first;
  while(it != boost::edges(*p_graph_).second) {
    if((*p_graph_)[*it].seq_.length() >= min_length)
      paths.push_back((*p_graph_)[*it].seq_);
    ++ it;
  }
  return;
}

void DeBruijnGraph::CondenseGraph()  {
  std::list<BoostEdge> source_edges;
  auto it = boost::vertices(*p_graph_).first;
  while(it != boost::vertices(*p_graph_).second) {
    /*****/
    //int d = boost::in_degree(*it, *p_graph_);
    //cout << "Check vertex:  " << (*p_graph_)[*it].seq_ << " indegree: " << d << endl;
    /*****/
    if(boost::in_degree(*it, *p_graph_) <= 0) {
      /*****/
      //cout << "Source vertex:  " << (*p_graph_)[*it].seq_ << endl;
      /*****/
      auto it_e = boost::out_edges(*it, *p_graph_).first;
      while(it_e != boost::out_edges(*it, *p_graph_).second) {
        source_edges.push_back(*it_e); ++ it_e;
      }
    }
    ++ it;
  }
  
  // condense the graph
  for(auto it = source_edges.begin(); it != source_edges.end(); ++ it) {
    Condense(*it);
  }
  
  // double check if all edges are filled out
  // if not, the edge is a cycle, break it
  auto it_e = boost::edges(*p_graph_).first;
  while(it_e != boost::edges(*p_graph_).second) {
    bool to_del = false; BoostEdge de;
    if((*p_graph_)[*it_e].seq_.length() <= 0)  {
        de = *it_e; to_del = true;
    }
    ++ it_e;
    if(to_del)  boost::remove_edge(de, *p_graph_);
  }
  return;
}

void DeBruijnGraph::Condense(const BoostEdge source_edge) {
  /*****/
  //cout << "Begin of condense graph" << endl;
  /*****/
  std::stack<BoostEdge> to_visit;
  to_visit.push(source_edge);
  while(!to_visit.empty()) {
    BoostEdge init_edge = to_visit.top(); BoostEdge current_edge = init_edge; to_visit.pop(); 
    BoostVertex head = boost::source(init_edge, *p_graph_);
    BoostVertex tail = head;  
    string p = (*p_graph_)[head].seq_.substr(0, kmer_size_);
    int total_multiplicity = (*p_graph_)[head].multiplicity_, total_vertices = 1;
    // construct the path sequence
    do {
      // update the vertex as traversed
      (*p_graph_)[tail].traversed_ = true;
      BoostVertex to_delete = tail;
      
      // define the new tail vertex
      tail = boost::target(current_edge, *p_graph_);
      
      /******/
      //cout << "Tail sequence: " << (*p_graph_)[tail].seq_ << endl;
      /******/
      
      // update the sequence and multiplicity
      p += (*p_graph_)[tail].seq_[kmer_size_ - 1];
      total_multiplicity += (*p_graph_)[tail].multiplicity_; total_vertices ++;
      
      /******/
      //cout << "before edge removed " << endl;
      //cout << "Tail sequence: " << (*p_graph_)[tail].seq_ << endl;
      /******/
      
      // remove the previous edge
      boost::remove_edge(current_edge, *p_graph_);
      
      /******/
      //cout << "edge removed " << endl;
      /******/
      
      // remove the previous tail vertex if it is not the head
      if(to_delete != head) {
        boost::remove_vertex(to_delete, *p_graph_);
        /******/
        //cout << "vertex removed " << endl;
        /******/
      }
      
      // quit condition 1: tail has been visited
      if((*p_graph_)[tail].traversed_) break;
      // quit condition 2: if tail has multiple in edges (remember of traversed edge has been removed)
      if(boost::in_degree(tail, *p_graph_) >= 1) break; 
      // quit condition 3: if tail has multiple out degree (if not, update the current_edge)
      if(boost::out_degree(tail, *p_graph_) == 1) 
        current_edge = *(boost::out_edges(tail, *p_graph_).first);
      else break;
    } while(1);
    
    // connect the head and the tail   
    std::pair<BoostEdge, bool> edge_new = add_edge(head, tail, *p_graph_);
    if(edge_new.second) {
      (*p_graph_)[edge_new.first].seq_ = p;
      (*p_graph_)[edge_new.first].multiplicity_ = total_multiplicity / total_vertices;
    }
    /*****/
    //cout << "condensed sequence:" << p.length() << "  " << p << endl;
    /*****/
    // push source vertices to the stack
    if(boost::out_degree(tail, *p_graph_) > 0 && !(*p_graph_)[tail].traversed_)  {
      auto it = boost::out_edges(tail, *p_graph_).first;
      while(it != boost::out_edges(tail, *p_graph_).second) {
        to_visit.push(*it); ++ it;
      }
    }
  }
  return;
}


bool cmp_bridge(const ShortBridgeType &a, const ShortBridgeType &b) {
  if(a.multiplicity > b.multiplicity)  {
    return true;
  } else if(a.multiplicity == b.multiplicity && a.len <= b.len) {
    return true;
  }
  return false;
}



void DeBruijnGraph::BridgeThreading(const int max_len)  {
  // get all the short bridges
  list<ShortBridgeType> bridges;
  auto it = boost::edges(*p_graph_).first;
  while(it != boost::edges(*p_graph_).second) {
    if((*p_graph_)[*it].seq_.length() <= max_len + kmer_size_) {
      ShortBridgeType b; 
      b.multiplicity = (*p_graph_)[*it].multiplicity_;
      b.len = (*p_graph_)[*it].seq_.length();
      b.source = boost::source(*it, *p_graph_);
      b.target = boost::target(*it, *p_graph_);
      bridges.push_back(b);  
    }
    ++ it;
  }
  
  // sort the bridges based on their multiplicity
  bridges.sort(cmp_bridge);
  // process the bridges
  for(auto it = bridges.begin(); it != bridges.end(); ++ it) {
    // make sure the edge exists
    pair<BoostEdge, bool> e = boost::edge(it->source, it->target, *p_graph_);
    // check for boundary conditions
    if(it->source == it->target)  { 
      boost::remove_edge(e.first, *p_graph_); continue;
    }  
    if(!e.second) continue;
    // do the threading
    if(boost::in_degree(it->source, *p_graph_) > 0 && boost::out_degree(it->source, *p_graph_) > 0)  {
      /*****/
      //cout << "Doing the threading: " << (*p_graph_)[e.first].seq_ << endl;
      /*****/
      auto it_pre = boost::in_edges(it->source, *p_graph_).first;
      while(it_pre != boost::in_edges(it->source, *p_graph_).second) {
        BoostVertex v_pre = boost::source(*it_pre, *p_graph_);
        auto it_suf = boost::out_edges(it->target, *p_graph_).first;
        while(it_suf != boost::out_edges(it->target, *p_graph_).second)  {
          BoostVertex v_suf = boost::target(*it_suf, *p_graph_);
          
          /*****/
          //cout << "==================" << endl;
          //cout << (*p_graph_)[*it_pre].seq_ << endl;
          //cout << (*p_graph_)[e.first].seq_ << endl;
          //cout << (*p_graph_)[*it_suf].seq_ << endl;
          /*****/
          
          // reconstruct the sequence
          string p = (*p_graph_)[*it_pre].seq_ 
              + (*p_graph_)[e.first].seq_.substr(kmer_size_)
              + (*p_graph_)[*it_suf].seq_.substr(kmer_size_);
          
          /*****/
          //cout << "new sequence:  " << p << endl;
          /*****/
          
          // connect the vertices by a new edge
          
          
          pair<BoostEdge, bool> ne = boost::add_edge(v_pre, v_suf, *p_graph_);
          if(ne.second)  {
            (*p_graph_)[ne.first].seq_ = p;
            (*p_graph_)[ne.first].multiplicity_ = (*p_graph_)[e.first].multiplicity_;
          } else 
            cerr << "Warning: DeBruijnGraph::BridgeThreading: failed to add edge to the graph." << endl;
                    
          ++ it_suf;
        }
        ++ it_pre;
      }
      /*****/
      //cout << "end of threading" << endl;
      /*****/
    } else if(boost::in_degree(it->source, *p_graph_) == 0 && 
        boost::out_degree(it->source, *p_graph_) > 0
    ) {
      BoostVertex v_pre = it->source;
      auto it_suf = boost::out_edges(it->target, *p_graph_).first;
      while(it_suf != boost::out_edges(it->target, *p_graph_).second) {
        BoostVertex v_suf = boost::target(*it_suf, *p_graph_);
        string p =  (*p_graph_)[e.first].seq_ + (*p_graph_)[*it_suf].seq_.substr(kmer_size_);
        pair<BoostEdge, bool> ne = boost::add_edge(v_pre, v_suf, *p_graph_);
        if(ne.second)  {
          (*p_graph_)[ne.first].seq_ = p;
          (*p_graph_)[ne.first].multiplicity_ = (*p_graph_)[e.first].multiplicity_;
        } else 
          cerr << "Warning: DeBruijnGraph::BridgeThreading: failed to add edge to the graph." << endl;
        ++ it_suf;
      }
    } else if(boost::in_degree(it->source, *p_graph_) > 0 && 
        boost::out_degree(it->source, *p_graph_) == 0
    ) {
      BoostVertex v_suf = it->target;
      auto it_pre = boost::in_edges(it->source, *p_graph_).first;
      while(it_pre != boost::in_edges(it->source, *p_graph_).second) {
        BoostVertex v_pre = boost::source(*it_pre, *p_graph_);
        string p = (*p_graph_)[*it_pre].seq_ + (*p_graph_)[e.first].seq_.substr(kmer_size_);
        pair<BoostEdge, bool> ne = boost::add_edge(v_pre, v_suf, *p_graph_);
        if(ne.second)  {
          (*p_graph_)[ne.first].seq_ = p;
          (*p_graph_)[ne.first].multiplicity_ = (*p_graph_)[e.first].multiplicity_;
        } else 
          cerr << "Warning: DeBruijnGraph::BridgeThreading: failed to add edge to the graph." << endl;
        ++ it_pre;
      }
    }
   
    /*****/
    //cout << "before edge removal" << endl;
    /*****/
    
    // remove the current edge
    
    boost::remove_edge(e.first, *p_graph_);
    /*****/
    //cout << "edge removed 1" << endl;
    /*****/
    // remove edges if the are already within the new graph
    if(boost::out_degree(it->source, *p_graph_) <= 0) { 
      auto itt = boost::in_edges(it->source, *p_graph_).first;
      while(itt != boost::in_edges(it->source, *p_graph_).second) {
        BoostEdge de = *itt; ++ itt;
        boost::remove_edge(de, *p_graph_);
      }
    }
    if(boost::in_degree(it->target, *p_graph_) <= 0)  { 
      auto itt = boost::out_edges(it->target, *p_graph_).first;
      while(itt != boost::out_edges(it->target, *p_graph_).second) {
        BoostEdge de = *itt; ++ itt;
        boost::remove_edge(de, *p_graph_);
      }
    }
    /*****/
    //cout << "edge removed 2" << endl;
    /*****/
  }
  RemoveOrphantVertices();
  return;
}



void DeBruijnGraph::PrintBranchingVertices() {
  for(auto it = vertices(*p_graph_).first; it != vertices(*p_graph_).second; ++ it) {
    if(in_degree(*it, *p_graph_) > 1 || out_degree(*it, *p_graph_) > 1)  {
      cout << "Node coverage: " << (*p_graph_)[*it].multiplicity_ << " in degree: " << in_degree(*it, *p_graph_) << " out degree: " << out_degree(*it, *p_graph_) << endl;
      if(in_degree(*it, *p_graph_) >= 1)  {
        for(auto it_e = in_edges(*it, *p_graph_).first; it_e != in_edges(*it, *p_graph_).second; ++ it_e)
          cout << "In edge coverage: " << (*p_graph_)[*it_e].multiplicity_ << endl;
      }
      if(out_degree(*it, *p_graph_) >= 1)  {
        for(auto it_e = out_edges(*it, *p_graph_).first; it_e != out_edges(*it, *p_graph_).second; ++ it_e)
          cout << "out edge coverage: " << (*p_graph_)[*it_e].multiplicity_ << endl;
      }
    }
  }
  return;
}

// obsolete function
/*

void DeBruijnGraph::GetShortBridges(
    const int max_len, const int flank_len,
    std::list<ShortBridgeType> &bridges
) {
  auto it = boost::edges(*p_graph_).first;
  while(it != boost::edges(*p_graph_).second) {
    if((*p_graph_)[*it].seq_.length() <= max_len && 
        boost::in_degree(boost::source(*it, *p_graph_), *p_graph_) >= 2 && 
        boost::out_degree(boost::target(*it, *p_graph_), *p_graph_) >= 2
    )  {
      ShortBridgeType b;
      b.multiplicity = (*p_graph_)[*it].multiplicity_;
      b.length = (*p_graph_)[*it].seq_.length();
      b.source = boost::source(*it, *p_graph_);
      b.target = boost::target(*it, *p_graph_);
      // check all source vertices
      auto it_pre = boost::in_edges(b.source, *p_graph_).first;
      while(it_pre != boost::in_edges(b.source, *p_graph_).second) {
        int l = (*p_graph_)[*it_pre].seq_.length();
        if(l >= kmer_size_ + flank_len)  {
          b.pre_vertex.push_back(boost::source(*it_pre, *p_graph_));
          b.pre_string.push_back((*p_graph_)[*it_pre].seq_.substr(l - kmer_size_ - flank_len, kmer_size_));
        }
        ++ it_pre;
      }
      // check all target vertices
      auto it_suf = boost::out_edges(b.target, *p_graph_).first;
      while(it_suf != boost::out_edges(b.target, *p_graph_).second) {
        int l = (*p_graph_)[*it_suf].seq_.length();
        if(l >= kmer_size_ + flank_len)  {
          b.suf_vertex.push_back(boost::target(*it_suf, *p_graph_));
          b.suf_string.push_back((*p_graph_)[*it_suf].seq_.substr(flank_len, kmer_size_));
        }
        ++ it_suf;
      }
      // record the bridge
      if(!b.pre_vertex.empty() && !b.suf_vertex.empty()) bridges.push_back(b);
      

      cout << "bridge sequence: " << (*p_graph_)[*it].seq_ << endl;
      cout << "pre strings: " <<  endl;
      for(auto it_p = b.pre_string.begin(); it_p != b.pre_string.end(); ++ it_p)
        cout << *it_p << endl; 
      cout << "suf strings: " <<  endl;
      for(auto it_p = b.suf_string.begin(); it_p != b.suf_string.end(); ++ it_p)
        cout << *it_p << endl; 

    }  
    ++ it;
  }
  return;
}


void DeBruijnGraph::TraverseUnbranched2(
  std::list<std::string> &paths
) {
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++ it) {
    if((*p_graph_)[*it].traversed_) continue;
    (*p_graph_)[*it].traversed_ = true;
    BoostVertex left_v, right_v;
    left_v = right_v = *it;
    string left_seq = "", right_seq = "";
    
    // right extension
    while(boost::out_degree(right_v, *p_graph_) == 1) {
      BoostEdge out = *(boost::out_edges(right_v, *p_graph_).first);
      right_v = boost::target(out, *p_graph_);
      right_seq += (*p_graph_)[right_v].seq_[kmer_size_ - 1];
      if((*p_graph_)[right_v].traversed_) break;
      (*p_graph_)[right_v].traversed_ = true;
    }
    // left extension
    
    while(boost::in_degree(left_v, *p_graph_) == 1) {
      BoostEdge in = *(boost::in_edges(left_v, *p_graph_).first);
      left_v = boost::source(in, *p_graph_);
      left_seq += (*p_graph_)[left_v].seq_[0];
      if((*p_graph_)[left_v].traversed_) break;
      (*p_graph_)[left_v].traversed_ = true;
    }
    left_seq = string(left_seq.rbegin(), left_seq.rend());  
    
    //paths.push_back(left_seq + (*p_graph_)[*it].seq_ + right_seq);
    cout << ">seq" << endl; 
    cout << left_seq + (*p_graph_)[*it].seq_ + right_seq << endl;
    //cout << (*p_graph_)[*it].seq_ << endl;
  }
  return;
}
*/


/*
void DeBruijnGraph::Construct2(
    const int num_seq, char ** const seq, 
    const FrequencyTable<KmerUnitType> &unit_freq, const int freq_cutoff
) {
  // initilize the object
  p_graph_ = new BoostGraph;
  p_info_map2_ = new HashBucket<KmerUnitType, KmerInfo>;
  initialized_ = true;
  // construct the graph
  int i, j;
  int kmer_id = 0;
  for(i = 0; i < num_seq; ++ i) {
    string s = seq[i];
    for(int j = 0; j < s.length() - kmer_size_; ++ j) {
      string mer = s.substr(j, kmer_size_);
      Kmer km(encoder_, mer);
      KmerInfo kinfo(km, 0); 
      KmerInfo *pv = p_info_map2_->Search(km.kmer_array_[km.size_ - 1], kinfo);
      if(pv == NULL)  p_info_map2_->Insert(km.kmer_array_[km.size_ - 1], kinfo);
      else  ++ pv->count_;
    }
  }
  return;
}
*/

