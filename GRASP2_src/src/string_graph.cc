#include "../include/string_graph.h"

using namespace std;

// import all information in the extension into the string graph
void StringGraph::ImportExtension(std::vector<ExtType> &extension) {
  if(!initialized_) p_graph_ = new BoostSTRGraph;
  initialized_ = true;
  int i, j;
  for(i = 0; i < extension.size(); ++ i) {
    // process the source read
    BoostSTRVertex v_source;
    auto it_source = node_hash_.find(extension[i].source_);
    if(it_source == node_hash_.end()) {
      // adding the new vertex with recorded ID
      STRVertexType node(extension[i].rid_);
      v_source = boost::add_vertex(node, *p_graph_);
      // record the read in node_hash_
      node_hash_[extension[i].source_] = v_source;
    } else  {
      v_source = it_source->second;
      // setting the ID if the corresponding vertex is already in the graph
      (*p_graph_)[v_source].SetID(extension[i].rid_);
    }
    
    //cout << "Vertex ID set: " << (*p_graph_)[v_source].rid_ << endl;
    
    // now process the target reads
    for(j = 0; j < extension[i].ir_position_.size(); ++ j) {
      BoostSTRVertex v_target;
      auto it_target = node_hash_.find(extension[i].ir_position_[j]);
      if(it_target == node_hash_.end()) {
        // adding the new vertex without ID
        STRVertexType node_foo;
        v_target = boost::add_vertex(node_foo, *p_graph_);
        // record the read in node_hash_
        node_hash_[extension[i].ir_position_[j]] = v_target;
      } else  {v_target = it_target->second;}
      // check the two vertices are not the same vertex, avoid self-cycles
      if(v_source == v_target) continue;
      // add corresponding edge between 
      pair<BoostSTREdge, bool> e_search = boost::edge(v_source, v_target, *p_graph_);
      if(!e_search.second)  {e_search = boost::add_edge(v_source, v_target, *p_graph_);} 
      if(e_search.second) {
        // set the corresponding length
        (*p_graph_)[e_search.first].SetLen(extension[i].ir_overlap_[j]);
      } else  {
        cout << "Error:: StringGraph::Construct: Failed to add edges between vertices!" << endl;
      }
    }
  }
  return;
}

void StringGraph::MultiComputeExtension(
    const int threads, const int min_overlap, 
    const int offset, const int n, char **seq, 
    BWT &bwt, BWT &rev_bwt, std::vector<ExtType> &extension
) {
  int i, j;
  // determine the ranges for each chunk
  vector<int> range;
  range.push_back(-1);
  int chunk_size = n / threads;
  for(i = 0; i < threads - 1; ++ i) range.push_back(range.back() + chunk_size);
  range.push_back(n - 1);
  // batch execution
  vector<vector<ExtType> > extension_multi(threads);
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < threads; ++ i) {
      ComputeExtension(
          min_overlap, range[i] + 1, range[i + 1] - range[i], &seq[range[i] + 1], 
          bwt, rev_bwt, extension_multi[i]
      );
    }
    #pragma omp taskwait
  }
  for(i = 0; i < threads; ++ i) {
    for(j = 0; j < extension_multi[i].size(); ++ j) {
      extension.push_back(extension_multi[i][j]); 
    }
  }
  return;
}

void StringGraph::ComputeExtension(
    const int min_overlap, const int offset, const int n, char **seq, 
    BWT &bwt, BWT &rev_bwt, std::vector<ExtType> &extension
) {
  BWTSearch bwt_searcher;
  IvInfo search_info(&bwt, &rev_bwt);  
  int i, j;
  char *r_seq = new char [99999];
  for(i = 0; i < n; ++ i) {
    // we need to search the reverse sequence to identify its unique rev_BWT index
    int n = strlen(seq[i]);
    for(j = 0; j < n; ++ j) r_seq[j] = seq[i][n - j - 1];
    r_seq[n] = '\0';
    AlignType pos;
    if(!bwt_searcher.IsContainedRead(r_seq, rev_bwt, pos)) {
      // recording information regarding the source read
      ExtType ext; ext.source_ = pos.begin; ext.rid_ = offset + i;
      // finding overlapping reads
      bwt_searcher.SearchBeginIntervals(seq[i], min_overlap, search_info); 
      // finding irreducible extensions
      bwt_searcher.FindIrreducible(search_info, ext.ir_position_, ext.ir_overlap_);  
      // resetting the search information (list)
      search_info.Reset();   
      // record the output
      extension.push_back(ext);
    }
  }
  delete [] r_seq;
  return;
}

void StringGraph::RemoveOrphantVertices() {
  // remove orphant vertices
  auto it_v = boost::vertices(*p_graph_).first;
  while(it_v != boost::vertices(*p_graph_).second) {
    BoostSTRVertex v = *it_v; ++ it_v;
    if(degree(v, *p_graph_) <= 0) boost::remove_vertex(v, *p_graph_);
  }
  return;
}

int StringGraph::RemoveTipsBeforeCondense()  {
  auto it_v = boost::vertices(*p_graph_).first;
  vector<BoostSTRVertex> to_delete;
  while(it_v != boost::vertices(*p_graph_).second) {
    BoostSTRVertex v = *it_v; ++ it_v;
    if(boost::out_degree(v, *p_graph_) == 0 && boost::in_degree(v, *p_graph_) == 1) {
      // the in_edge of that vertex
      BoostSTREdge e = *(boost::in_edges(v, *p_graph_).first); int l = (*p_graph_)[e].len_;
      // the source vertex for the in_edge
      BoostSTRVertex s = boost::source(e, *p_graph_);
      if(boost::out_degree(s, *p_graph_) > 1)  {
        // if the source contains multiple out_edges, only reserve the one with longest overlap
        auto it_e = boost::out_edges(s, *p_graph_).first;
        while(it_e != boost::out_edges(s, *p_graph_).second) {
          if((*p_graph_)[*it_e].len_ > l || 
              boost::out_degree(boost::target(*it_e, *p_graph_), *p_graph_) > 0
          ) {
            to_delete.push_back(v); break;
          }
          ++ it_e;
        }
      }
    } else if(boost::out_degree(v, *p_graph_) == 0 && boost::in_degree(v, *p_graph_) > 1) {
      // if multiple paths end at one vertex, remove the vertex
      to_delete.push_back(v);
    } else if(boost::in_degree(v, *p_graph_) == 0 && boost::out_degree(v, *p_graph_) == 1) {
      // the out_edge of that vertex
      BoostSTREdge e = *(boost::out_edges(v, *p_graph_).first); int l = (*p_graph_)[e].len_;
      // the target vertex for the out_edge
      BoostSTRVertex t = boost::target(e, *p_graph_);
      if(boost::in_degree(t, *p_graph_) > 1)  {
        // if the target contains multiple in_edges, only reserve the one with longest overlap
        auto it_e = boost::in_edges(t, *p_graph_).first;
        while(it_e != boost::in_edges(t, *p_graph_).second) {
          if((*p_graph_)[*it_e].len_ > l ||
              boost::in_degree(boost::source(*it_e, *p_graph_), *p_graph_) > 0    
          )  {
            to_delete.push_back(v); break;
          }
          ++ it_e;
        }
      }
    } else if(boost::in_degree(v, *p_graph_) == 0 && boost::out_degree(v, *p_graph_) > 1) {
      // if multiple paths begin at one vertex, remove the vertex
      to_delete.push_back(v);
    }
  }
  // remove the edges and the vertices
  BoostSTREdge de;
  for(auto it = to_delete.begin(); it != to_delete.end(); ++ it) {
    auto it_ie = boost::in_edges(*it, *p_graph_).first;
    while(it_ie != boost::in_edges(*it, *p_graph_).second) {
      de = *it_ie; ++ it_ie;
      boost::remove_edge(de, *p_graph_);
    }
    auto it_oe = boost::out_edges(*it, *p_graph_).first;
    while(it_oe != boost::out_edges(*it, *p_graph_).second) {
      de = *it_oe; ++ it_oe;
      boost::remove_edge(de, *p_graph_);;
    }
    //boost::remove_vertex(*it, *p_graph_);
  }
  return to_delete.size();
}

// check if there exists some un-labeled vertex or un-labeld path
// delete those un-labeld vertex or paths and clear the graph
void StringGraph::CheckGraph(void)  {
  auto it_v = boost::vertices(*p_graph_).first;
  while(it_v != boost::vertices(*p_graph_).second) {
    BoostSTRVertex v = *it_v; ++ it_v;
    if((*p_graph_)[v].rid_ < 0) {
      // remove associated edges first
      auto it_ie = boost::in_edges(v, *p_graph_).first;
      while(it_ie != boost::in_edges(v, *p_graph_).second) {
        BoostSTREdge ie = *it_ie; ++ it_ie;
        boost::remove_edge(ie, *p_graph_);
      }
      auto it_oe = boost::out_edges(v, *p_graph_).first;
      while(it_oe != boost::out_edges(v, *p_graph_).second) {
        BoostSTREdge oe = *it_oe; ++ it_oe;
        boost::remove_edge(oe, *p_graph_);
      }
      // remove the vertex
      boost::remove_vertex(v, *p_graph_);
      //cout << "Erroneous vertex removed" << endl;
    }
  }
  
  auto it_e = boost::edges(*p_graph_).first;
  while(it_e != boost::edges(*p_graph_).second) {
    BoostSTREdge e = *it_e; ++ it_e;
    if((*p_graph_)[e].len_ < 0 && (*p_graph_)[e].seq_.empty()) {
      boost::remove_edge(e, *p_graph_);
      //cout << "Erroneous edge removed" << endl;
    }
  } 
  return;
}

void StringGraph::CheckSelfCycle(void)  {
  vector<BoostSTREdge> to_delete;
  auto it = boost::edges(*p_graph_).first;
  while(it != boost::edges(*p_graph_).second) {
    if(boost::source(*it, *p_graph_) == boost::target(*it, *p_graph_)) 
      to_delete.push_back(*it);
    ++ it;
  }
  for(int i = 0; i < to_delete.size(); ++ i) {
    boost::remove_edge(to_delete[i], *p_graph_);
  }
  return;
}

void StringGraph::CondenseGraph(char** seq)  {
  std::list<BoostSTREdge> source_edges;
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
    Condense(seq, *it);
  }
  
  // double check if all edges are filled out
  // if not, the edge is a cycle, break it
  auto it_e = boost::edges(*p_graph_).first;
  while(it_e != boost::edges(*p_graph_).second) {
    bool to_del = false; BoostSTREdge de;
    if((*p_graph_)[*it_e].seq_.length() <= 0)  {
        de = *it_e; to_del = true;
    }
    ++ it_e;
    if(to_del)  boost::remove_edge(de, *p_graph_);
  }
  return;
}

void StringGraph::Condense(char** seq, const BoostSTREdge source_edge) {
  std::stack<BoostSTREdge> to_visit;
  to_visit.push(source_edge);
  while(!to_visit.empty()) {
    BoostSTREdge init_edge = to_visit.top(); BoostSTREdge current_edge = init_edge; to_visit.pop(); 
    BoostSTRVertex head = boost::source(init_edge, *p_graph_);
    BoostSTRVertex tail = boost::target(init_edge, *p_graph_);  
    if(head == tail) continue;
    else tail = head;
    string p = seq[(*p_graph_)[head].rid_];
    int total_vertices = 1;
    vector<int> path_info; path_info.push_back((*p_graph_)[head].rid_);
    // construct the path sequence
    do {
      // update the vertex as traversed
      (*p_graph_)[tail].SetTraversedTag(true);
      BoostSTRVertex to_delete = tail;     
      // define the new tail vertex
      tail = boost::target(current_edge, *p_graph_);   
      // record the path information
      path_info.push_back((*p_graph_)[current_edge].len_); 
      path_info.push_back((*p_graph_)[tail].rid_);  
      // update the sequence and multiplicity
      p += &seq[(*p_graph_)[tail].rid_][(*p_graph_)[current_edge].len_];
      ++ total_vertices;      
      // remove the previous edge
      boost::remove_edge(current_edge, *p_graph_);
      // remove the previous tail vertex if it is not the head
      if(to_delete != head) {
        boost::remove_vertex(to_delete, *p_graph_);
      }
      // quit condition 1: tail has been visited
      if((*p_graph_)[tail].IsTraversed()) break;
      // quit condition 2: if tail has multiple in edges (remember of traversed edge has been removed)
      if(boost::in_degree(tail, *p_graph_) >= 1) break; 
      // quit condition 3: if tail has multiple out degree (if not, update the current_edge)
      if(boost::out_degree(tail, *p_graph_) == 1) 
        current_edge = *(boost::out_edges(tail, *p_graph_).first);
      else break;
    } while(1);
    
    // connect the head and the tail   
    if(head != tail)  {
      std::pair<BoostSTREdge, bool> edge_new = add_edge(head, tail, *p_graph_);
      if(edge_new.second) {
        (*p_graph_)[edge_new.first].SetSeq(p);
        (*p_graph_)[edge_new.first].path_info_ = path_info;
      }
    }
    // push source vertices to the stack
    if(boost::out_degree(tail, *p_graph_) > 0 && !(*p_graph_)[tail].IsTraversed())  {
      auto it = boost::out_edges(tail, *p_graph_).first;
      while(it != boost::out_edges(tail, *p_graph_).second) {
        to_visit.push(*it); ++ it;
      }
    }
  }
  return;
}


/*********************************************
void StringGraph::ClusteringSeqs(const int identity)  {
  int i, j;
  AlignBatch aligner;
  auto it = boost::vertices(*p_graph_).first;
  while(it != boost::vertices(*p_graph_).second) {
    if(boost::out_degree(*it, *p_graph_) > 1) {

      cout << "Branching case detected: " << endl;

      vector<string> seqs;
      auto it_e = boost::out_edges(*it, *p_graph_).first;
      while(it_e != boost::out_edges(*it, *p_graph_).second) {
        seqs.push_back((*p_graph_)[*it_e].seq_); 
        cout << (*p_graph_)[*it_e].seq_ << endl;
        int num_out = boost::out_degree(boost::target(*it_e, *p_graph_), *p_graph_);
        cout << "num out degree:  " << num_out << endl;
        ++ it_e;
      }
      cout << "Num of outgoing sequences: " << seqs.size() << endl;
      
      for(i = 0; i < seqs.size() - 1; ++ i) {
        for(j = i + 1; j < seqs.size(); ++ j) {
          cout << "alignment begins" << endl;
          int l = seqs[i].length() < seqs[j].length() ? seqs[i].length() : seqs[j].length();
          int b = (l * (100 - identity) / 100) + 1;
          int ed = aligner.AlignGlobalSingleEditDist(seqs[i].substr(0, l), seqs[j].substr(0,l), 2 * b);
          cout << "alignment: " << i << " " << j << " " << ed << endl;
        }
      }
    }
    ++ it;
  }
  return;
}
**********************************************/

void StringGraph::TraverseUnbranched(
  std::list<std::string> &paths, const int min_length
) {
  std::list<BoostSTRVertex> sources;
  auto it = boost::edges(*p_graph_).first;
  while(it != boost::edges(*p_graph_).second) {
    if((*p_graph_)[*it].seq_.length() >= min_length)
      paths.push_back((*p_graph_)[*it].seq_);
    ++ it;
  }
  return;
}

void StringGraph::WriteGraph(
    BioAlphabet &alphabet, 
    char** seq, const std::string &file_name
) {
  ofstream out_file;
  out_file.open(file_name, ios::out);
  if(!out_file.is_open())  {
    cout << "StringGraph::WriteGraph: Error in writing indexing file " << file_name << "; Abort." << endl;
  }
  KmerUnitcoder coder(alphabet, 6);
  MinimizerSort m_sorter;
  
  // writing all edges
  int e_size = boost::num_edges(*p_graph_);
  char **pheader = new char* [e_size];
  char **pseq = new char* [e_size];
  int e_index = 0;
  auto ite = boost::edges(*p_graph_).first;
  while(ite != boost::edges(*p_graph_).second) {
    stringstream out_header;
    BoostSTRVertex s = boost::source(*ite, *p_graph_);
    BoostSTRVertex t = boost::target(*ite, *p_graph_);
    out_header << ">" << (*p_graph_)[s].rid_ << ":" << strlen(seq[(*p_graph_)[s].rid_]) 
        << ":" << (*p_graph_)[t].rid_ << ":" << strlen(seq[(*p_graph_)[t].rid_]);
    for(int i = 0; i < (*p_graph_)[*ite].path_info_.size(); ++ i) {
      out_header << ":" << (*p_graph_)[*ite].path_info_[i];
    }
    pheader[e_index] = new char [out_header.str().length() + 1];
    strcpy(pheader[e_index], out_header.str().c_str());
    pseq[e_index] = new char [(*p_graph_)[*ite].seq_.length() + 1];
    strcpy(pseq[e_index], (*p_graph_)[*ite].seq_.c_str()); 
    ++ e_index; ++ ite;
  }  
  m_sorter.SortSeqs(coder, e_index, pheader, pseq);
  for(int i = 0; i < e_index; ++ i) {
    out_file << pheader[i] << endl;
    out_file << pseq[i] << endl;
    delete [] pheader[i]; delete [] pseq[i];
  }
  delete [] pheader; delete [] pseq;
  // writing all orphant reads
  int v_index = 0;
  auto itv = boost::vertices(*p_graph_).first;
  while(itv != boost::vertices(*p_graph_).second) {
    if(degree(*itv, *p_graph_) <= 0)  {
      ++ v_index;
    }
    ++ itv;
  }
  pheader = new char* [v_index]; pseq = new char* [v_index];
  
  v_index = 0;
  itv = boost::vertices(*p_graph_).first;
  while(itv != boost::vertices(*p_graph_).second) {
    if(degree(*itv, *p_graph_) <= 0)  {
      stringstream out_header;
      out_header << ">" << (*p_graph_)[*itv].rid_ 
          << ":" << strlen(seq[(*p_graph_)[*itv].rid_]) 
          << ":" << (*p_graph_)[*itv].rid_ << ":" 
          << strlen(seq[(*p_graph_)[*itv].rid_]);
      pheader[v_index] = new char [out_header.str().length() + 1];
      strcpy(pheader[v_index], out_header.str().c_str());
      pseq[v_index] = new char [strlen(seq[(*p_graph_)[*itv].rid_]) + 1];
      strcpy(pseq[v_index], seq[(*p_graph_)[*itv].rid_]);
      ++ v_index;
    }
    ++ itv;
  }
  m_sorter.SortSeqs(coder, v_index, pheader, pseq);
  for(int i = 0; i < v_index; ++ i) {
    out_file << pheader[i] << endl;
    out_file << pseq[i] << endl;
    delete [] pheader[i]; delete [] pseq[i];
  }
  delete [] pheader; delete [] pseq;

  out_file.close();
  return;
}

void StringGraph::LoadGraph(
    const std::string &file_name, 
    std::vector<int> &orphan_rid, std::vector<std::string> &orphan_seq
) {
  if(!initialized_) p_graph_ = new BoostSTRGraph;
  initialized_ = true;
  string line;
  ifstream in_file;
  in_file.open(file_name, ios::in);
  if(!in_file.is_open())  {
    cout << "StringGraph::LoadGraph: Error in loading indexing file " << file_name << "; Abort." << endl;
  }
  // define the vertex hash (key is read ID; value is reference to the vertex in graph)
  unordered_map<int, BoostSTRVertex> vertex_hash;
  // loading the files line by line
  while(getline(in_file, line))  {
    if(line[0] == '>')  {   // a proper header
      // load the sequence information
      string seq; getline(in_file, seq);     
      // extract the header information
      vector<int> range;
      range.push_back(0);
      for(int i = 0; i < line.length(); ++ i) {
        if(line[i] == ':') range.push_back(i);
      }
      range.push_back(line.length());
      // the range should contain exactly 5 entries that corresponds to 4 fields in the header
      // if not, discard the current sequence
      if(range.size() < 5) continue;
      // decode the read header information
      int sid = stoi(line.substr(range[0] + 1, range[1] - range[0] - 1)),
          slen = stoi(line.substr(range[1] + 1, range[2] - range[1] - 1)),
          tid = stoi(line.substr(range[2] + 1, range[3] - range[2] - 1)),
          tlen = stoi(line.substr(range[3] + 1, range[4] - range[3] - 1));   
           
      // if the current sequence corresponds to an orphan vertex
      if(sid == tid)  {
        orphan_rid.push_back(sid); orphan_seq.push_back(seq);
        continue;
      }
      // add or locate the two vertices in the graph
      BoostSTRVertex v_source, v_target;
      auto it_source = vertex_hash.find(sid);
      if(it_source == vertex_hash.end()) {
        STRVertexType node(sid); node.len_ = slen;
        v_source = boost::add_vertex(node, *p_graph_);
        vertex_hash[sid] = v_source;
      } else  {
        v_source = it_source->second;  
      }
      auto it_target = vertex_hash.find(tid);
      if(it_target == vertex_hash.end()) {
        STRVertexType node(tid); node.len_ = tlen;
        v_target = boost::add_vertex(node, *p_graph_);
        vertex_hash[tid] = v_target;
      } else  {
        v_target = it_target->second;
      }
      // add the edge between these vertices
      pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_target, *p_graph_);
      if(e_search.second) {
        // set the corresponding length
        (*p_graph_)[e_search.first].SetSeq(seq);
        // adding back the path information
        for(int i = 4; i < range.size() - 1; ++ i) {         
          (*p_graph_)[e_search.first].path_info_.push_back(
              stoi(line.substr(range[i] + 1, range[i + 1] - range[i] - 1))
          );
        }        

      } else  {
        cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!" << endl;
      }
    }
  }
  in_file.close();
  return;
}

int StringGraph::RecordEdgeSeqs(std::vector<std::string> &seqs) {
  int n = 0;
  auto ite = boost::edges(*p_graph_).first;
  while(ite != boost::edges(*p_graph_).second) {
    // setting the read ID
    (*p_graph_)[*ite].sid_ = n; 
    // copying the sequence
    seqs.push_back((*p_graph_)[*ite].seq_);
    ++ n; ++ ite;
  }
  return n;
}

void StringGraph::ComputeGraphEdgeMapping(std::vector<BoostSTREdge> &graph_edge) {
  int n = boost::num_edges(*p_graph_);
  graph_edge.resize(n);
  auto ite = boost::edges(*p_graph_).first;
  while(ite != boost::edges(*p_graph_).second) {
    graph_edge[(*p_graph_)[*ite].sid_] = *ite;
    ++ ite;
  }
  return;
}


bool SortSeqAlnInfoType(const SeqAlnInfoType &i1, const SeqAlnInfoType &i2) {
  if(i1.score > i2.score || (i1.score == i2.score && (i1.q2 - i1.q1) < (i2.q2 - i2.q1)))
    return true;
  return false;  
}

bool SortByScore(
  const std::pair<std::string, int> &s1, 
  const std::pair<std::string, int> &s2
) {
  if(s1.second > s2.second || 
      (s1.second == s2.second && s1.first.length() > s2.first.length())
  )  return true;
  return false;
}

void StringGraph::GetHighScoringPaths(
    const int cutoff, const int query_len,
    std::vector<std::string> &edge_seqs,
    std::vector<BoostSTREdge> &graph_edge,
    std::vector<int> &edge_ID, std::vector<int> &score, 
    std::vector<std::pair<int, int> > &q_interval,
    const double sim_cutoff,
    std::vector<std::string> &high_scoring_seqs
) {
  if(edge_ID.size() != score.size())  {
    cout << "Error: StringGraph::GetHighScoringPaths: inconsistent number of sequence IDs and alignment scores, abort." << endl;
    exit(1);
  }
  int i, j;
  // first sort the sequences based on their alignment scores
  list<SeqAlnInfoType> seed_seqs;
  unordered_map<int, int> aligned_seqs;
  for(i = 0; i < edge_ID.size(); ++ i)  {
    SeqAlnInfoType ai; 
    ai.id = edge_ID[i]; ai.score = score[i];
    ai.q1 = q_interval[i].first; ai.q2 = q_interval[i].second;
    seed_seqs.push_back(ai);
    if(score[i] > 0 && 
        (aligned_seqs.find(edge_ID[i]) == aligned_seqs.end() || 
         aligned_seqs[edge_ID[i]] < score[i])
    )  {
      aligned_seqs[edge_ID[i]] = score[i];
    }
  }
  // the sort the alignment information
  seed_seqs.sort(SortSeqAlnInfoType);
  
  //for(auto it = seed_seqs.begin(); it != seed_seqs.end(); ++ it) {
  //  cout << ">" << it->score << endl; 
  //}
  
  // traverse the graph using each alignment as seed
  unordered_map<int, int> visited;
  for(auto it = seed_seqs.begin(); it != seed_seqs.end(); ++ it) {
    // skip the sequences that have been traversed
    if(visited.find(it->id) != visited.end())  continue;
    if(it->id >= graph_edge.size()) {
      high_scoring_seqs.push_back(edge_seqs[it->id]);
      continue;
    }
    if(it->score < cutoff) break;
    // use the current edge sequence as the seed
    BoostSTREdge seed_edge = graph_edge[it->id];
    // extend in both directions
    //cout << "left extension called" << endl;
    vector<pair<string, int> > ext_seqs_left, ext_seqs_right;
    ProgressiveExtendLeft(
        seed_edge, it->q1, edge_seqs, 
        aligned_seqs, ext_seqs_left, visited
    );
    //cout << "left extension done" << endl;
    
    vector<pair<string, int> > nr_left_seq;
    if(ext_seqs_left.size() > 2)  {
      sort(ext_seqs_left.begin(), ext_seqs_left.end(), SortByScore);
      vector<bool> keep_seq(ext_seqs_left.size(), true);
      for(i = 0; i < ext_seqs_left.size() - 1; ++ i) {
        if(!keep_seq[i]) continue;
        for(j = i + 1; j < ext_seqs_left.size(); ++ j) {
          double sim = EstSeqSimilarity(
              ext_seqs_left[i].first, ext_seqs_left[j].first, false
          );
          if(sim >= sim_cutoff) keep_seq[j] = false;
        }
      }
      for(i = 0; i < ext_seqs_left.size(); ++ i) {
        if(keep_seq[i]) nr_left_seq.push_back(ext_seqs_left[i]);
      }
    } else  {nr_left_seq = ext_seqs_left;}
    
    //cout << "right extension called" << endl;
    ProgressiveExtendRight(
        seed_edge, query_len - it->q2, edge_seqs, 
        aligned_seqs, ext_seqs_right, visited
    );
    //cout << "right extension done" << endl;
    
    vector<pair<string, int> > nr_right_seq;
    if(ext_seqs_right.size() > 2)  {
      sort(ext_seqs_right.begin(), ext_seqs_right.end(), SortByScore);
      vector<bool> keep_seq(ext_seqs_right.size(), true);
      for(i = 0; i < ext_seqs_right.size() - 1; ++ i) {
        if(!keep_seq[i]) continue;
        for(j = i + 1; j < ext_seqs_right.size(); ++ j) {
          double sim = EstSeqSimilarity(
              ext_seqs_right[i].first, ext_seqs_right[j].first, true
          );
          if(sim >= sim_cutoff) keep_seq[j] = false;
        }
      }
      for(i = 0; i < ext_seqs_right.size(); ++ i) {
        if(keep_seq[i]) nr_right_seq.push_back(ext_seqs_right[i]);
      }
    } else  {nr_right_seq = ext_seqs_right;}
 
    // merge the extensions   
    //cout << "merge extension called" << endl;
    MergeExtensions(
        edge_seqs[it->id], 
        nr_left_seq, nr_right_seq, 
        high_scoring_seqs
    );
    //cout << "merge extension done" << endl;
    
    //**************************
    //cout << "flanking sizes:  " << nr_left_seq.size() << "  " << nr_right_seq.size() << endl;
    //cout << "core sequence: " << edge_seqs[it->id] << endl;
    //for(auto its = high_scoring_seqs.begin(); its != high_scoring_seqs.end(); ++ its) {
    //  cout << "???: " << it->id << "  " << *its << endl; 
    //}
  }
  return;
}

void StringGraph::GetHighScoringPaths(
    const int cutoff, const int query_len,
    std::vector<std::string> &edge_seqs,
    std::vector<BoostSTREdge> &graph_edge,
    const int num_seqs, int *edge_ID, int *score, 
    int *q_interval, const double sim_cutoff,
    std::vector<std::string> &high_scoring_seqs
) {
  int i, j;
  // first sort the sequences based on their alignment scores
  list<SeqAlnInfoType> seed_seqs;
  unordered_map<int, int> aligned_seqs;
  for(i = 0; i < num_seqs; ++ i)  {
    SeqAlnInfoType ai; 
    ai.id = edge_ID[i]; ai.score = score[i];
    ai.q1 = q_interval[2 * i]; ai.q2 = q_interval[2 * i + 1];
    seed_seqs.push_back(ai);
    if(score[i] > 0 && 
        (aligned_seqs.find(edge_ID[i]) == aligned_seqs.end() || 
         aligned_seqs[edge_ID[i]] < score[i])
    )  {
      aligned_seqs[edge_ID[i]] = score[i];
    }
  }
  // the sort the alignment information
  seed_seqs.sort(SortSeqAlnInfoType);
  
  //for(auto it = seed_seqs.begin(); it != seed_seqs.end(); ++ it) {
  //  cout << ">" << it->score << endl; 
  //}
  
  // traverse the graph using each alignment as seed
  unordered_map<int, int> visited;
  for(auto it = seed_seqs.begin(); it != seed_seqs.end(); ++ it) {
    // skip the sequences that have been traversed
    if(visited.find(it->id) != visited.end())  continue;
    if(it->id >= graph_edge.size()) {
      high_scoring_seqs.push_back(edge_seqs[it->id]);
      continue;
    }
    if(it->score < cutoff) break;
    // use the current edge sequence as the seed
    BoostSTREdge seed_edge = graph_edge[it->id];
    // extend in both directions
    //cout << "left extension called" << endl;
    vector<pair<string, int> > ext_seqs_left, ext_seqs_right;
    ProgressiveExtendLeft(
        seed_edge, it->q1, edge_seqs, 
        aligned_seqs, ext_seqs_left, visited
    );
    //cout << "left extension done" << endl;
    
    vector<pair<string, int> > nr_left_seq;
    if(ext_seqs_left.size() > 2)  {
      sort(ext_seqs_left.begin(), ext_seqs_left.end(), SortByScore);
      vector<bool> keep_seq(ext_seqs_left.size(), true);
      for(i = 0; i < ext_seqs_left.size() - 1; ++ i) {
        if(!keep_seq[i]) continue;
        for(j = i + 1; j < ext_seqs_left.size(); ++ j) {
          double sim = EstSeqSimilarity(
              ext_seqs_left[i].first, ext_seqs_left[j].first, false
          );
          if(sim >= sim_cutoff) keep_seq[j] = false;
        }
      }
      for(i = 0; i < ext_seqs_left.size(); ++ i) {
        if(keep_seq[i]) nr_left_seq.push_back(ext_seqs_left[i]);
      }
    } else  {nr_left_seq = ext_seqs_left;}
    
    //cout << "right extension called" << endl;
    ProgressiveExtendRight(
        seed_edge, query_len - it->q2, edge_seqs, 
        aligned_seqs, ext_seqs_right, visited
    );
    //cout << "right extension done" << endl;
    
    vector<pair<string, int> > nr_right_seq;
    if(ext_seqs_right.size() > 2)  {
      sort(ext_seqs_right.begin(), ext_seqs_right.end(), SortByScore);
      vector<bool> keep_seq(ext_seqs_right.size(), true);
      for(i = 0; i < ext_seqs_right.size() - 1; ++ i) {
        if(!keep_seq[i]) continue;
        for(j = i + 1; j < ext_seqs_right.size(); ++ j) {
          double sim = EstSeqSimilarity(
              ext_seqs_right[i].first, ext_seqs_right[j].first, true
          );
          if(sim >= sim_cutoff) keep_seq[j] = false;
        }
      }
      for(i = 0; i < ext_seqs_right.size(); ++ i) {
        if(keep_seq[i]) nr_right_seq.push_back(ext_seqs_right[i]);
      }
    } else  {nr_right_seq = ext_seqs_right;}
 
    // merge the extensions   
    //cout << "merge extension called" << endl;
    MergeExtensions(
        edge_seqs[it->id], 
        nr_left_seq, nr_right_seq, 
        high_scoring_seqs
    );
    //cout << "merge extension done" << endl;
    
    //**************************
    //cout << "flanking sizes:  " << nr_left_seq.size() << "  " << nr_right_seq.size() << endl;
    //cout << "core sequence: " << edge_seqs[it->id] << endl;
    //for(auto its = high_scoring_seqs.begin(); its != high_scoring_seqs.end(); ++ its) {
    //  cout << "???: " << it->id << "  " << *its << endl; 
    //}
  }
  return;
}

double StringGraph::EstSeqSimilarity(
    const std::string &seq1, const std::string &seq2, const bool is_left
) {
  int n = seq1.length() < seq2.length() ? seq1.length() : seq2.length();
  int same = 0, total = 0;
  if(is_left)  {
    for(int i = 0; i < n; ++ i) {
      if(seq1[i] == seq2[i]) ++ same;
      ++ total;
    }
  } else  {
    if(seq1.length() > seq2.length()) {  
      int k = seq1.length() - seq2.length();
      for(int i = 0; i < n; ++ i) {
        if(seq1[k + i] == seq2[i]) ++ same;
        ++ total;
      }
    } else  {
      int k = seq2.length() - seq1.length();
      for(int i = 0; i < n; ++ i) {
        if(seq1[i] == seq2[k + i]) ++ same;
        ++ total;
      }
    }
  }
  return (double) same / total;
}

void StringGraph::ProgressiveExtendLeft(
    BoostSTREdge &source_edge, 
    const int expected_len, std::vector<std::string> &seqs,
    std::unordered_map<int, int> &aligned_seqs, 
    std::vector<std::pair<std::string, int> > &expanded_seqs,
    std::unordered_map<int, int> &visited_seqs
) {
  if(expected_len <= 0 || visited_seqs.find((*p_graph_)[source_edge].sid_) != visited_seqs.end()) return;
  // the list contains a set of edges that constitute the path
  stack<list<BoostSTREdge> > path_edges;
  // the length of the current path
  stack<int> sum_length;
  // initialized the traversal
  // handle edge ID
  list<BoostSTREdge> init_path; init_path.push_back(source_edge); path_edges.push(init_path);
  // handle length
  sum_length.push(seqs[(*p_graph_)[source_edge].sid_].length());
  // BFS
  while(!path_edges.empty()) {
    list<BoostSTREdge> current_path = path_edges.top(); path_edges.pop();
    int current_length = sum_length.top(); sum_length.pop();
    BoostSTREdge end_edge = current_path.back();
    // notice we are extending left, so get the source vertex
    BoostSTRVertex end_vertex = source(end_edge, *p_graph_);
    int ovl = (*p_graph_)[end_vertex].len_;
    int num_extensions = 0;
    auto itie = in_edges(end_vertex, *p_graph_).first;
    while(itie != in_edges(end_vertex, *p_graph_).second) { 
      int eid = (*p_graph_)[*itie].sid_; 
      // check if the edge has been aligned
      if(visited_seqs.find(eid) != visited_seqs.end())  { 
          ++ itie; continue; 
      }
      // record that the current edge is visited
      visited_seqs[eid] = aligned_seqs.find(eid) == aligned_seqs.end() ? 
          0 : aligned_seqs[eid];
      // construct information for the next stage       
      int next_len = current_length + seqs[eid].length() - ovl;
      list<BoostSTREdge> next_path = current_path; 
      next_path.push_back(*itie);     
      
      // wether to record the current extension information
      if(next_len >= expected_len)  {
        string concat_seq = "";
        int path_score = 0;
        // record that the vertices were stored from C-terminal to N-terminal
        // so traverse it in reversed order
        for(auto itpe = next_path.rbegin(); itpe != next_path.rend(); ++ itpe) {
          string s = seqs[(*p_graph_)[*itpe].sid_];
          BoostSTRVertex v = source(*itpe, *p_graph_);
          if(concat_seq == "")  { concat_seq = s; } 
          else  {
            int pre_len = (*p_graph_)[v].len_; 
            concat_seq += s.substr(pre_len, s.length() - pre_len);  
          }
          path_score += visited_seqs[(*p_graph_)[*itpe].sid_];
        }
        expanded_seqs.push_back(make_pair(concat_seq, path_score));
      } else  {
        // if the length is shorter than the expected length, push the current
        // information to the stack
        path_edges.push(next_path);
        sum_length.push(next_len);
      }
      ++ itie; ++ num_extensions;
    }
    // if the current path is a terminal, record it
    if(num_extensions <= 0 && current_path.size() > 1)  {
      string concat_seq = "";
      int path_score = 0;
      // record that the vertices were stored from C-terminal to N-terminal
      // so traverse it in reversed order
      for(auto itpe = current_path.rbegin(); itpe != current_path.rend(); ++ itpe) {
        string s = seqs[(*p_graph_)[*itpe].sid_];
        BoostSTRVertex v = source(*itpe, *p_graph_);
        if(concat_seq == "")  { concat_seq = s; } 
        else  {
          int pre_len = (*p_graph_)[v].len_; 
          concat_seq += s.substr(pre_len, s.length() - pre_len);
        }
        path_score += visited_seqs[(*p_graph_)[*itpe].sid_];
      }
      expanded_seqs.push_back(make_pair(concat_seq, path_score));
    }
  }
  return;
}
  
void StringGraph::ProgressiveExtendRight(
    BoostSTREdge &source_edge, 
    const int expected_len, std::vector<std::string> &seqs,
    std::unordered_map<int, int> &aligned_seqs, 
    std::vector<std::pair<std::string, int> > &expanded_seqs,
    std::unordered_map<int, int> &visited_seqs
) {
    if(expected_len <= 0 || visited_seqs.find((*p_graph_)[source_edge].sid_) != visited_seqs.end()) return;
  // the list contains a set of edges that constitute the path
  stack<list<BoostSTREdge> > path_edges;
  // the length of the current path
  stack<int> sum_length;
  // initialized the traversal
  // handle edge ID
  list<BoostSTREdge> init_path; init_path.push_back(source_edge); path_edges.push(init_path);
  // handle length
  sum_length.push(seqs[(*p_graph_)[source_edge].sid_].length());
  // BFS
  while(!path_edges.empty()) {
    list<BoostSTREdge> current_path = path_edges.top(); path_edges.pop();
    int current_length = sum_length.top(); sum_length.pop();
    BoostSTREdge end_edge = current_path.back();
    // notice we are extending right, so get the target vertex
    BoostSTRVertex end_vertex = target(end_edge, *p_graph_);
    int ovl = (*p_graph_)[end_vertex].len_;
    int num_extensions = 0;
    auto itoe = out_edges(end_vertex, *p_graph_).first;
    while(itoe != out_edges(end_vertex, *p_graph_).second) { 
      int eid = (*p_graph_)[*itoe].sid_; 
      // check if the edge has been aligned
      if(visited_seqs.find(eid) != visited_seqs.end())  { 
          ++ itoe; continue; 
      }
      // record that the current edge is visited
      visited_seqs[eid] = aligned_seqs.find(eid) == aligned_seqs.end() ? 
          0 : aligned_seqs[eid];
      // construct information for the next stage       
      int next_len = current_length + seqs[eid].length() - ovl;
      list<BoostSTREdge> next_path = current_path; 
      next_path.push_back(*itoe);     
      // wether to record the current extension information
      if(next_len >= expected_len)  {
        string concat_seq = ""; int path_score = 0;
        // record that the vertices were stored from N-terminal to C-terminal
        // so traverse it in forward order
        for(auto itpe = next_path.begin(); itpe != next_path.end(); ++ itpe) {
          string s = seqs[(*p_graph_)[*itpe].sid_];
          BoostSTRVertex v = target(*itpe, *p_graph_);
          if(concat_seq == "")  { concat_seq = s; } 
          else  {
            int post_len = (*p_graph_)[v].len_; 
            concat_seq += s.substr(post_len, s.length() - post_len);  
          }
          path_score += visited_seqs[(*p_graph_)[*itpe].sid_];
        }
        expanded_seqs.push_back(std::make_pair(concat_seq, path_score));
      } else  {
        // if the length is shorter than the expected length, push the current
        // information to the stack
        path_edges.push(next_path);
        sum_length.push(next_len);
      }
      ++ itoe; ++ num_extensions;
    }
    if(num_extensions <= 0 && current_path.size() > 1)  {
      string concat_seq = "";
      int path_score = 0;
      // record that the vertices were stored from N-terminal to C-terminal
      // so traverse it in forward order
      for(auto itpe = current_path.begin(); itpe != current_path.end(); ++ itpe) {
        string s = seqs[(*p_graph_)[*itpe].sid_];
        BoostSTRVertex v = target(*itpe, *p_graph_);
        if(concat_seq == "")  { concat_seq = s; } 
        else  { 
          int post_len = (*p_graph_)[v].len_; 
          concat_seq += s.substr(post_len, s.length() - post_len);    
        }
        path_score += visited_seqs[(*p_graph_)[*itpe].sid_];
      }
      expanded_seqs.push_back(std::make_pair(concat_seq, path_score));
    }
  }
  return;
}
  
void StringGraph::MergeExtensions(
    std::string &source_seq, 
    std::vector<std::pair<std::string, int> > &left_seqs, 
    std::vector<std::pair<std::string, int> > &right_seqs,
    std::vector<std::string> &merged_seqs
) {
  int i, j;
  // if both left and right extensions have size 0, then only take the source sequence
  if(left_seqs.size() == 0 && right_seqs.size() == 0)  {
    merged_seqs.push_back(source_seq);
  } else if(left_seqs.size() == 0) {
    // if right extension is empty, take all right extensions
    for(j = 0; j < right_seqs.size(); ++ j) {
      merged_seqs.push_back(right_seqs[j].first);
    }
  } else if(right_seqs.size() == 0) {
    // if left extension is empty, take all left extensions
    for(i = 0; i < left_seqs.size(); ++ i) {
      merged_seqs.push_back(left_seqs[i].first);
    }
  } else  {
    // if both left and right extensions are non-empty, try all combinations
    for(i = 0; i < left_seqs.size(); ++ i) {  
      merged_seqs.push_back(
        left_seqs[i].first.substr(0, left_seqs[i].first.length() - source_seq.length()) + right_seqs[0].first
      );
    }
    for(j = 0; j < right_seqs.size(); ++ j) {
      merged_seqs.push_back(
        left_seqs[0].first.substr(0, left_seqs[0].first.length() - source_seq.length()) + right_seqs[j].first
      );
    }
  }
  //*****************************************
  //cout << "Merge phase: **************************" << endl;
  //for(auto it = left_seqs.begin(); it != left_seqs.end(); ++ it) {
  //  cout << ">left: " << it->second << endl;
  //}
  //for(auto it = right_seqs.begin(); it != right_seqs.end(); ++ it) {
  //  cout << ">right: " << it->second << endl;
  //}
  return;
}

