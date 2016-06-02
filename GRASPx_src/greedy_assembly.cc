#include "greedy_assembly.h"

using namespace std;

GreedyAssembly::GreedyAssembly(void)  {
  return;
}

GreedyAssembly::~GreedyAssembly(void) {
  return;
}

void GreedyAssembly::BuildRIDLink(
    std::list<ReadType>& candidate_reads, 
    std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link
) {
  for(auto it = candidate_reads.begin(); it != candidate_reads.end(); ++ it) {
    rid_link[it->rid].push_back(it);
  }
  return;
}

bool GreedyAssembly::IsPreRegionOverlap(
    int r1_begin, int r1_end, int r2_begin, int r2_end
) {
  if(r2_begin < r1_begin && r1_begin < r2_end && r2_end < r1_end)  {
    return true;
  }
  return false;
}

bool GreedyAssembly::IsAftRegionOverlap(
    int r1_begin, int r1_end, int r2_begin, int r2_end
) {
  if(r1_begin < r2_begin && r2_begin < r1_end && r1_end < r2_end)  {
    return true;
  }
  return false;
}

ReadConnectType GreedyAssembly::GetBestConnectReadFW(
    ReachableReads& index_obj, ReadConnectType& ext_read, 
    std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link
) {
  double opt_evalue = (double) _MAX_VALUE;
  ReadConnectType opt_read; 
  opt_read.connect_terminate = true;
  list<OverlapType> possible_ext = index_obj.fw_read_ext_[ext_read.rid];
  for(auto it = possible_ext.begin(); it != possible_ext.end(); ++ it) {
    // check if the read is recruited for the query
    if(rid_link.find(it->rid) == rid_link.end())  {
      continue;
    }
    //cout << "^^^^^^:  " << it->rid << " checked" << endl;
    // if yes, search for alignment score
    for(auto it_r = rid_link[it->rid].begin(); it_r != rid_link[it->rid].end(); ++ it_r) {
      ReadListIterType opt_it = *it_r;
      //cout << " ^^^^: evalue: " << opt_it->aln_e_value << " " << opt_evalue << endl;
      //cout << " ^^^^: positions: " << ext_read.q_begin << " " << ext_read.q_end << "  " << opt_it->q_begin << " "  << opt_it->q_end << endl;
      //cout << " ^^^^: check: " << (bool) IsAftRegionOverlap(ext_read.q_begin, ext_read.q_end, opt_it->q_begin, opt_it->q_end) << endl;
      if(opt_it->aln_e_value < opt_evalue && 
          IsAftRegionOverlap(ext_read.q_begin, ext_read.q_end, opt_it->q_begin, opt_it->q_end)
      )  {
        //cout << " ^^^^: info recorded" << endl; 
        opt_evalue = opt_it->aln_e_value;
        // record read information
        opt_read.rid = opt_it->rid;
        opt_read.aln_e_value = opt_it->aln_e_value;
        opt_read.q_begin = opt_it->q_begin;
        opt_read.q_end = opt_it->q_end;
        // fill out connection-specific information
        opt_read.connect_terminate = false;
        opt_read.overlap_len_prev = ext_read.overlap_len_aft = it->len;
        opt_read.overlap_len_aft = 0;
      }
    }
  }
  return opt_read;
}

ReadConnectType GreedyAssembly::GetBestConnectReadRE(
    ReachableReads& index_obj, ReadConnectType& ext_read, 
    std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link
) {
  double opt_evalue = (double) _MAX_VALUE;
  ReadConnectType opt_read; 
  opt_read.connect_terminate = true;
  list<OverlapType> possible_ext = index_obj.re_read_ext_[ext_read.rid];
  for(auto it = possible_ext.begin(); it != possible_ext.end(); ++ it) {
    // check if the read is recruited for the query
    if(rid_link.find(it->rid) == rid_link.end())  {
      continue;
    }
    //cout << "^^^^^^:  " << it->rid << " checked" << endl;
    // if yes, search for alignment score
    for(auto it_r = rid_link[it->rid].begin(); it_r != rid_link[it->rid].end(); ++ it_r) {
      ReadListIterType opt_it = *it_r;
      //cout << " ^^^^: evalue: " << opt_it->aln_e_value << endl;
      //cout << " ^^^^: positions: " << ext_read.q_begin << " " << ext_read.q_end << "  " << opt_it->q_begin << " "  << opt_it->q_end << endl;
      if(opt_it->aln_e_value < opt_evalue &&  
           IsPreRegionOverlap(ext_read.q_begin, ext_read.q_end, opt_it->q_begin, opt_it->q_end)
      ) {
        opt_evalue = opt_it->aln_e_value;
        // record read information
        opt_read.rid = opt_it->rid;
        opt_read.aln_e_value = opt_it->aln_e_value;
        opt_read.q_begin = opt_it->q_begin;
        opt_read.q_end = opt_it->q_end;
        // fill out connection-specific information
        opt_read.connect_terminate = false;
        opt_read.overlap_len_aft = ext_read.overlap_len_prev = it->len;
        opt_read.overlap_len_prev = 0;
      }
    }
  }
  return opt_read;
}

void GreedyAssembly::GreedyExtend(
    ReachableReads& index_obj, 
    std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link,
    ReadType seed_read, double evalue_cutoff,
    std::deque<ReadConnectType>& contig_path
) {
  // initilize the contig using the seed read
  ReadConnectType rc;
  rc.rid = seed_read.rid;
  rc.q_begin = seed_read.q_begin, rc.q_end = seed_read.q_end;
  rc.aln_e_value = seed_read.aln_e_value;
  rc.overlap_len_prev = rc.overlap_len_aft = 0;
  contig_path.push_back(rc);
  // extend to both direction
  FWAssemble(index_obj, rid_link, evalue_cutoff, contig_path);
  REAssemble(index_obj, rid_link, evalue_cutoff, contig_path);
  return;
}

void GreedyAssembly::FWAssemble(
    ReachableReads& index_obj, 
    std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link, 
    double evalue_cutoff, std::deque<ReadConnectType>& contig_path
) {
  int num_low_ext = 0;
  while(num_low_ext < 3) {
    // extend
    ReadConnectType ext_read = GetBestConnectReadFW(
      index_obj, contig_path.back(), rid_link
    );
    // check goodness of the extension
    if(ext_read.connect_terminate)  {
      //cout << "extension terminate" << endl;
      break;
    } else  {
      //cout << "path extended: " << ext_read.rid << endl;
      if(ext_read.aln_e_value > evalue_cutoff)  {
        ++ num_low_ext;
      } else  {
        num_low_ext = 0;
      }
      contig_path.push_back(ext_read);
    }
  }
  // rollback to the prev good read
  while(num_low_ext > 0) {
    //cout << "path poped" << endl;
    contig_path.pop_back();
    -- num_low_ext;
  }
  return;
}

void GreedyAssembly::REAssemble(
    ReachableReads& index_obj, 
    std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link, 
    double evalue_cutoff, std::deque<ReadConnectType>& contig_path
) {
  int num_low_ext = 0;
  while(num_low_ext < 3) {
    // extend
    ReadConnectType ext_read = GetBestConnectReadRE(
      index_obj, contig_path.front(), rid_link
    );
    // check goodness of the extension
    if(ext_read.connect_terminate)  {
      break;
    } else  {
      if(ext_read.aln_e_value > evalue_cutoff)  {
        ++ num_low_ext;
      } else  {
        num_low_ext = 0;
      }
      contig_path.push_front(ext_read);
    }
  }
  // rollback to the prev good read
  while(num_low_ext > 0) {
    contig_path.pop_front();
    -- num_low_ext;
  }
  return;
}

std::string GreedyAssembly::SpellContigSequence(
    SequenceBuild& seq_obj, std::deque<ReadConnectType>& contig_path
) {
  string spell_seq = "";
  int overlap = 0;
  for(auto it = contig_path.begin(); it != contig_path.end(); ++ it) {
    spell_seq += string(seq_obj.sequence_[it->rid] + overlap);
    overlap = it->overlap_len_aft;
  }
  return spell_seq;
}
