#include "unitiger.h"

using namespace std;

void Unitiger::GetUnitigs(
    SequenceBuild& seq_obj, 
    std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
    std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read,
    std::list<std::string>& unitigs
) {
  /*
  for(auto it = re_ext_read.begin(); it != re_ext_read.end(); ++ it) {
    cout << it->first << endl;
  }
  return;
  */
  unordered_map<RIDType, bool> traversed_reads;
  for(auto it = fw_ext_read.begin(); it != fw_ext_read.end(); ++ it) {
    if(traversed_reads.find(it->first) != traversed_reads.end())  {
      continue;
    }
    traversed_reads[it->first] = true;
    RIDType seed_read = it->first;
    RIDType current_read = it->first;
    string unitig_seq = seq_obj.sequence_[(unsigned int) it->first];
    //cout << "read id: " << (unsigned int) current_read << endl;
    //cout << unitig_seq << endl; 
    // right extension
    bool connected = false;
    while(fw_ext_read.find(current_read) != fw_ext_read.end()) {
      auto it_l = fw_ext_read.find(current_read);
      OverlapType ol = it_l->second.front();
      if(it_l->second.size() != 1 || traversed_reads.find(ol.rid) != traversed_reads.end())  {
        traversed_reads[current_read] = true;
        break;
      } 
      // extend the sequence
      //cout << unitig_seq << endl;
      string next_seq = seq_obj.sequence_[ol.rid];
      //cout << next_seq << endl;
      unitig_seq += next_seq.substr(ol.len, next_seq.length() - ol.len);
      connected = true;
      // prepare next loop
      current_read = ol.rid;
      traversed_reads[current_read] = true;
      //cout << "end of right loop" << endl;
    }
    // left extension
    current_read = seed_read;
    while(re_ext_read.find(current_read) != re_ext_read.end()) {
      auto it_l = re_ext_read.find(current_read);
      OverlapType ol = it_l->second.front();
      if(it_l->second.size() != 1 || traversed_reads.find(ol.rid) != traversed_reads.end())  {
        traversed_reads[current_read] = true;
        break;
      } 
      // extend the sequence
      //cout << unitig_seq << endl;
      string next_seq = seq_obj.sequence_[ol.rid];
      //cout << next_seq << endl;
      unitig_seq = next_seq.substr(0, next_seq.length() - ol.len) + unitig_seq;
      connected = true;
      // prepare next loop
      current_read = ol.rid;
      traversed_reads[current_read] = true;
      //cout << "end of left loop" << endl;
    }
    if(connected)  {
      //cout << "pushed:  " << unitig_seq << endl;
      unitigs.push_back(unitig_seq);
    }
    //cout << "end of big loop" << endl;
  }
  return;
}
