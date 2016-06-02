#include "unitiger.h"
#include "sequence_build.h"
#include "database_index.h"

#include <iostream>
#include <string>
#include <unordered_map>
#include <list>

using namespace std;

int main(int argc, char **argv)  {
  string seq_file = argv[1];
  string index_file = argv[2];
  // load in sequence
  SequenceBuild seq_obj(seq_file);
  //cout << "Load sequence finished" << endl;
  // load in read extension index
  unordered_map<RIDType, list<OverlapType> > fw_ext, re_ext;
  DatabaseIndex db_obj;
  db_obj.LoadReadExt(index_file, fw_ext, re_ext);
  //cout << "Load index finished" << endl;
  // retrieve unitigs
  Unitiger ut_obj;
  list<string> unitigs;
  ut_obj.GetUnitigs(seq_obj, fw_ext, re_ext, unitigs);
  // print out
  int n = 0;
  for(auto it = unitigs.begin(); it != unitigs.end(); ++ it) {
    cout << ">contig" << n << endl << *it << endl;
    ++ n;
  }
  return 0;
}
