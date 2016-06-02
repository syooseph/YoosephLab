#include "sequence_build.h"
#include "database_index.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>
#include <list>

#ifndef _UNITIGER_
#define _UNITIGER_

class Unitiger  {
 public:
  Unitiger() {return;}
  ~Unitiger() {return;}
  void GetUnitigs(
      SequenceBuild& seq_obj, 
      std::unordered_map<RIDType, std::list<OverlapType> >& fw_ext_read,
      std::unordered_map<RIDType, std::list<OverlapType> >& re_ext_read,
      std::list<std::string>& unitigs
  );
};

#endif
