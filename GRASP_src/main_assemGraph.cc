#include <iostream>
#include <vector>
#include <string>

#include "assembly_graph.h"

using namespace std;

int main()  {
  AssemblyGraph testG;
  VertexProperty v1, v2;
  v1.foo = 2;
  v2.foo = 3;
  BoostVertex u = testG.AddVertex(v1);
  BoostVertex v = testG.AddVertex(v2);
  
  EdgeProperty e;
  e.foo = 5;
  BoostEdge ed = testG.AddEdge(u, v, e);
  
  VertexProperty vcheck = testG.AccessVertex(v);
  EdgeProperty echeck = testG.AccessEdge(ed);
  
  testG.TraverseDFS();
  
  return 0;
}
