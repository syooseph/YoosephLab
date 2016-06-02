#include "seq_align.h"
#include "scoring_function.h"

#include <iostream>
#include <string>

typedef int AlignScoreType;

using namespace std;

int main()  {
  string seqA = "MKRSHYFIAVPLTSEAKQAISRFSGDASSSLPFRTWVHEEDYHITLAFLGDVPPGKMAPLCEAMAAVAAKSAPFSLALAGLGTFGERTAPRIFWQGVKEEAALNELRRDVYEACLSLGFSLDRRPFAPHITIARKWQGEAPFQPEALRSLPAASTVFSVPEIVLYRTNMEKTPKYETIAAFPLLGAPDGRTGEGMGQLLKLRDYISRYETDVYHYVPEFIRLKQW";
  string seqB = "MRAFIAIDVNESVRDSLVRAQDYIGSKEAKIKFVERENLHITLKFLGEITEEQAEEIKNILKKIAEKYKKHEVKVKGIGVFPNPNYIRVIWAGIENDEIIREMAREIEDELAKLGFKKEGNFVAHITLGRVKFVKDKLGLTMKLKELANEDFGSFVVDAIELKKSTLTPKGPIYETLARFELSE";
  cout << seqA << endl;
  cout << seqB << endl;
  ScoringFunction<int> score_scheme(PROTEIN, BLOSUM62, -1, -11);
  SeqAlign<AlignScoreType> map_reads(seqA, seqB, &score_scheme, SEMIGLOBAL, 40, 40, false);
  map_reads.Align();
  cout << "after Align, before TraceBack" << endl;
  map_reads.TraceBack();
  vector<int> aligned_region = map_reads.GetAlignmentRegion();
  cout << "aligned_region:  " << aligned_region[0] << " " << aligned_region[1] << " " << aligned_region[2] << " " << aligned_region[3] << endl;
  map_reads.PrintAlignment();
  
  //unsigned int i, j;
  //for(i = 0; i <= seqA.length(); ++ i) {
  //  for(j = 0; j <= seqB.length(); ++ j) {
  //    cout << (unsigned int) map_reads.AccessDPScoreMatrix(i, j, TRACE) << "\t";
  //  }
  //  cout << endl;
  //}
  return 0;
}
