#include "scoring_function.h"
#include "seq_align.h"
#include <iostream>
#include <string>

using namespace std;

typedef int AlignScoreType;

int main()  {
  string seqA = "RVIWAGIENDEIIREMA";
  string seqB = "QVALAGVDHKEIIDELT";
  ScoringFunction<AlignScoreType> score_scheme(PROTEIN, BLOSUM62, -1, -11);
  SeqAlign<AlignScoreType> align_init(seqA, seqB, &score_scheme, GLOBAL, 5, 5, false);
  align_init.Align();
  align_init.TraceBack();
  align_init.PrintAlignment();
  double bit_score = score_scheme.ComputeBitScore(align_init.GetAlignmentScore());
  cout << "bit score: " << bit_score << endl;
  return 0;
}
